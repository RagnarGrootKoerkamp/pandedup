//! Greedy matching of unitig ends through the compacted de Bruijn graph.

use crate::{log_file_stats, mss::MssKey, timing::StageTiming};
use packed_seq::{PackedSeqVec, SeqVec, complement_char};
use std::io::{BufWriter, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use tracing::info;

const DEAD: u32 = u32::MAX;

#[derive(Clone, Copy)]
#[repr(C, packed(1))]
struct End<K> {
    kmer: K,
    id: u32,
}

const _: () = {
    assert!(std::mem::size_of::<End<u64>>() == 12);
    assert!(std::mem::size_of::<End<u128>>() == 20);
};

/// CSR adjacency indexed by an oriented outgoing tail. An edge's destination
/// is the oriented unitig entered at its head; its low bit is the orientation.
struct Graph {
    edges: Vec<u32>,
    offsets: Vec<u32>,
}

impl Graph {
    fn outgoing(&self, id: u32) -> &[u32] {
        &self.edges[self.offsets[id as usize] as usize..self.offsets[id as usize + 1] as usize]
    }
}

fn graph<K: MssKey>(k: usize, seq: &PackedSeqVec, ranges: &[Range<usize>]) -> Graph {
    let overlap = k - 1;
    let mut heads = Vec::with_capacity(2 * ranges.len());
    let mut tails = Vec::with_capacity(2 * ranges.len());
    info!("Building graph of {} unitigs..", ranges.len());
    for (index, range) in ranges.iter().enumerate() {
        let id = (index as u32) * 2;
        let (head, tail, rc_head, rc_tail): (K, K, K, K) = (
            K::read_kmer(seq, overlap, range.start),
            K::read_kmer(seq, overlap, range.end - overlap),
            K::read_revcomp_kmer(seq, overlap, range.start),
            K::read_revcomp_kmer(seq, overlap, range.end - overlap),
        );
        heads.push(End { kmer: head, id });
        heads.push(End {
            kmer: rc_tail,
            id: id + 1,
        });
        tails.push(End { kmer: tail, id });
        tails.push(End {
            kmer: rc_head,
            id: id + 1,
        });
    }
    info!("Sorting heads");
    heads.sort_unstable_by_key(|end| end.kmer);
    info!("Sorting tails");
    tails.sort_unstable_by_key(|end| end.kmer);

    info!("Matching unitig ends");
    let mut pairs = Vec::new();
    let (mut h, mut t) = (0, 0);
    while h < heads.len() && t < tails.len() {
        // Copy keys out of packed entries before borrowing them for comparison.
        let head_key = heads[h].kmer;
        let tail_key = tails[t].kmer;
        match head_key.cmp(&tail_key) {
            std::cmp::Ordering::Less => h += 1,
            std::cmp::Ordering::Greater => t += 1,
            std::cmp::Ordering::Equal => {
                let key = heads[h].kmer;
                let mut he = h;
                let mut te = t;
                while he < heads.len() {
                    let next = heads[he].kmer;
                    if next != key {
                        break;
                    }
                    he += 1;
                }
                while te < tails.len() {
                    let next = tails[te].kmer;
                    if next != key {
                        break;
                    }
                    te += 1;
                }
                for tail in &tails[t..te] {
                    for head in &heads[h..he] {
                        // The reverse orientation of the same endpoint does not
                        // traverse a unitig or introduce a new graph edge.
                        if tail.id != (head.id ^ 1) {
                            pairs.push((tail.id, head.id));
                        }
                    }
                }
                h = he;
                t = te;
            }
        }
    }
    info!("Sorting {} edges", pairs.len());
    pairs.sort_unstable();
    // let old_len = pairs.len();
    // info!("Dedup edges");
    // pairs.dedup();
    // info!("Dedup edges => {} left", pairs.len());
    // assert_eq!(old_len, pairs.len(), "duplicate graph edges");
    assert!(pairs.len() < u32::MAX as usize, "too many graph edges");
    info!("Building CSR adjacency");
    let mut offsets = vec![0; 2 * ranges.len() + 1];
    info!(
        "Filling offsets: {} MB",
        std::mem::size_of_val(offsets.as_slice()) / (1024 * 1024)
    );
    for &(from, _) in &pairs {
        offsets[from as usize + 1] += 1;
    }
    for i in 1..offsets.len() {
        offsets[i] += offsets[i - 1];
    }
    info!("Mapping pairs");
    Graph {
        edges: pairs.into_iter().map(|(_, to)| to).collect(),
        offsets,
    }
}

#[derive(Clone, Copy)]
#[repr(C, packed(1))]
struct Link {
    target: u32,
    overlap: u8,
}

const _: () = assert!(std::mem::size_of::<Link>() == 5);

impl Link {
    const EMPTY: Self = Self {
        target: DEAD,
        overlap: 0,
    };
    fn is_empty(self) -> bool {
        self.target == DEAD
    }
}

#[derive(Clone, Copy)]
struct ActivePath {
    start: u32,
    current: u32,
}

fn tail_slot(id: u32) -> usize {
    (id as usize & !1) + usize::from(id & 1 == 0)
}
fn head_slot(id: u32) -> usize {
    (id as usize & !1) + usize::from(id & 1 != 0)
}

fn connect(links: &mut [Link], from: u32, to: u32, overlap: usize) {
    let overlap = overlap as u8;
    links[tail_slot(from)] = Link {
        target: to,
        overlap,
    };
    links[head_slot(to)] = Link {
        target: from ^ 1,
        overlap,
    };
}

/// Link every unitig end, preferring paths through the graph whose added
/// sequence is shorter than a zero-overlap join.
fn match_ends(k: usize, ranges: &[Range<usize>], graph: &Graph) -> Vec<Link> {
    let n = ranges.len();
    info!(
        "Allocating {} links: {} MB",
        2 * n,
        std::mem::size_of::<Link>() * 2 * n / (1024 * 1024)
    );
    let mut links = vec![Link::EMPTY; 2 * n];
    info!(
        "Allocating {} receiving ids: {} MB",
        2 * n,
        std::mem::size_of::<u32>() * 2 * n / (1024 * 1024)
    );
    let mut receiving: Vec<u32> = (0..2 * n as u32).collect();
    let mut buckets: Vec<Vec<ActivePath>> = vec![Vec::new(); k];
    eprintln!("Pushing to initial bucket");
    for id in 0..2 * n as u32 {
        buckets[0].push(ActivePath {
            start: id,
            current: id,
        });
    }

    // At distance k there is no shared base left, so a graph path cannot
    // improve on the zero-overlap fallback.
    for distance in 0..k {
        info!(
            "Distance {distance}: {} active paths, {} receiving ends",
            buckets[distance].len(),
            receiving.len()
        );
        let mut candidates = Vec::<ActivePath>::new();
        for path in std::mem::take(&mut buckets[distance]) {
            if !links[tail_slot(path.start)].is_empty() {
                continue;
            }
            for &next in graph.outgoing(path.current) {
                candidates.push(ActivePath {
                    start: path.start,
                    current: next,
                });
            }
        }
        // The receiving ids and candidate tail ids have the same orientation
        // encoding. Sorting permits a single merge pass over both vectors.
        candidates.sort_unstable_by_key(|path| (path.current, path.start));
        let mut r = 0;
        for candidate in &mut candidates {
            while r < receiving.len() && (receiving[r] == DEAD || receiving[r] < candidate.current)
            {
                r += 1;
            }
            if r == receiving.len() || receiving[r] != candidate.current {
                continue;
            }
            if !links[tail_slot(candidate.start)].is_empty() {
                candidate.start = DEAD;
                continue;
            }
            if tail_slot(candidate.start) == head_slot(candidate.current) {
                continue;
            }
            if links[head_slot(candidate.current)].is_empty() {
                connect(
                    &mut links,
                    candidate.start,
                    candidate.current,
                    k - 1 - distance,
                );
                receiving[r] = DEAD;
                candidate.start = DEAD;
            } else {
                receiving[r] = DEAD;
            }
        }
        receiving.retain(|&id| id != DEAD && links[head_slot(id)].is_empty());

        for candidate in candidates {
            if candidate.start == DEAD || !links[tail_slot(candidate.start)].is_empty() {
                continue;
            }
            let weight = ranges[(candidate.current / 2) as usize].len() - (k - 1);
            let next_distance = distance.saturating_add(weight);
            if next_distance >= k {
                continue;
            }
            buckets[next_distance].push(ActivePath {
                start: candidate.start,
                current: candidate.current,
            });
        }
        info!(
            "After distance {distance}, bucket queue sizes: {:?}",
            buckets.iter().map(Vec::len).collect::<Vec<_>>()
        );
    }

    // Every free outgoing end and free receiving end must be paired so that
    // reconstruction consists entirely of cycles.
    let outgoing: Vec<_> = (0..2 * n as u32)
        .filter(|&id| links[tail_slot(id)].is_empty())
        .collect();
    // Each free physical end also has a receiving orientation (id ^ 1).
    // Pair physical ends once each, then enter the second in that orientation.
    assert_eq!(outgoing.len() % 2, 0, "odd number of unmatched ends");
    for pair in outgoing.chunks_exact(2) {
        connect(&mut links, pair[0], pair[1] ^ 1, 0);
    }
    links
}

fn append(
    output: &mut Vec<u8>,
    seq: &PackedSeqVec,
    range: &Range<usize>,
    reverse: bool,
    overlap: usize,
    k: usize,
) {
    let mut bases = seq.slice(range.clone()).unpack();
    if reverse {
        bases.reverse();
        for base in &mut bases {
            *base = complement_char(*base);
        }
    }
    if !output.is_empty() {
        let mask = (k - 1 - overlap).min(bases.len() - overlap);
        for base in &mut bases[overlap..overlap + mask] {
            *base = base.to_ascii_lowercase();
        }
    }
    output.extend_from_slice(&bases[overlap..]);
}

/// Build a masked superstring from unitigs. Newly created crossing k-mers are
/// lowercase, as in `mss::masked_superstring`.
pub fn masked_superstring<K: MssKey>(
    k: usize,
    seq: PackedSeqVec,
    ranges: Vec<Range<usize>>,
) -> Vec<u8> {
    assert!(k > 0 && k <= K::BITS / 2);
    let (ranges, short): (Vec<_>, Vec<_>) = ranges.into_iter().partition(|r| r.len() >= k);
    assert!(ranges.len() <= (u32::MAX as usize / 2), "too many unitigs");
    let graph = graph::<K>(k, &seq, &ranges);

    let links = match_ends(k, &ranges, &graph);

    info!("Reconstruct output");
    let mut output = Vec::new();
    let mut done = vec![false; ranges.len()];
    for start in 0..ranges.len() {
        if done[start] {
            continue;
        }
        let mut id = (start as u32) * 2;
        let mut overlap = 0;
        loop {
            let index = (id / 2) as usize;
            assert!(!done[index], "link cycle revisits unitig before closing");
            done[index] = true;
            append(&mut output, &seq, &ranges[index], id & 1 != 0, overlap, k);
            let link = links[tail_slot(id)];
            assert!(!link.is_empty(), "unmatched unitig end");
            if (link.target / 2) as usize == start {
                break;
            }
            id = link.target;
            overlap = link.overlap as usize;
        }
    }
    for range in &short {
        if !range.is_empty() {
            append(&mut output, &seq, range, false, 0, k);
        }
    }
    output
}

/// Read unitigs and write one masked superstring in FASTA format.
pub fn run(input: &Path, output: Option<&Path>, k: usize) -> PathBuf {
    let timing = StageTiming::start();
    info!("Reading unitigs..");
    let (seq, ranges) = PackedSeqVec::from_fastx(input);
    let input_bases = ranges.iter().map(Range::len).sum();
    log_file_stats("Read", input, Some((ranges.len(), input_bases))).unwrap();
    let superstring = if k <= 32 {
        masked_superstring::<u64>(k, seq, ranges)
    } else {
        masked_superstring::<u128>(k, seq, ranges)
    };
    let output = output
        .map(Path::to_path_buf)
        .unwrap_or_else(|| input.with_extension("matchtigs.msfa"));
    let mut writer = BufWriter::new(std::fs::File::create(&output).unwrap());
    writeln!(writer, ">matchtigs-masked-superstring").unwrap();
    writer.write_all(&superstring).unwrap();
    writer.write_all(b"\n").unwrap();
    drop(writer);
    log_file_stats("Wrote", &output, Some((1, superstring.len()))).unwrap();
    info!(
        "matchtigs: {}, output {} bases ({})",
        timing.finish(),
        superstring.len(),
        output.display()
    );
    output
}
