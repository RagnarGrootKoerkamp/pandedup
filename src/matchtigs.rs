//! Greedy matching of unitig ends through the compacted de Bruijn graph.

use crate::{log_file_stats, mss::MssKey, timing::StageTiming};
use packed_seq::{PackedSeqVec, SeqVec, complement_char};
use std::io::{BufWriter, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use tracing::info;
use voracious_radix_sort::{RadixSort, Radixable};

const DEAD: u32 = u32::MAX;

#[derive(Clone, Copy)]
#[repr(C, packed(1))]
#[doc(hidden)]
pub struct End<K> {
    kmer: K,
    id: u32,
}

const _: () = {
    assert!(std::mem::size_of::<End<u64>>() == 12);
    assert!(std::mem::size_of::<End<u128>>() == 20);
};

impl Radixable<u64> for End<u64> {
    type Key = u64;
    fn key(&self) -> Self::Key {
        self.kmer
    }
}

impl Radixable<u128> for End<u128> {
    type Key = u128;
    fn key(&self) -> Self::Key {
        self.kmer
    }
}

impl<K: MssKey> PartialEq for End<K> {
    fn eq(&self, other: &Self) -> bool {
        let left = self.kmer;
        let right = other.kmer;
        left == right
    }
}

impl<K: MssKey> PartialOrd for End<K> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let left = self.kmer;
        let right = other.kmer;
        Some(left.cmp(&right))
    }
}

#[derive(Clone, Copy, PartialEq, PartialOrd)]
struct Edge {
    from: u32,
    to: u32,
}

impl Radixable<u64> for Edge {
    type Key = u64;
    fn key(&self) -> Self::Key {
        (u64::from(self.from) << 32) | u64::from(self.to)
    }
}

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

fn graph<K: MssKey>(k: usize, seq: &PackedSeqVec, ranges: &[Range<usize>]) -> Graph
where
    End<K>: Radixable<K, Key = K>,
{
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
    heads.voracious_mt_sort(rayon::current_num_threads());
    info!("Sorting tails");
    tails.voracious_mt_sort(rayon::current_num_threads());

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
                            pairs.push(Edge {
                                from: tail.id,
                                to: head.id,
                            });
                        }
                    }
                }
                h = he;
                t = te;
            }
        }
    }
    info!("Sorting {} edges", pairs.len());
    pairs.voracious_mt_sort(rayon::current_num_threads());
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
    for edge in &pairs {
        offsets[edge.from as usize + 1] += 1;
    }
    for i in 1..offsets.len() {
        offsets[i] += offsets[i - 1];
    }
    info!("Mapping pairs");
    Graph {
        edges: pairs.into_iter().map(|edge| edge.to).collect(),
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

const UNREACHABLE: u8 = u8::MAX;

fn traversal_weight(ranges: &[Range<usize>], node: u32, k: usize) -> usize {
    ranges[(node / 2) as usize].len() - (k - 1)
}

/// Nearest still-free outgoing end for each oriented unitig tail.
struct DistanceField {
    distance: Vec<u8>,
    owner: Vec<u32>,
    parent: Vec<u32>,
    unsettled: Vec<bool>,
    affected: Vec<u32>,
}

impl DistanceField {
    fn new(nodes: usize) -> Self {
        Self {
            distance: vec![0; nodes],
            owner: (0..nodes as u32).collect(),
            parent: vec![DEAD; nodes],
            unsettled: vec![false; nodes],
            affected: Vec::new(),
        }
    }

    fn relax(&mut self, node: u32, distance: usize, owner: u32, parent: u32, k: usize) -> bool {
        let index = node as usize;
        if distance >= k || !self.unsettled[index] || distance >= self.distance[index] as usize {
            return false;
        }
        let newly_reached = self.distance[index] == UNREACHABLE;
        self.distance[index] = distance as u8;
        self.owner[index] = owner;
        self.parent[index] = parent;
        newly_reached
    }

    /// Delete one source and repair exactly the part of its shortest-path tree
    /// that depended on it. Incoming edges use the reverse-complement graph.
    fn remove_source(&mut self, source: u32, graph: &Graph, ranges: &[Range<usize>], k: usize) {
        assert_eq!(
            self.owner[source as usize], source,
            "removed source lost its root label"
        );
        assert_eq!(
            self.distance[source as usize], 0,
            "active source has nonzero distance"
        );
        self.affected.clear();
        self.affected.push(source);
        self.unsettled[source as usize] = true;
        let mut cursor = 0;
        while cursor < self.affected.len() {
            let node = self.affected[cursor];
            cursor += 1;
            for &child in graph.outgoing(node) {
                let index = child as usize;
                if self.owner[index] == source
                    && self.parent[index] == node
                    && !self.unsettled[index]
                {
                    self.unsettled[index] = true;
                    self.affected.push(child);
                }
            }
            self.distance[node as usize] = UNREACHABLE;
            self.owner[node as usize] = DEAD;
            self.parent[node as usize] = DEAD;
        }

        let mut pending = 0usize;
        let mut first_distance = k;
        for index in 0..self.affected.len() {
            let node = self.affected[index];
            let weight = traversal_weight(ranges, node, k);
            for &reverse_predecessor in graph.outgoing(node ^ 1) {
                let predecessor = reverse_predecessor ^ 1;
                let from = predecessor as usize;
                if !self.unsettled[from] && self.owner[from] != DEAD {
                    let distance = (self.distance[from] as usize).saturating_add(weight);
                    if self.relax(node, distance, self.owner[from], predecessor, k) {
                        pending += 1;
                    }
                    first_distance = first_distance.min(distance);
                }
            }
        }

        // Every edge adds at least one base, so one scan of the affected
        // nodes per distance settles all labels at that distance. A node
        // occupies one slot in `affected`, even if many paths relax it.
        for distance in first_distance..k {
            for position in 0..self.affected.len() {
                let node = self.affected[position];
                let index = node as usize;
                if !self.unsettled[index] || self.distance[index] as usize != distance {
                    continue;
                }
                self.unsettled[index] = false;
                pending -= 1;
                let owner = self.owner[index];
                for &next in graph.outgoing(node) {
                    if self.relax(
                        next,
                        distance.saturating_add(traversal_weight(ranges, next, k)),
                        owner,
                        node,
                        k,
                    ) {
                        pending += 1;
                    }
                }
            }
            if pending == 0 {
                break;
            }
        }
        for &node in &self.affected {
            self.unsettled[node as usize] = false;
        }
    }
}

/// Reusable scratch for the rare case where a receiver's nearest owner is
/// its own physical end, which cannot be linked to itself.
struct ReverseSearch {
    distance: Vec<u8>,
    touched: Vec<u32>,
}

impl ReverseSearch {
    fn new(nodes: usize) -> Self {
        Self {
            distance: vec![UNREACHABLE; nodes],
            touched: Vec::new(),
        }
    }

    fn push(&mut self, node: u32, distance: usize, k: usize) {
        if distance >= k || distance >= self.distance[node as usize] as usize {
            return;
        }
        if self.distance[node as usize] == UNREACHABLE {
            self.touched.push(node);
        }
        self.distance[node as usize] = distance as u8;
    }

    fn nearest(
        &mut self,
        target: u32,
        limit: usize,
        field: &DistanceField,
        graph: &Graph,
        ranges: &[Range<usize>],
        k: usize,
    ) -> Option<(u32, usize)> {
        for &reverse_predecessor in graph.outgoing(target ^ 1) {
            self.push(reverse_predecessor ^ 1, 0, k);
        }
        let mut result = None;
        let mut best_distance = limit;
        for distance in 0..limit.min(k) {
            if distance >= best_distance {
                break;
            }
            let mut cursor = 0;
            while cursor < self.touched.len() {
                let node = self.touched[cursor];
                cursor += 1;
                if self.distance[node as usize] as usize != distance {
                    continue;
                }
                let source = field.owner[node as usize];
                let combined = distance + field.distance[node as usize] as usize;
                if source != DEAD && source != target ^ 1 && combined < best_distance {
                    best_distance = combined;
                    result = Some((source, combined));
                }
                let next_distance = distance.saturating_add(traversal_weight(ranges, node, k));
                if next_distance < best_distance && next_distance < k {
                    for &reverse_predecessor in graph.outgoing(node ^ 1) {
                        self.push(reverse_predecessor ^ 1, next_distance, k);
                    }
                }
            }
        }
        for &node in &self.touched {
            self.distance[node as usize] = UNREACHABLE;
        }
        self.touched.clear();
        result
    }
}

fn tail_slot(id: u32) -> usize {
    (id ^ 1) as usize
}
fn head_slot(id: u32) -> usize {
    id as usize
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

/// Find the cheapest still-free source for one receiving head. Distances in the
/// field are to outgoing tails, so entering the receiver costs nothing here.
fn nearest_receiver(
    target: u32,
    field: &DistanceField,
    reverse_search: &mut ReverseSearch,
    graph: &Graph,
    ranges: &[Range<usize>],
    k: usize,
) -> Option<(u32, usize)> {
    let forbidden = target ^ 1;
    let mut best = None;
    let mut forbidden_distance = usize::MAX;
    for &reverse_predecessor in graph.outgoing(target ^ 1) {
        let predecessor = reverse_predecessor ^ 1;
        let index = predecessor as usize;
        let source = field.owner[index];
        if source == DEAD {
            continue;
        }
        let distance = field.distance[index] as usize;
        if source == forbidden {
            forbidden_distance = forbidden_distance.min(distance);
        } else if best.is_none_or(|(_, old_distance)| distance < old_distance) {
            best = Some((source, distance));
        }
    }
    let limit = best.map_or(k, |(_, distance)| distance);
    if forbidden_distance <= limit {
        reverse_search
            .nearest(target, limit, field, graph, ranges, k)
            .or(best)
    } else {
        best
    }
}

/// Greedy distance-ordered matching with one shortest-path label per graph
/// node. At each distance, expand the outgoing edges of all nodes on that
/// frontier, match the receiving ends they reach, then repair labels after
/// deleting the matched sources. Repeat until that distance is exhausted.
fn match_ends(k: usize, ranges: &[Range<usize>], graph: &Graph) -> Vec<Link> {
    let nodes = 2 * ranges.len();
    let mut links = vec![Link::EMPTY; nodes];
    let mut field = DistanceField::new(nodes);
    let mut reverse_search = ReverseSearch::new(nodes);
    let mut seen = vec![false; nodes];
    let mut matched_sources = Vec::new();
    let mut linked = 0usize;

    for distance in 0..k {
        let mut passes = 0usize;
        loop {
            passes += 1;
            seen.fill(false);
            matched_sources.clear();
            let mut frontier = 0usize;
            for node in 0..nodes as u32 {
                if field.distance[node as usize] as usize != distance {
                    continue;
                }
                frontier += 1;
                // The final edge enters the receiver at zero cost. Charging
                // its unitig would include the last unitig in the distance.
                for &target in graph.outgoing(node) {
                    if seen[target as usize] || !links[head_slot(target)].is_empty() {
                        continue;
                    }
                    seen[target as usize] = true;
                    let Some((source, actual_distance)) =
                        nearest_receiver(target, &field, &mut reverse_search, graph, ranges, k)
                    else {
                        continue;
                    };
                    assert!(actual_distance >= distance, "receiver distance decreased");
                    if actual_distance != distance || !links[tail_slot(source)].is_empty() {
                        continue;
                    }
                    connect(&mut links, source, target, k - 1 - distance);
                    matched_sources.push(source);
                    linked += 1;
                }
            }
            info!(
                "Distance {distance}, pass {passes}: {frontier} frontier nodes, {} links",
                matched_sources.len()
            );
            if matched_sources.is_empty() {
                break;
            }
            // Both physical ends of each link are independent sources.
            // Their removal can reveal another match at this same distance.
            for &source in &matched_sources {
                let target = links[tail_slot(source)].target;
                field.remove_source(source, graph, ranges, k);
                field.remove_source(target ^ 1, graph, ranges, k);
            }
        }
        info!("After distance {distance}: {linked} links");
    }

    let outgoing: Vec<_> = (0..nodes as u32)
        .filter(|&id| links[tail_slot(id)].is_empty())
        .collect();
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
) -> Vec<u8>
where
    End<K>: Radixable<K, Key = K>,
{
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matching_revisits_a_distance_after_removing_sources() {
        // Each edge has its reverse-complement counterpart. The first match
        // consumes source 0 and its reverse end 3. Those sources initially
        // hide the direct match between source 6 and receiver 4.
        let graph = Graph {
            edges: vec![2, 4, 1, 7, 1, 7, 2, 4],
            offsets: vec![0, 2, 2, 2, 4, 4, 6, 8, 8],
        };
        let k = 5;
        let ranges = vec![0..k; 4];
        let links = match_ends(k, &ranges, &graph);
        let direct = links
            .iter()
            .filter(|link| link.overlap == (k - 1) as u8)
            .count();
        assert_eq!(direct, 4);
        let first_target = links[tail_slot(0)].target;
        let second_target = links[tail_slot(6)].target;
        assert_eq!(first_target, 2);
        assert_eq!(second_target, 4);
    }
}
