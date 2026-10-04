use super::*;
use packed_seq::{PackedSeqVec, SeqVec};
use rand::{RngExt, SeedableRng, rngs::StdRng};
use std::collections::{HashMap, HashSet};
use std::io::Cursor;

#[test]
fn split_contigs_on_ambiguous_bases() {
    let mut contig = Vec::from(b"NacGTnRYACGTN");
    assert_eq!(
        split_contig(&mut contig),
        [b"ACGT".to_vec(), b"ACGT".to_vec()]
    );
    let mut contig = Vec::from(b"NNRY");
    assert!(split_contig(&mut contig).is_empty());
}

#[test]
fn processing_does_not_join_across_ambiguous_bases() {
    let reader = MemoryReader::new(vec![b"ACGNNttgc".to_vec()]);
    let writer = Mutex::new(Vec::new());
    let stats = process(&reader, &writer, 3, 8, Some(1), false, false, 8);
    assert_eq!((stats.input_records, stats.input_bp), (2, 7));
    assert_eq!((stats.output_contigs, stats.output_bp), (2, 7));

    let mut output_reader =
        needletail::parse_fastx_reader(Cursor::new(writer.into_inner().unwrap())).unwrap();
    let mut output = Vec::new();
    while let Some(record) = output_reader.next() {
        output.push(record.unwrap().seq().into_owned());
    }
    output.sort();
    assert_eq!(output, [b"ACG".to_vec(), b"TTGC".to_vec()]);
}

struct MemoryReader {
    sequences: Mutex<Vec<Vec<u8>>>,
    idx: AtomicUsize,
}

impl MemoryReader {
    fn new(sequences: Vec<Vec<u8>>) -> Self {
        Self {
            sequences: Mutex::new(sequences),
            idx: AtomicUsize::new(0),
        }
    }
}

impl InputReader for MemoryReader {
    fn next_sample(&self) -> Option<(usize, Box<dyn Iterator<Item = Vec<u8>> + '_>)> {
        let sequence = self.sequences.lock().unwrap().pop()?;
        let idx = self.idx.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        Some((idx, Box::new(std::iter::once(sequence))))
    }
}

fn mutate(base: &[u8], rate: f64, rng: &mut impl RngExt) -> Vec<u8> {
    let n_changes = (base.len() as f64 * rate).ceil() as usize;
    let mut positions = FxHashSet::default();
    while positions.len() < n_changes {
        positions.insert(rng.random_range(0..base.len()));
    }

    let mut sequence = base.to_vec();
    for position in positions {
        let current = b"ACGT"
            .iter()
            .position(|&x| x == sequence[position])
            .unwrap();
        let replacement = rng.random_range(0..3);
        sequence[position] = b"ACGT"[(current + replacement + 1) % 4];
    }
    sequence
}

fn encode_kmer(kmer: &[u8]) -> u128 {
    kmer.iter().fold(0, |encoded, &base| {
        (encoded << 2)
            | match base {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => panic!("unexpected base {base}"),
            }
    })
}

fn encode_reverse_complement(kmer: &[u8]) -> u128 {
    kmer.iter().rev().fold(0, |encoded, &base| {
        (encoded << 2)
            | match base {
                b'A' => 3,
                b'C' => 2,
                b'G' => 1,
                b'T' => 0,
                _ => panic!("unexpected base {base}"),
            }
    })
}

fn encode_kmer_for_spectrum(kmer: &[u8], canonical: bool) -> u128 {
    let value = encode_kmer(kmer);
    if canonical {
        value.min(encode_reverse_complement(kmer))
    } else {
        value
    }
}

fn kmer_values<'a>(
    sequences: impl IntoIterator<Item = &'a [u8]>,
    k: usize,
    canonical: bool,
) -> HashSet<u128> {
    sequences
        .into_iter()
        .flat_map(|seq| {
            seq.windows(k)
                .map(move |kmer| encode_kmer_for_spectrum(kmer, canonical))
        })
        .collect()
}

fn decode_kmer(mut value: u128, k: usize) -> Vec<u8> {
    let mut kmer = vec![b'A'; k];
    for base in kmer.iter_mut().rev() {
        *base = b"ACGT"[(value & 3) as usize];
        value >>= 2;
    }
    kmer
}

fn locate_kmer(
    sequences: &[Vec<u8>],
    k: usize,
    value: u128,
    canonical: bool,
) -> Option<(usize, usize)> {
    sequences.iter().enumerate().find_map(|(sequence, bases)| {
        bases
            .windows(k)
            .position(|kmer| encode_kmer_for_spectrum(kmer, canonical) == value)
            .map(|offset| (sequence, offset))
    })
}

fn ascii_kmer(value: u128, k: usize) -> String {
    String::from_utf8(decode_kmer(value, k)).unwrap()
}

fn report_spectrum_failure(
    input: &[Vec<u8>],
    output: &[Vec<u8>],
    k: usize,
    w: usize,
    canonical: bool,
    expected: &HashSet<u128>,
    actual: &HashSet<u128>,
) -> ! {
    eprintln!("k-mer set mismatch for k={k}, w={w}");
    eprintln!("input sequences:");
    for (index, sequence) in input.iter().enumerate() {
        eprintln!("  {index}: {}", String::from_utf8_lossy(sequence));
    }
    eprintln!("output sequences:");
    for (index, sequence) in output.iter().enumerate() {
        eprintln!("  {index}: {}", String::from_utf8_lossy(sequence));
    }
    let mut input_spectrum: Vec<_> = expected.iter().map(|&value| ascii_kmer(value, k)).collect();
    input_spectrum.sort();
    let mut output_spectrum: Vec<_> = actual.iter().map(|&value| ascii_kmer(value, k)).collect();
    output_spectrum.sort();
    eprintln!("input spectrum: {input_spectrum:?}");
    eprintln!("output spectrum: {output_spectrum:?}");
    eprintln!("missing k-mers:");
    for &value in expected.difference(actual) {
        eprintln!(
            "  {} at input {:?}",
            ascii_kmer(value, k),
            locate_kmer(input, k, value, canonical)
        );
    }
    eprintln!("extra k-mers:");
    for &value in actual.difference(expected) {
        eprintln!(
            "  {} at output {:?}",
            ascii_kmer(value, k),
            locate_kmer(output, k, value, canonical)
        );
    }
    panic!("k-mer spectrum mismatch");
}

// Compact a directed de Bruijn graph by walking every maximal nonbranching
// path. Keep this independent of the production graph builder: the test input
// starts as sequences, not as precomputed unitigs.
fn naive_unitigs(sequences: &[Vec<u8>], k: usize) -> Vec<Vec<u8>> {
    fn vertex_id(vertices: &mut HashMap<Vec<u8>, usize>, label: &[u8]) -> usize {
        let next = vertices.len();
        *vertices.entry(label.to_vec()).or_insert(next)
    }

    let mut kmers: Vec<Vec<u8>> = sequences
        .iter()
        .flat_map(|sequence| sequence.windows(k).map(|window| window.to_vec()))
        .collect::<HashSet<_>>()
        .into_iter()
        .collect();
    kmers.sort_unstable();

    let mut vertices = HashMap::new();
    let edges: Vec<_> = kmers
        .into_iter()
        .map(|kmer| {
            let from = vertex_id(&mut vertices, &kmer[..k - 1]);
            let to = vertex_id(&mut vertices, &kmer[1..]);
            (from, to, kmer)
        })
        .collect();
    let mut incoming = vec![0usize; vertices.len()];
    let mut outgoing = vec![Vec::new(); vertices.len()];
    for (edge, &(from, to, _)) in edges.iter().enumerate() {
        outgoing[from].push(edge);
        incoming[to] += 1;
    }

    fn walk(
        first: usize,
        edges: &[(usize, usize, Vec<u8>)],
        incoming: &[usize],
        outgoing: &[Vec<usize>],
        visited: &mut [bool],
    ) -> Vec<u8> {
        let mut unitig = edges[first].2.clone();
        let mut edge = first;
        loop {
            visited[edge] = true;
            let to = edges[edge].1;
            if incoming[to] != 1 || outgoing[to].len() != 1 {
                break;
            }
            let next = outgoing[to][0];
            if visited[next] {
                break;
            }
            unitig.push(*edges[next].2.last().unwrap());
            edge = next;
        }
        unitig
    }

    let mut visited = vec![false; edges.len()];
    let mut unitigs = Vec::new();
    for vertex in 0..vertices.len() {
        if incoming[vertex] == 1 && outgoing[vertex].len() == 1 {
            continue;
        }
        for &edge in &outgoing[vertex] {
            if !visited[edge] {
                unitigs.push(walk(edge, &edges, &incoming, &outgoing, &mut visited));
            }
        }
    }
    // The remaining edges are isolated 1-in/1-out cycles.
    for edge in 0..edges.len() {
        if !visited[edge] {
            unitigs.push(walk(edge, &edges, &incoming, &outgoing, &mut visited));
        }
    }
    assert!(visited.into_iter().all(|done| done));
    unitigs.sort_unstable();
    unitigs
}

#[test]
fn matchtigs_on_mutated_copies_compacted_into_unitigs() {
    let k = 15;
    let mut rng = StdRng::seed_from_u64(0x9b39_7c64_602e_8d41);
    let parent: Vec<u8> = (0..350).map(|_| b"ACGT"[rng.random_range(0..4)]).collect();
    let mut copies = Vec::with_capacity(1_000);
    copies.push(parent.clone());
    for _ in 1..1_000 {
        let mut copy = parent.clone();
        let position = rng.random_range(0..copy.len());
        let original = b"ACGT"
            .iter()
            .position(|&base| base == copy[position])
            .unwrap();
        copy[position] = b"ACGT"[(original + 1 + rng.random_range(0..3)) % 4];
        copies.push(copy);
    }

    let unitigs = naive_unitigs(&copies, k);
    let input_bases: usize = unitigs.iter().map(Vec::len).sum();
    eprintln!(
        "{} copies compacted into {} input unitigs, {} total bases",
        copies.len(),
        unitigs.len(),
        input_bases,
    );
    assert!((800..=1_200).contains(&unitigs.len()));
    let mut expected: Vec<_> = unitigs
        .iter()
        .flat_map(|unitig| {
            unitig
                .windows(k)
                .map(|window| encode_kmer_for_spectrum(window, true))
        })
        .collect();
    expected.sort_unstable();
    expected.dedup();
    let copy_spectrum = kmer_values(copies.iter().map(Vec::as_slice), k, true);
    assert_eq!(
        expected.len(),
        copy_spectrum.len(),
        "naive compaction changed the copy spectrum"
    );
    assert!(expected.iter().all(|kmer| copy_spectrum.contains(kmer)));

    let mut packed = PackedSeqVec::default();
    let ranges: Vec<_> = unitigs
        .iter()
        .map(|unitig| packed.push_ascii(unitig))
        .collect();
    // Exercise both packed key implementations through graph recovery,
    // matching, and reconstruction on the same independently built unitigs.
    for (key_type, output) in [
        (
            "u64",
            matchtigs::masked_superstring::<u64>(k, packed.clone(), ranges.clone()),
        ),
        (
            "u128",
            matchtigs::masked_superstring::<u128>(k, packed, ranges),
        ),
    ] {
        eprintln!(
            "{key_type} output: 1 masked superstring, {} total bases",
            output.len()
        );
        let mut actual: Vec<_> = output
            .windows(k)
            .filter(|window| window[k - 1].is_ascii_uppercase())
            .map(|window| {
                let uppercase: Vec<_> = window.iter().map(u8::to_ascii_uppercase).collect();
                encode_kmer_for_spectrum(&uppercase, true)
            })
            .collect();
        actual.sort_unstable();
        actual.dedup();
        let missing: Vec<_> = expected
            .iter()
            .filter(|kmer| actual.binary_search(kmer).is_err())
            .copied()
            .collect();
        let extra: Vec<_> = actual
            .iter()
            .filter(|kmer| expected.binary_search(kmer).is_err())
            .copied()
            .collect();
        eprintln!(
            "{key_type} canonical k-mers: {} input, {} output, {} missing, {} extra",
            expected.len(),
            actual.len(),
            missing.len(),
            extra.len(),
        );
        assert!(
            missing.is_empty() && extra.is_empty(),
            "{key_type} masked matchtigs spectrum differs: missing {:?}, extra {:?}",
            &missing[..missing.len().min(5)],
            &extra[..extra.len().min(5)],
        );
    }
}

#[test]
fn random_sequences_preserve_the_kmer_set() {
    tracing_subscriber::fmt::init();
    let mut rng = StdRng::seed_from_u64(0x4d595df4d0f33173);

    for &n in &[10, 30, 100, 300, 1_000, 3_000] {
        eprintln!("n={n}");
        for &rate in &[0.1, 0.01, 0.001] {
            for &c in &[1, 2, 5, 10, 30, 100, 300] {
                let base: Vec<_> = (0..n).map(|_| b"ACGT"[rng.random_range(0..4)]).collect();
                let mut sequences = Vec::new();
                for _ in 0..c {
                    sequences.push(mutate(&base, rate, &mut rng));
                }

                for &k in &[3, 7, 15, 31, 63] {
                    for canonical in [false, true] {
                        let expected =
                            kmer_values(sequences.iter().map(Vec::as_slice), k, canonical);
                        for &w in &[8, 24, 100, 500] {
                            for reference in [false, true] {
                                for threads in [1, 3] {
                                    let reader = MemoryReader::new(sequences.clone());
                                    let writer = Mutex::new(Vec::new());
                                    process(
                                        &reader,
                                        &writer,
                                        k,
                                        w,
                                        Some(threads),
                                        reference,
                                        canonical,
                                        8,
                                    );

                                    let output = writer.into_inner().unwrap();
                                    let mut output_reader =
                                        needletail::parse_fastx_reader(Cursor::new(output))
                                            .unwrap();
                                    let mut output_sequences = Vec::new();
                                    while let Some(record) = output_reader.next() {
                                        output_sequences.push(record.unwrap().seq().into_owned());
                                    }
                                    let actual = kmer_values(
                                        output_sequences.iter().map(Vec::as_slice),
                                        k,
                                        canonical,
                                    );
                                    if actual != expected {
                                        report_spectrum_failure(
                                            &sequences,
                                            &output_sequences,
                                            k,
                                            w,
                                            canonical,
                                            &expected,
                                            &actual,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
