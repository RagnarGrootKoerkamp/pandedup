use super::*;
use rand::{rngs::StdRng, RngExt, SeedableRng};
use std::collections::HashSet;
use std::io::Cursor;

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

fn kmer_values<'a>(sequences: impl IntoIterator<Item = &'a [u8]>, k: usize) -> HashSet<u128> {
    sequences
        .into_iter()
        .flat_map(|seq| seq.windows(k).map(encode_kmer))
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

fn locate_kmer(sequences: &[Vec<u8>], k: usize, value: u128) -> Option<(usize, usize)> {
    sequences.iter().enumerate().find_map(|(sequence, bases)| {
        bases
            .windows(k)
            .position(|kmer| encode_kmer(kmer) == value)
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
            locate_kmer(input, k, value)
        );
    }
    eprintln!("extra k-mers:");
    for &value in actual.difference(expected) {
        eprintln!(
            "  {} at output {:?}",
            ascii_kmer(value, k),
            locate_kmer(output, k, value)
        );
    }
    panic!("k-mer spectrum mismatch");
}

#[test]
fn random_sequences_preserve_the_kmer_set() {
    let mut rng = StdRng::seed_from_u64(0x4d595df4d0f33173);

    for &n in &[10, 30, 100, 1_000, 10_000] {
        eprintln!("n={n}");
        for &rate in &[0.1, 0.01, 0.001] {
            for &c in &[1, 2, 5, 10, 100, 1000] {
                let base: Vec<_> = (0..n).map(|_| b"ACGT"[rng.random_range(0..4)]).collect();
                let mut sequences = Vec::new();
                for _ in 0..c {
                    sequences.push(mutate(&base, rate, &mut rng));
                }

                for &k in &[3, 5, 7, 15, 31, 63, 64] {
                    let expected = kmer_values(sequences.iter().map(Vec::as_slice), k);
                    for &w in &[5, 10, 25, 50, 100, 200, 500] {
                        let args = Args {
                            input: PathBuf::new(),
                            output: None,
                            k,
                            w,
                            threads: Some(1),
                            reference: false,
                            canonical: false,
                            mini_k: 8,
                        };
                        let reader = MemoryReader::new(sequences.clone());
                        let seen: &[_; 256] =
                            &std::array::from_fn(|_i| RwLock::new(FxHashSet::default()));
                        let global_stats = Mutex::new(Stats::default());
                        let writer = Mutex::new(Vec::new());
                        let reference = RwLock::new((Vec::new(), FxHashMap::default()));

                        while process_sample(
                            &args,
                            &reader,
                            seen,
                            &global_stats,
                            &writer,
                            &reference,
                        )
                        .is_some()
                        {}

                        let output = writer.into_inner().unwrap();
                        let mut output_reader =
                            needletail::parse_fastx_reader(Cursor::new(output)).unwrap();
                        let mut output_sequences = Vec::new();
                        while let Some(record) = output_reader.next() {
                            output_sequences.push(record.unwrap().seq().into_owned());
                        }
                        let actual = kmer_values(output_sequences.iter().map(Vec::as_slice), k);
                        if actual != expected {
                            report_spectrum_failure(
                                &sequences,
                                &output_sequences,
                                k,
                                w,
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
