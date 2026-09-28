//! Run without arguments.
//! This investigates over-reported kmers in our deduped version of the HPRCv2
//! agc, because Deacon finds more minimizers in our deduped version that in the original.
//! Collects the smallest 0.1% of the 64-mers in hprcv2.k64.fa.zst.
//! Then, it iterates over the hprcv2.agc and removes kmers once it has seen them.
//! In the end, only over-reported kmers remain; ideally there are none of those...
use std::{
    env,
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
    sync::{
        RwLock,
        atomic::{AtomicUsize, Ordering},
    },
};

use fxhash::FxHashSet;
use ragc_core::{Decompressor, DecompressorConfig};

const K: usize = 64;
// The bottom 0.1% of the 2^128 possible 64-mer values.
const VALUE_LIMIT: u128 = u128::MAX / 1000;
const MERGE_BATCH_SIZE: usize = 16_384;

fn visit_sequence_kmers(sequence: &[u8], mut visit: impl FnMut(u128)) {
    let mut value = 0u128;
    let mut valid_bases = 0usize;
    for &base in sequence {
        let digit = match base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => {
                value = 0;
                valid_bases = 0;
                continue;
            }
        };
        // For K=64, shifting discards the previous leading base automatically.
        value = (value << 2) | digit;
        valid_bases += 1;
        if valid_bases >= K && value < VALUE_LIMIT {
            visit(value);
        }
    }
}

fn add_sequence_kmers(sequence: &[u8], kmers: &mut FxHashSet<u128>) {
    visit_sequence_kmers(sequence, |value| {
        kmers.insert(value);
    });
}

fn remove_present_kmers(batch: &mut FxHashSet<u128>, remaining: &RwLock<FxHashSet<u128>>) {
    if batch.is_empty() {
        return;
    }
    {
        let existing = remaining.read().unwrap();
        batch.retain(|value| existing.contains(value));
    }
    if !batch.is_empty() {
        let mut remaining = remaining.write().unwrap();
        for value in batch.drain() {
            remaining.remove(&value);
        }
    }
}

fn remove_agc_kmers(path: &Path, output: FxHashSet<u128>) -> FxHashSet<u128> {
    let decompressor = Decompressor::open(
        path.to_str().expect("AGC path must be valid UTF-8"),
        DecompressorConfig { verbosity: 0 },
    )
    .unwrap();
    let samples = decompressor.list_samples();
    let next_sample = AtomicUsize::new(0);
    let remaining = RwLock::new(output);
    let threads = num_cpus::get_physical().max(1).min(samples.len());
    eprintln!(
        "reading {} AGC samples with {threads} threads",
        samples.len()
    );

    // A Decompressor is not Sync; each worker gets its own archive reader.
    let readers: Vec<_> = (0..threads)
        .map(|_| decompressor.clone_for_thread().unwrap())
        .collect();
    std::thread::scope(|scope| {
        for mut reader in readers {
            let samples = &samples;
            let next_sample = &next_sample;
            let remaining = &remaining;
            scope.spawn(move || {
                let mut batch = FxHashSet::default();
                loop {
                    let i = next_sample.fetch_add(1, Ordering::Relaxed);
                    let Some(sample) = samples.get(i) else { break };
                    for contig in reader.list_contigs(sample).unwrap() {
                        let mut sequence = reader.get_contig(sample, &contig).unwrap();
                        sequence
                            .iter_mut()
                            .for_each(|base| *base = b"ACGT"[(*base as usize) % 4]);
                        visit_sequence_kmers(&sequence, |value| {
                            batch.insert(value);
                            if batch.len() >= MERGE_BATCH_SIZE {
                                remove_present_kmers(&mut batch, remaining);
                            }
                        });
                    }
                    remove_present_kmers(&mut batch, remaining);
                    let unmatched = remaining.read().unwrap().len();
                    eprintln!("finished sample {i}: {sample}; unmatched kmers: {unmatched}");
                }
            });
        }
    });
    remaining.into_inner().unwrap()
}

fn fasta_kmers(path: &Path) -> FxHashSet<u128> {
    let decoder =
        zstd::stream::read::Decoder::new(BufReader::new(File::open(path).unwrap())).unwrap();
    let mut reader = needletail::parse_fastx_reader(decoder).unwrap();
    let mut kmers = FxHashSet::default();
    while let Some(record) = reader.next() {
        add_sequence_kmers(&record.unwrap().seq(), &mut kmers);
    }
    eprintln!(
        "finished reading {}: kmers: {}",
        path.display(),
        kmers.len()
    );
    kmers
}

fn main() {
    let mut args = env::args_os().skip(1);
    let agc = PathBuf::from(args.next().unwrap_or_else(|| "hprcv2.agc".into()));
    let fasta = PathBuf::from(args.next().unwrap_or_else(|| "hprcv2.k64.fa.zst".into()));

    eprintln!("reading {}", fasta.display());
    let output = fasta_kmers(&fasta);
    let output_count = output.len();
    eprintln!("reading {}", agc.display());
    let remaining = remove_agc_kmers(&agc, output);

    println!("output kmers retained: {}", output_count);
    println!("also present in AGC:   {}", output_count - remaining.len());
    println!("extra in output:      {}", remaining.len());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rolling_kmers_respect_sequence_boundaries() {
        let mut kmers = FxHashSet::default();
        let low = b"A".repeat(K);
        add_sequence_kmers(&low, &mut kmers);
        assert_eq!(kmers.len(), 1);
        assert!(kmers.contains(&0));

        // The invalid base splits the input into two fragments shorter than K.
        let mut split = b"A".repeat(K - 1);
        split.push(b'N');
        split.extend_from_slice(&b"A".repeat(K - 1));
        let mut split_kmers = FxHashSet::default();
        add_sequence_kmers(&split, &mut split_kmers);
        assert!(split_kmers.is_empty());

        // High-valued kmers are excluded by the 1% filter.
        let mut high_kmers = FxHashSet::default();
        add_sequence_kmers(&b"T".repeat(K), &mut high_kmers);
        assert!(high_kmers.is_empty());

        let mut mixed = b"A".repeat(K);
        mixed.extend_from_slice(b"CGTAN");
        mixed.extend_from_slice(&b"A".repeat(K));
        mixed.extend_from_slice(b"CGTA");
        let expected: FxHashSet<u128> = mixed
            .windows(K)
            .filter_map(|window| {
                window.iter().try_fold(0u128, |value, &base| {
                    let digit = match base {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        b'T' => 3,
                        _ => return None,
                    };
                    Some((value << 2) | digit)
                })
            })
            .filter(|&value| value < VALUE_LIMIT)
            .collect();
        let mut actual = FxHashSet::default();
        add_sequence_kmers(&mixed, &mut actual);
        assert_eq!(actual, expected);
    }

    #[test]
    fn subtraction_keeps_only_output_kmers_absent_from_agc() {
        let remaining = RwLock::new(FxHashSet::from_iter([0, 1, 2]));
        let mut agc_batch = FxHashSet::from_iter([1, 2, 3]);
        remove_present_kmers(&mut agc_batch, &remaining);
        assert!(agc_batch.is_empty());
        assert_eq!(remaining.into_inner().unwrap(), FxHashSet::from_iter([0]));
    }
}
