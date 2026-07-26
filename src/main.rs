#![allow(unused)]
use std::{ops::Range, sync::atomic::AtomicUsize};

use gxhash::gxhash128;
use seq_hash::{packed_seq::u32x8, KmerHasher};
use simd_minimizers::packed_seq::AsciiSeq;

fn main() {
    let path = "/home/philae/git/eth/data/hprcv2.agc";
    use ragc_core::{Decompressor, DecompressorConfig};

    // Open an archive
    let config = DecompressorConfig::default();
    let decompressor = Decompressor::open(path, config).unwrap();

    // List available samples
    let samples = decompressor.list_samples();
    drop(decompressor);
    // println!("Found {} samples", samples.len());
    // eprintln!("samples: {samples:?}");

    let mut seen = std::collections::HashSet::new();
    let k = 40usize; // overlap
    let mini_k = 8; // k for minimizers. smaller to make it more stable
    let w = 100; // at most w apart
    eprintln!("k: {k}   w: {w}");
    let l = k + w - 1; // syncmer length

    let mut num_phrases = 0usize;
    let mut phrases = 0usize;
    let mut skipped = 0usize;
    let mut short = 0usize;
    let mut short_bp = 0;
    let mut prefix_bp = 0;
    let mut suffix_bp = 0;
    let mut output_bp = 0;
    let mut skipped_bp = 0;
    let mut input_bp = 0;
    let mut num_ranges = 0usize;

    let next = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        let (write, read) = std::sync::mpsc::sync_channel(4);
        for _t in 0..4 {
            let write = write.clone();
            let config = DecompressorConfig::default();
            let mut decompressor = Decompressor::open(path, config).unwrap();
            let samples = &samples;
            let next = &next;
            // let scheme = simd_minimizers::closed_syncmers(k, w);
            // let scheme = simd_minimizers::minimizers(k, w);
            let scheme = simd_minimizers::minimizers(mini_k, w);
            // let scheme = seq_hash::NtHasher::<false>::new(mini_k);
            scope.spawn(move || loop {
                let idx = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if idx >= samples.len() {
                    break;
                }
                let sample = &samples[idx];
                let contigs = decompressor.list_contigs(sample).unwrap();

                let start = std::time::Instant::now();
                let seqs: Vec<_> = contigs
                    .iter()
                    .map(|contig| decompressor.get_contig(&sample, &contig).unwrap())
                    .collect();
                let len: usize = seqs.iter().map(|seq| seq.len()).sum();
                let mid = std::time::Instant::now();
                let seqs_and_poss: Vec<_> = seqs
                    .into_iter()
                    .map(|seq| {
                        let mut positions = vec![];
                        scheme.run(AsciiSeq(&seq), &mut positions);

                        // For PFP
                        // {
                        //     let threshold = (u32::MAX as u64 * 2 / (w as u64 + 1)) as u32;
                        //     let simd_threshold = u32x8::splat(threshold);
                        //     let hashes = scheme
                        //         .hash_kmers_simd(AsciiSeq(&seq), 1)
                        //         .map(|h| h.simd_lt(simd_threshold))
                        //         .collect();
                        //     for i in 0..hashes.len() {
                        //         if hashes[i] > 0 {
                        //             positions.push(i as u32);
                        //         }
                        //     }
                        // }

                        // sanity check
                        // for &[p, q] in positions.array_windows() {
                        //     assert!(p < q, "{p} < {q}");
                        //     assert!(q <= p + w as u32, "{q} <= {p} + {w}");
                        //     assert!(
                        //         p as usize + mini_k <= seq.len(),
                        //         "{p} + {w} <= {}",
                        //         seq.len()
                        //     );
                        //     // assert!(p as usize + w <= seq.len(), "{p} + {w} <= {}", seq.len());
                        //     // assert!(p as usize + l <= seq.len(), "{p} + {l} <= {}", seq.len());
                        // }
                        // assert!(*positions.last().unwrap() as usize + l >= seq.len());

                        (seq, positions)
                    })
                    .collect();
                let end = std::time::Instant::now();

                eprintln!(
                    "push sample {idx:>3} ({:3.1} Gbp): read: {:?} minis: {:?}",
                    len as f32 / 1e9,
                    mid - start,
                    end - mid
                );
                write.send(seqs_and_poss).unwrap();
            });
        }
        drop(write);

        // Read from the queue in the current thread.

        for (si, sample) in read.into_iter().enumerate() {
            let mut new_ranges = 0;
            let mut new_bp = 0;
            let mut new_phrases = 0;
            let start = std::time::Instant::now();
            for (seq, positions) in sample {
                input_bp += seq.len();
                if positions.is_empty() {
                    short += 1;
                    short_bp += seq.len();
                    continue;
                }

                // Emit the prefix.
                let mut active = 0..positions[0] as usize + k - 1;
                let mut push = |range: Range<usize>| {
                    assert!(range.end >= active.end);
                    if range.start <= active.end {
                        active.end = range.end;
                    } else {
                        new_ranges += 1;
                        num_ranges += 1;
                        new_bp += active.len();
                        output_bp += active.len();
                        skipped_bp += range.start - active.end;
                        active = range;
                    }
                };
                let prefix = &seq[..positions[0] as usize + k - 1];
                prefix_bp += prefix.len();
                // Emit the syncmers.
                for &[p, q] in positions.array_windows::<2>() {
                    let p = p as usize;
                    let q = q as usize;
                    num_phrases += 1;
                    // minimizer parse
                    let phrase = &seq[p..(q + k).min(seq.len())];
                    let hash = gxhash128(phrase, 0);
                    if seen.insert(hash) {
                        phrases += 1;
                        new_phrases += 1;
                        push(p..(q + k).min(seq.len()));
                    } else {
                        skipped += 1;
                    }
                }
                // Emit the suffix.
                let suffix_range = positions[positions.len() - 1] as usize..seq.len();
                let suffix = &seq[suffix_range.clone()];
                suffix_bp += suffix.len();
                push(suffix_range);
                num_ranges += 1;
                new_ranges += 1;
                output_bp += active.len();
                new_bp += active.len();
            }
            eprintln!("process sample {si}: {:?}", start.elapsed());
            eprintln!(
                "  new bp:           {:>8.3} Mbp ({:3.1} bp/range)",
                new_bp as f32 / 1e6,
                new_bp as f32 / new_ranges as f32
            );
            eprintln!(
                "  new phrases:        {:>8.3} M   ({:3.1} /range)",
                new_phrases as f32 / 1e6,
                new_phrases as f32 / new_ranges as f32
            );
            eprintln!(
                "  total phrases:   {:>8.3} M   ({:3.1}%)",
                phrases as f32 / 1e6,
                100.0 * phrases as f32 / num_phrases as f32
            );
            // eprintln!("syncmers skipped: {skipped:>9}");
            // eprintln!("short:   {short}");
            // eprintln!("num_syncmers: {num_syncmers}");
            // eprintln!("prefix_bp:  {prefix_bp}");
            // eprintln!("suffix_bp:  {suffix_bp}");
            eprintln!("  num_ranges:       {:>8.3} M", num_ranges as f32 / 1e6);
            // eprintln!("input_bp:   {:>8.3} Gbp", input_bp as f32 / 1e9);
            eprintln!(
                "  output_bp:        {:>8.3} Gbp ({:3.1}%)",
                output_bp as f32 / 1e9,
                100.0 * output_bp as f32 / input_bp as f32
            );
            // eprintln!("skipped_bp: {:>8.3} Gbp", skipped_bp as f32 / 1e9);
            // eprintln!("short_bp:   {short_bp}");
            eprintln!();
        }
    });
}
