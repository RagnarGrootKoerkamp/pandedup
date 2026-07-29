#![allow(unused)]
use std::{
    io::{BufWriter, Write},
    ops::Range,
    sync::{atomic::AtomicUsize, Mutex, RwLock},
};

use gxhash::gxhash128;
use seq_hash::{packed_seq::u32x8, KmerHasher};
use simd_minimizers::packed_seq::AsciiSeq;

const THREADS: usize = 6;

fn main() {
    let path = "/home/philae/git/eth/data/hprcv2.agc";
    use ragc_core::{Decompressor, DecompressorConfig};

    // Open an archive
    let config = DecompressorConfig::default();
    let mut decompressor = Decompressor::open(path, config).unwrap();

    // List available samples
    let samples = decompressor.list_samples();
    let samples: Vec<_> = samples
        .into_iter()
        .map(|sample| {
            let contigs = decompressor.list_contigs(&sample).unwrap();
            (sample, contigs)
        })
        .collect();
    // println!("Found {} samples", samples.len());
    // eprintln!("samples: {samples:?}");

    let mut seen = RwLock::new(std::collections::HashSet::new());
    let k = 40usize; // overlap
    let mini_k = 8; // k for minimizers. smaller to make it more stable
    let w = 100; // at most w apart
    eprintln!("k: {k}   w: {w}");
    let l = k + w - 1; // syncmer length

    let mut total_filtered_phrases = 0usize;
    let mut taken_phrases = 0usize;
    let mut skipped = 0usize;
    let mut short = 0usize;
    let mut short_bp = 0;
    let mut prefix_bp = 0;
    let mut suffix_bp = 0;
    let mut output_bp = 0;
    let mut skipped_bp = 0;
    let mut input_bp = 0;
    let mut num_ranges = 0usize;

    // // output text
    // let mut output = vec![];
    // // end position of each sequence in output.
    // let mut ends = vec![0];
    let mut output = BufWriter::with_capacity(
        1 << 20,
        std::fs::File::create(format!("deduped_k{k}_w{w}_mk{mini_k}.fa")).unwrap(),
    );

    let next = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        let (write, read) = std::sync::mpsc::sync_channel(THREADS);
        for _t in 0..THREADS {
            let write = write.clone();
            let decompressor = &decompressor;
            let samples = &samples;
            let next = &next;
            let seen = &seen;
            // let scheme = simd_minimizers::closed_syncmers(k, w);
            // let scheme = simd_minimizers::minimizers(k, w);
            let scheme = simd_minimizers::minimizers(mini_k, w);
            // let scheme = seq_hash::NtHasher::<false>::new(mini_k);
            scope.spawn(move || loop {
                let idx = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if idx >= samples.len() {
                    break;
                }

                let (sample, contigs) = &samples[idx];

                let start = std::time::Instant::now();
                let seqs: Vec<_> = contigs
                    .iter()
                    .map(|contig| {
                        let mut contig = decompressor
                            .get_contig(&sample, &contig)
                            .unwrap();
                        contig
                            .iter_mut()
                            .for_each(|b| *b = b"ACGT"[(*b as usize) % 4]);
                        contig
                    })
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

                // Pre-filter already seen phrases.
                let seen = seen.read().unwrap();
                let locking = std::time::Instant::now();
                let seqs_and_phrases = seqs_and_poss
                    .into_iter()
                    .map(|(seq, positions)| {
                        let mut phrases = vec![];
                        if positions.is_empty() {
                            return (seq, phrases);
                        }
                        // Prefix phrase.
                        if positions[0] > 0 {
                            phrases.push((0, (positions[0] as usize + k).min(seq.len())));
                        }
                        for &[p, q] in positions.array_windows::<2>() {
                            let p = p as usize;
                            let q = (q as usize + k).min(seq.len());
                            let phrase = &seq[p..q];
                            let hash = gxhash128(phrase, 0);
                            if !seen.contains(&hash) {
                                phrases.push((p, q));
                            }
                        }
                        // Suffix phrase.
                        if *positions.last().unwrap() as usize + k < seq.len() {
                            phrases.push((*positions.last().unwrap() as usize, seq.len()));
                        }
                        (seq, phrases)
                    })
                    .collect::<Vec<_>>();
                drop(seen);

                let end2 = std::time::Instant::now();

                eprintln!(
                    "push sample {idx:>3} ({:3.1} Gbp): read: {:?} minis: {:?} read-lock: {:?} filter: {:?}",
                    len as f32 / 1e9,
                    mid - start,
                    end - mid,
                    locking-end,
                    end2 - locking
                );
                write.send(seqs_and_phrases).unwrap();
            });
        }
        drop(write);

        // Read from the queue in the current thread.

        for (si, sample) in read.into_iter().enumerate() {
            let mut new_ranges = 0;
            let mut new_bp = 0;
            let mut new_phrases = 0;
            let mut seen = seen.write().unwrap();
            let start = std::time::Instant::now();
            for (seq, phrases) in sample {
                input_bp += seq.len();
                if phrases.is_empty() {
                    short += 1;
                    short_bp += seq.len();
                    continue;
                }

                let mut active = 0..0;

                let mut push = |range: Range<usize>| {
                    assert!(range.end >= active.end);
                    if range.start <= active.end {
                        active.end = range.end;
                    } else {
                        output.write_all(b">\n");
                        output.write_all(&seq[active.clone()]).unwrap();
                        output.write_all(b"\n");

                        new_ranges += 1;
                        num_ranges += 1;
                        new_bp += active.len();
                        output_bp += active.len();
                        skipped_bp += range.start - active.end;
                        active = range;
                    }
                };

                // Emit the syncmers.
                for (p, q) in phrases {
                    total_filtered_phrases += 1;
                    let phrase = &seq[p..q];
                    let hash = gxhash128(phrase, 0);
                    if seen.insert(hash) {
                        taken_phrases += 1;
                        new_phrases += 1;
                        push(p..q);
                    } else {
                        skipped += 1;
                    }
                }

                output.write_all(b">\n");
                output.write_all(&seq[active.clone()]).unwrap();
                output.write_all(b"\n");

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
                taken_phrases as f32 / 1e6,
                100.0 * taken_phrases as f32 / total_filtered_phrases as f32
            );

            // eprintln!("syncmers skipped: {skipped:>9}");
            // eprintln!("short:   {short}");
            // eprintln!("num_syncmers: {num_syncmers}");
            // eprintln!("prefix_bp:  {prefix_bp}");
            // eprintln!("suffix_bp:  {suffix_bp}");
            eprintln!("  num_ranges:        {:>8.3} M", num_ranges as f32 / 1e6);
            // eprintln!("input_bp:   {:>8.3} Gbp", input_bp as f32 / 1e9);
            eprintln!(
                "  output_bp:         {:>8.3} Gbp ({:3.1}%)",
                output_bp as f32 / 1e9,
                100.0 * output_bp as f32 / input_bp as f32
            );
            // eprintln!("skipped_bp: {:>8.3} Gbp", skipped_bp as f32 / 1e9);
            // eprintln!("short_bp:   {short_bp}");
            eprintln!();
        }
    });
}
