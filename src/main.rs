use std::ops::Range;

use gxhash::gxhash128;
use simd_minimizers::packed_seq::AsciiSeq;

fn main() {
    let path = "/home/philae/git/eth/data/hprcv2.agc";
    use ragc_core::{Decompressor, DecompressorConfig};

    // Open an archive
    let config = DecompressorConfig::default();
    let mut decompressor = Decompressor::open(path, config).unwrap();

    // List available samples
    let samples = decompressor.list_samples();
    // println!("Found {} samples", samples.len());
    // eprintln!("samples: {samples:?}");

    let mut seen = std::collections::HashSet::new();
    let k = 40; // overlap
    let w = 100; // at most w apart
    eprintln!("k: {k}   w: {w}");
    let l = k + w - 1; // syncmer length
    let syncmers = simd_minimizers::closed_syncmers(k, w);

    let mut num_syncmers = 0;
    let mut taken = 0;
    let mut skipped = 0;
    let mut short = 0;
    let mut short_bp = 0;
    let mut prefix_bp = 0;
    let mut suffix_bp = 0;
    let mut output_bp = 0;
    let mut skipped_bp = 0;
    let mut input_bp = 0;
    let mut num_ranges = 0;

    let mut cnt = 0;
    let mut syncmer_poss = vec![];
    for (si, sample) in samples.iter().enumerate() {
        let mut new_ranges = 0;
        let mut new_bp = 0;
        let mut new_taken = 0;
        eprintln!("sample {si}: {}", sample);
        let contigs = decompressor.list_contigs(sample).unwrap();
        for contig in &contigs {
            cnt += 1;
            let seq = decompressor.get_contig(&sample, &contig).unwrap();
            input_bp += seq.len();

            // eprintln!("{cnt:>3} ({:>10}): contig {contig}", seq.len());

            // Minimizer-based PFP instead?
            syncmer_poss.clear();
            syncmers.run(AsciiSeq(&seq), &mut syncmer_poss);

            if syncmer_poss.is_empty() {
                short += 1;
                short_bp += seq.len();
                continue;
            }

            // sanity check
            for &[p, q] in syncmer_poss.array_windows() {
                assert!(p < q, "{p} < {q}");
                assert!(q <= p + w as u32, "{q} <= {p} + {w}");
                // assert!(p as usize + w <= seq.len(), "{p} + {w} <= {}", seq.len());
                assert!(p as usize + l <= seq.len(), "{p} + {l} <= {}", seq.len());
            }
            // assert!(*syncmer_poss.last().unwrap() as usize + l >= seq.len());

            // Emit the prefix.
            let mut active = 0..syncmer_poss[0] as usize + k - 1;
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
            let prefix = &seq[..syncmer_poss[0] as usize + k - 1];
            prefix_bp += prefix.len();
            // Emit the syncmers.
            for &p in syncmer_poss.iter() {
                let p = p as usize;
                num_syncmers += 1;
                let syncmer = &seq[p..p + l];
                // output_bp += syncmer.len();
                let hash = gxhash128(syncmer, 0);
                if seen.insert(hash) {
                    // eprintln!("take {p} hash={hash}");
                    taken += 1;
                    new_taken += 1;
                    push(p..p + l);
                } else {
                    // eprintln!("skip {p} hash={hash}");
                    skipped += 1;
                }
            }
            // Emit the suffix.
            let suffix_range = syncmer_poss[syncmer_poss.len() - 1] as usize + w..seq.len();
            let suffix = &seq[suffix_range.clone()];
            suffix_bp += suffix.len();
            push(suffix_range);
            num_ranges += 1;
            new_ranges += 1;
            output_bp += active.len();
            new_bp += active.len();
        }
        eprintln!(
            "  new bp:           {:>8.3} Mbp ({:3.1} bp/range)",
            new_bp as f32 / 1e6,
            new_bp as f32 / new_ranges as f32
        );
        eprintln!(
            "  new taken:        {:>8.3} M   ({:3.1} /range)",
            new_taken as f32 / 1e6,
            new_taken as f32 / new_ranges as f32
        );
        eprintln!(
            "  syncmers taken:   {:>8.3} M   ({:3.1}%)",
            taken as f32 / 1e6,
            100.0 * taken as f32 / num_syncmers as f32
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
}
