use clap::Parser;
use either::Either;
use ragc_core::{Decompressor, DecompressorConfig};
use std::{
    io::{BufWriter, Write},
    ops::Range,
    path::PathBuf,
    sync::{atomic::AtomicUsize, Mutex, RwLock},
};

use gxhash::gxhash128;
use seq_hash::AntiLexHasher;
use simd_minimizers::packed_seq::AsciiSeq;

/// Build a non-minimal SPSS (spectrum-preserving string set, or k-mer spectrum) from an .agc file.
#[derive(clap::Parser)]
struct Args {
    /// Input .agc file.
    input: PathBuf,
    /// Output path. Defaults to `input.dedup.fa.zst`.
    #[clap(short)]
    output: Option<PathBuf>,
    /// Build a k-mer spectrum.
    #[clap(short, default_value = "64")]
    k: usize,
    /// Window-size for minimizer-phrases. Small w shrink the output but need more memory.
    #[clap(short, default_value = "100")]
    w: usize,
    /// Number of threads to use. Defaults to the number of logical cores.
    #[clap(short = 'j', long)]
    threads: Option<usize>,

    /// Dedup across reverse-complements.
    #[clap(long)]
    canonical: bool,

    /// Minimizer length for phrases.
    #[clap(long, default_value = "8")]
    mini_k: usize,
}

#[derive(Default, Clone, Copy, derive_more::AddAssign)]
struct Stats {
    input_bp: usize,
    total_phrases: usize,
    output_bp: usize,
    unique_phrases: usize,
    output_contigs: usize,
}

fn main() {
    let Args {
        input,
        output,
        k,
        w,
        canonical,
        mini_k,
        threads,
    } = Args::parse();

    eprintln!("k: {k}   w: {w}");
    // Open an archive
    let config = DecompressorConfig { verbosity: 0 };
    let mut decompressor = Decompressor::open(&input.to_string_lossy(), config).unwrap();

    // List available samples
    let samples = decompressor.list_samples();
    let samples: Vec<_> = samples
        .into_iter()
        .map(|sample| {
            let contigs = decompressor.list_contigs(&sample).unwrap();
            (sample, contigs)
        })
        .collect();

    println!("Found {} samples", samples.len());

    // TODO: zstd output
    if let Some(output) = &output {
        assert!(
            output.extension().unwrap() == "zst",
            "Output file must have .zst extension"
        );
    }
    let stats_and_writer = {
        let output_path = output.unwrap_or_else(|| input.with_extension("dedup.fa"));
        let buf_writer =
            BufWriter::with_capacity(1 << 20, std::fs::File::create(output_path).unwrap());
        let writer = zstd::Encoder::new(buf_writer, 0).unwrap().auto_finish();
        Mutex::new((Stats::default(), writer))
    };

    let hasher = AntiLexHasher::<false>::new(mini_k);
    let scheme = if canonical {
        Either::Left(simd_minimizers::canonical_minimizers(mini_k, w))
    } else {
        // let scheme = simd_minimizers::closed_syncmers(k, w);
        // let scheme = simd_minimizers::minimizers(k, w);
        // Either::Right(simd_minimizers::minimizers(mini_k, w))
        Either::Right(simd_minimizers::minimizers(mini_k, w).hasher(&hasher))
    };

    let seen: [_; 256] = std::array::from_fn(|_i| RwLock::new(std::collections::HashSet::new()));

    let next = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        let threads = threads.unwrap_or_else(|| num_cpus::get_physical());
        for _t in 0..threads {
            let decompressor = &decompressor;
            let samples = &samples;
            let next = &next;
            let seen = &seen;
            let scheme = &scheme;
            let stats_and_writer = &stats_and_writer;
            scope.spawn(|| loop {
                let idx = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if idx >= samples.len() {
                    break;
                }
                let mut local_stats = Stats::default();

                let (sample, contigs) = &samples[idx];

                let start = std::time::Instant::now();
                let seqs: Vec<_> = contigs
                    .iter()
                    .map(|contig| {
                        let mut contig = decompressor.get_contig(&sample, &contig).unwrap();
                        contig
                            .iter_mut()
                            .for_each(|b| *b = b"ACGT"[(*b as usize) % 4]);
                        contig
                    })
                    .collect();
                local_stats.input_bp = seqs.iter().map(|seq| seq.len()).sum();
                let mid = std::time::Instant::now();
                let (seqs, poss): (Vec<_>, Vec<_>) = seqs
                    .into_iter()
                    .map(|seq| {
                        let mut positions = vec![];
                        match scheme {
                            Either::Left(scheme) => {
                                drop(scheme.run(AsciiSeq(&seq), &mut positions))
                            }
                            Either::Right(scheme) => {
                                drop(scheme.run(AsciiSeq(&seq), &mut positions))
                            }
                        }

                        (seq, positions)
                    })
                    .unzip();
                let end = std::time::Instant::now();

                let rc_seqs: Vec<_> = if canonical {
                    seqs.iter()
                        .map(|seq| {
                            let mut rc_seq = vec![];
                            rc_seq.extend(seq.iter().rev().map(|bp| 3 - bp));
                            rc_seq
                        })
                        .collect()
                } else {
                    vec![]
                };

                let mut phrases = vec![];
                for (i, positions) in poss.into_iter().enumerate() {
                    if positions.is_empty() {
                        continue;
                    }

                    let seq = &seqs[i];
                    let empty = vec![];
                    let rc_seq = rc_seqs.get(i).unwrap_or(&empty);

                    let with_hash = |p, q| {
                        let phrase = &seq[p..q];
                        let hash = gxhash128(phrase, 0);
                        if canonical {
                            let rc_phrase = &rc_seq[seq.len() - q..seq.len() - p];
                            let rc_hash = gxhash128(rc_phrase, 0);
                            (i, p, q, hash + rc_hash)
                        } else {
                            (i, p, q, hash)
                        }
                    };

                    // Prefix phrase.
                    if positions[0] > 0 {
                        phrases.push(with_hash(0, (positions[0] as usize + k).min(seq.len())));
                    }
                    for &[p, q] in positions.array_windows::<2>() {
                        // Skip non-forward minimizer pairs.
                        if q < p {
                            continue;
                        }
                        let p = p as usize;
                        let q = (q as usize + k).min(seq.len());
                        let (i, p, q, hash) = with_hash(p, q);
                        phrases.push((i, p, q, hash));
                    }
                    // Suffix phrase.
                    if (*positions.last().unwrap() as usize) + k < seq.len() {
                        phrases.push(with_hash(*positions.last().unwrap() as usize, seq.len()));
                    }
                }
                local_stats.total_phrases = phrases.len();

                let end2 = std::time::Instant::now();

                let mut order = (0..256).map(|_| vec![]).collect::<Vec<_>>();
                for (i, p) in phrases.iter().enumerate() {
                    let part = p.3 as u8;
                    order[part as usize].push(i as u32);
                }
                let end3 = std::time::Instant::now();
                for part in 0..=255 {
                    // Read-only filter phase
                    {
                        let seen = seen[part as usize].read().unwrap();
                        for &idx in &order[part as usize] {
                            let (i, _p, _q, hash) = &mut phrases[idx as usize];
                            assert!(*hash as u8 == part);
                            if seen.contains(hash) {
                                *i = usize::MAX;
                            }
                        }
                    }

                    // Try to write missing
                    {
                        let mut seen = seen[part as usize].write().unwrap();
                        for &idx in &order[part as usize] {
                            let (i, _p, _q, hash) = &mut phrases[idx as usize];
                            if *i != usize::MAX {
                                if !seen.insert(*hash) {
                                    *i = usize::MAX;
                                }
                            }
                        }
                    }
                }
                let end4 = std::time::Instant::now();
                phrases.retain(|(i, _p, _q, _hash)| *i != usize::MAX);
                local_stats.unique_phrases = phrases.len();
                let end5 = std::time::Instant::now();

                eprintln!(
                    "push sample {idx:>3} ({:3.1} Gbp): \
                     read: {:5.2?}s minis: {:5.2?}s phrases: \
                     {:5.2?}s sort: {:5.2?}s lookups: {:5.2?}s sort: {:5.2?}s",
                    (local_stats.input_bp as f32) / 1e9,
                    (mid - start).as_secs_f32(),
                    (end - mid).as_secs_f32(),
                    (end2 - end).as_secs_f32(),
                    (end3 - end2).as_secs_f32(),
                    (end4 - end3).as_secs_f32(),
                    (end5 - end4).as_secs_f32(),
                );

                let mut stats_and_writer = stats_and_writer.lock().unwrap();
                let (global_stats, writer) = &mut *stats_and_writer;

                // Process
                let start = std::time::Instant::now();
                for phrases in phrases.chunk_by(|l, r| l.0 == r.0) {
                    let seq = &seqs[phrases[0].0];

                    let mut active = 0..0;

                    let mut push = |range: Range<usize>| {
                        if range.start <= active.end {
                            // End can decrease for non-forward canonical minimizers.
                            active.end = active.end.max(range.end);
                        } else {
                            writer.write_all(b">\n").unwrap();
                            writer.write_all(&seq[active.clone()]).unwrap();
                            writer.write_all(b"\n").unwrap();

                            local_stats.output_contigs += 1;
                            local_stats.output_bp += active.len();

                            active = range;
                        }
                    };

                    // Emit the syncmers.
                    for &(_i, p, q, _hash) in phrases {
                        push(p..q);
                    }

                    push(usize::MAX..usize::MAX);
                    assert!(active.len() == 0);
                }

                *global_stats += local_stats;

                eprintln!("process sample {idx}: {:?}", start.elapsed());
                eprintln!(
                    "  new bp:            {:>8.3} Mbp ({:3.1} bp/contig)",
                    local_stats.output_bp as f32 / 1e6,
                    local_stats.output_bp as f32 / local_stats.output_contigs as f32
                );
                eprintln!(
                    "  new phrases:       {:>8.3} M   ({:3.1} /contig; {:4.1}%)",
                    local_stats.unique_phrases as f32 / 1e6,
                    local_stats.unique_phrases as f32 / local_stats.output_contigs as f32,
                    local_stats.unique_phrases as f32 / local_stats.total_phrases as f32 * 100.0
                );
                eprintln!(
                    "  unique phrases:    {:>8.3} M   ({:3.1}%)",
                    global_stats.unique_phrases as f32 / 1e6,
                    100.0 * global_stats.unique_phrases as f32 / global_stats.total_phrases as f32
                );

                eprintln!(
                    "  num_contigs:       {:>8.3} M",
                    global_stats.output_contigs as f32 / 1e6
                );
                eprintln!(
                    "  output_bp:         {:>8.3} Gbp ({:3.1}%)",
                    global_stats.output_bp as f32 / 1e9,
                    100.0 * global_stats.output_bp as f32 / global_stats.input_bp as f32
                );
                eprintln!();
            });
        }
    });
}
