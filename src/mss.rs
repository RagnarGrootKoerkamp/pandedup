//! Compute a masked superstring using kmercamel's bigreedy algorithm.
//!
//! See: https://github.com/OndrejSladky/kmercamel, and the paper:
//!
//! Ondřej Sladký, Pavel Veselý, and Karel Břinda:
//! Masked superstrings as a unified framework for textual k-mer set representations.
//! bioRxiv 2023.02.01.526717, 2023. https://doi.org/10.1101/2023.02.01.526717

use seq_hash::packed_seq::{self, Seq};
use std::sync::{
    Mutex,
    atomic::{AtomicUsize, Ordering},
};
use tracing::info;
use voracious_radix_sort::{RadixSort, Radixable};

#[derive(Clone, Copy)]
struct SortEntry {
    key: u128,
    index: usize,
    reverse: usize,
}

impl PartialEq for SortEntry {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key
    }
}

impl PartialOrd for SortEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.key.cmp(&other.key))
    }
}

impl Radixable<u128> for SortEntry {
    type Key = u128;

    fn key(&self) -> Self::Key {
        self.key
    }
}

#[derive(Clone, Copy)]
struct MergePart {
    start: [usize; 4],
    end: [usize; 4],
    len: usize,
}

/// Merge 4 slices of `input` into `output`, using the given `mask` to compare keys.
fn merge_head_part(input: &[SortEntry], output: &mut [SortEntry], mask: u128, part: MergePart) {
    let mut cursors = part.start;
    for slot in output {
        let mut best_base = 0;
        let mut best_key = u128::MAX;
        for base in 0..4 {
            if cursors[base] < part.end[base] {
                let key = input[cursors[base]].key & mask;
                if key <= best_key {
                    best_base = base;
                    best_key = key;
                }
            }
        }
        let mut entry = input[cursors[best_base]];
        entry.key = best_key;
        *slot = entry;
        cursors[best_base] += 1;
    }
}

/// Drop the leading base from sorted head keys by merging the four sorted
/// ranges sharing that base. `overlap` is the length of the new keys.
fn resort_heads(
    heads: &mut Vec<SortEntry>,
    scratch: &mut Vec<SortEntry>,
    overlap: usize,
    threads: usize,
) {
    if heads.is_empty() {
        return;
    }
    let shift = 2 * overlap;
    let mask = (1u128 << shift) - 1;
    let mut bounds = [0; 5];
    let mut end = 0;
    for base in 0..4 {
        while end < heads.len() && (heads[end].key >> shift) == base as u128 {
            end += 1;
        }
        bounds[base + 1] = end;
    }
    debug_assert_eq!(end, heads.len());

    let chunks = threads.max(1).min(heads.len());
    let width = (1u128 << shift).div_ceil(chunks as u128);
    let mut cursors = [bounds[0], bounds[1], bounds[2], bounds[3]];
    let mut parts = Vec::with_capacity(chunks);
    for chunk in 1..=chunks {
        let upper = width * chunk as u128;
        let start = cursors;
        for base in 0..4 {
            cursors[base] += heads[cursors[base]..bounds[base + 1]]
                .partition_point(|entry| entry.key & mask < upper);
        }
        let len = (0..4).map(|base| cursors[base] - start[base]).sum();
        parts.push(MergePart {
            start,
            end: cursors,
            len,
        });
    }
    debug_assert_eq!(
        parts.iter().map(|part| part.len).sum::<usize>(),
        heads.len()
    );

    scratch.resize(
        heads.len(),
        SortEntry {
            key: 0,
            index: 0,
            reverse: 0,
        },
    );
    let input = heads.as_slice();
    std::thread::scope(|scope| {
        let mut remaining = scratch.as_mut_slice();
        for part in parts {
            let (output, rest) = remaining.split_at_mut(part.len);
            remaining = rest;
            if !output.is_empty() {
                scope.spawn(move || merge_head_part(input, output, mask, part));
            }
        }
    });
    std::mem::swap(heads, scratch);
}

struct Contig {
    /// The head of the contig, in fwd and rc direction
    kmer_in: [u128; 2],
    /// The tail of the contig, in fwd and rc direction
    kmer_out: [u128; 2],
    /// Neighbors at the forward head (0) and forward tail (1).
    nbs: [Option<Link>; 2],
}

struct Link {
    index: usize,
    reverse: bool,
    overlap: u8,
}

/// Input: unitigs, simplitigs=pathtigs, eulertigs (or (greedy) matchtigs).
///
/// Output: a superstring that contains all the input contigs, with newly
/// introduced "stitching" kmers masked out by lowercasing them.
///
/// If the input contains duplicate kmers, in case of (greedy) matchtigs, those
/// will be preserved in the output.
pub fn masked_superstring(k: usize, contigs: &Vec<Vec<u8>>) -> Vec<u8> {
    info!("masked_superstring: k={}, contigs={}", k, contigs.len());
    assert!((1..=64).contains(&k));

    // Each contig is stored as a fwd and rc kmer for the head, and one for the tail.
    // [[fwd in, fwd out], [rc in, rc out]]
    let mut ends: Vec<Contig> = vec![];
    for c in contigs {
        // TODO handle short contigs
        if c.len() < k {
            continue;
        }
        let seq = packed_seq::AsciiSeq(c);
        let head_fw = seq.read_kmer_u128(k, 0);
        let head_rc = seq.read_revcomp_kmer_u128(k, 0);
        let tail_fw = seq.read_kmer_u128(k, c.len() - k);
        let tail_rc = seq.read_revcomp_kmer_u128(k, c.len() - k);
        ends.push(Contig {
            kmer_in: [head_fw, tail_rc],
            kmer_out: [tail_fw, head_rc],
            nbs: [None, None],
        });
    }

    let mut total_merged = 0;
    let mut total_len = contigs.iter().map(|c| c.len()).sum::<usize>();
    let threads = std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1);
    // A physical end has an entry in each sorted vector. Keep its used bit
    // separately so both entries can be removed without reading `ends` again.
    let word_bits = usize::BITS as usize;
    let used: Vec<_> = (0..(2 * ends.len()).div_ceil(word_bits))
        .map(|_| AtomicUsize::new(0))
        .collect();

    // Sorting full tail k-mers also sorts every shorter prefix used below.
    let mut tails = Vec::with_capacity(2 * ends.len());
    for (index, contig) in ends.iter().enumerate() {
        for reverse in 0..2 {
            tails.push(SortEntry {
                key: contig.kmer_out[reverse],
                index,
                reverse,
            });
        }
    }
    info!("sorting {} tails..", tails.len());
    tails.voracious_mt_sort(threads);

    let first_overlap = k - 1;
    let first_head_mask = (1u128 << (2 * first_overlap)) - 1;
    let mut heads = Vec::with_capacity(2 * ends.len());
    for (index, contig) in ends.iter().enumerate() {
        for reverse in 0..2 {
            heads.push(SortEntry {
                key: contig.kmer_in[reverse] & first_head_mask,
                index,
                reverse,
            });
        }
    }
    info!("sorting {} heads..", heads.len());
    heads.voracious_mt_sort(threads);
    let mut head_scratch = Vec::with_capacity(heads.len());

    for overlap in (0..=first_overlap).rev() {
        info!("overlap: {}", overlap);
        let tail_shift = 2 * (k - overlap);
        if overlap < first_overlap {
            resort_heads(&mut heads, &mut head_scratch, overlap, threads);
        }

        info!("merging..");
        let chunks = threads.max(1);
        let width = (1u128 << (2 * overlap)).div_ceil(chunks as u128);
        let mut head_start = 0;
        let mut tail_start = 0;
        let mut ranges = Vec::with_capacity(chunks);
        for chunk in 1..=chunks {
            let upper = width * chunk as u128;
            let head_end =
                head_start + heads[head_start..].partition_point(|head| head.key < upper);
            let tail_end = tail_start
                + tails[tail_start..].partition_point(|tail| {
                    tail.key.checked_shr(tail_shift as u32).unwrap_or(0) < upper
                });
            ranges.push((head_start..head_end, tail_start..tail_end));
            head_start = head_end;
            tail_start = tail_end;
        }
        debug_assert_eq!(head_start, heads.len());
        debug_assert_eq!(tail_start, tails.len());

        // Each mutex owns a disjoint mutable slice. Lock both slices in index
        // order when a link crosses shards, so workers cannot deadlock.
        let shard_len = ends.len().div_ceil(chunks.saturating_mul(16)).max(1);
        let shards: Vec<_> = ends.chunks_mut(shard_len).map(Mutex::new).collect();
        let merged = std::thread::scope(|scope| {
            let mut workers = Vec::new();
            let mut remaining_heads = heads.as_mut_slice();
            let mut remaining_tails = tails.as_mut_slice();
            for (head_range, tail_range) in ranges {
                let (head_slice, rest) = remaining_heads.split_at_mut(head_range.len());
                remaining_heads = rest;
                let (tail_slice, rest) = remaining_tails.split_at_mut(tail_range.len());
                remaining_tails = rest;
                if head_range.is_empty() || tail_range.is_empty() {
                    continue;
                }
                let shards = &shards;
                let used = &used;
                workers.push(scope.spawn(move || {
                    let mut merged = 0;
                    let mut i = 0;
                    let mut j = 0;
                    while i < head_slice.len() && j < tail_slice.len() {
                        let head = head_slice[i];
                        let tail = tail_slice[j];
                        let tail_key = tail.key.checked_shr(tail_shift as u32).unwrap_or(0);
                        if head.key < tail_key {
                            i += 1;
                            continue;
                        }
                        if head.key > tail_key {
                            j += 1;
                            continue;
                        }
                        let tail_slot = 1 - tail.reverse;
                        if head.index == tail.index && head.reverse == tail_slot {
                            j += 1;
                            continue;
                        }

                        let head_shard = head.index / shard_len;
                        let tail_shard = tail.index / shard_len;
                        let (head_used, tail_used) = if head_shard == tail_shard {
                            let mut shard = shards[head_shard].lock().unwrap();
                            let head_local = head.index % shard_len;
                            let tail_local = tail.index % shard_len;
                            let head_used = shard[head_local].nbs[head.reverse].is_some();
                            let tail_used = shard[tail_local].nbs[tail_slot].is_some();
                            if !head_used && !tail_used {
                                shard[head_local].nbs[head.reverse] = Some(Link {
                                    index: tail.index,
                                    reverse: tail.reverse != 0,
                                    overlap: overlap as u8,
                                });
                                shard[tail_local].nbs[tail_slot] = Some(Link {
                                    index: head.index,
                                    reverse: head.reverse != 0,
                                    overlap: overlap as u8,
                                });
                            }
                            (head_used, tail_used)
                        } else {
                            let (lower, upper) = shards.split_at(head_shard.max(tail_shard));
                            let mut lower = lower[head_shard.min(tail_shard)].lock().unwrap();
                            let mut upper = upper[0].lock().unwrap();
                            let (head_end, tail_end) = if head_shard < tail_shard {
                                (
                                    &mut lower[head.index % shard_len],
                                    &mut upper[tail.index % shard_len],
                                )
                            } else {
                                (
                                    &mut upper[head.index % shard_len],
                                    &mut lower[tail.index % shard_len],
                                )
                            };
                            let head_used = head_end.nbs[head.reverse].is_some();
                            let tail_used = tail_end.nbs[tail_slot].is_some();
                            if !head_used && !tail_used {
                                head_end.nbs[head.reverse] = Some(Link {
                                    index: tail.index,
                                    reverse: tail.reverse != 0,
                                    overlap: overlap as u8,
                                });
                                tail_end.nbs[tail_slot] = Some(Link {
                                    index: head.index,
                                    reverse: head.reverse != 0,
                                    overlap: overlap as u8,
                                });
                            }
                            (head_used, tail_used)
                        };
                        if !head_used && !tail_used {
                            let head_port = 2 * head.index + head.reverse;
                            let tail_port = 2 * tail.index + tail_slot;
                            used[head_port / word_bits]
                                .fetch_or(1usize << (head_port % word_bits), Ordering::Relaxed);
                            used[tail_port / word_bits]
                                .fetch_or(1usize << (tail_port % word_bits), Ordering::Relaxed);
                            head_slice[i].index = usize::MAX;
                            tail_slice[j].index = usize::MAX;
                            merged += 1;
                            i += 1;
                            j += 1;
                        } else {
                            if head_used {
                                head_slice[i].index = usize::MAX;
                                i += 1;
                            }
                            if tail_used {
                                tail_slice[j].index = usize::MAX;
                                j += 1;
                            }
                        }
                    }
                    merged
                }));
            }
            workers
                .into_iter()
                .map(|worker| worker.join().unwrap())
                .sum::<usize>()
        });
        drop(shards);
        if merged != 0 {
            eprintln!("Retain..");
            heads.retain(|head| head.index != usize::MAX);
            tails.retain(|tail| tail.index != usize::MAX);
        }
        total_merged += merged;
        total_len -= merged * overlap;
        info!(
            "overlap {overlap} merged {:>9} total merged {total_merged:>9} remaining {:>9} total len {total_len:>11}",
            merged,
            ends.len() - total_merged
        );
    }
    assert_eq!(total_merged, ends.len());
    // TODO break cycles

    unimplemented!("break cycles");

    let mut num_cycles = 0;
    let mut done = vec![false; ends.len()];
    for i in 0..ends.len() {
        if done[i] {
            continue;
        }
        let mut j = i;
        let mut reverse = false;
        loop {
            done[j] = true;
            let link = ends[j].nbs[!reverse as usize].as_ref().unwrap();
            j = link.index;
            reverse = link.reverse;
            if j == i {
                break;
            }
        }
        num_cycles += 1;
    }
    eprintln!("Number of cycles: {num_cycles}");

    todo!()
}
