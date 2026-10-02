//! Compute a masked superstring using kmercamel's bigreedy algorithm.
//!
//! See: https://github.com/OndrejSladky/kmercamel, and the paper:
//!
//! Ondřej Sladký, Pavel Veselý, and Karel Břinda:
//! Masked superstrings as a unified framework for textual k-mer set representations.
//! bioRxiv 2023.02.01.526717, 2023. https://doi.org/10.1101/2023.02.01.526717

use core::ops::{BitAnd, Mul, Range, Shl, Shr, Sub};
use packed_seq::{PackedSeqVec, SeqVec};
use rayon::prelude::*;
use std::sync::Mutex;
use tracing::{debug, info};
use voracious_radix_sort::{RadixKey, RadixSort, Radixable};

pub trait MssKey:
    RadixKey
    + Copy
    + Ord
    + Send
    + Sync
    + BitAnd<Output = Self>
    + Mul<Output = Self>
    + Shl<usize, Output = Self>
    + Shr<usize, Output = Self>
    + Sub<Output = Self>
{
    const BITS: usize;
    const ZERO: Self;
    const ONE: Self;
    const MAX: Self;
    fn from_usize(value: usize) -> Self;
    fn div_ceil(self, rhs: Self) -> Self;
    fn checked_shr(self, rhs: u32) -> Option<Self>;
    fn read_kmer(seq: &PackedSeqVec, k: usize, pos: usize) -> Self;
    fn read_revcomp_kmer(seq: &PackedSeqVec, k: usize, pos: usize) -> Self;
}

impl MssKey for u64 {
    const BITS: usize = 64;
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const MAX: Self = u64::MAX;
    fn from_usize(value: usize) -> Self {
        value as u64
    }
    fn div_ceil(self, rhs: Self) -> Self {
        u64::div_ceil(self, rhs)
    }
    fn checked_shr(self, rhs: u32) -> Option<Self> {
        u64::checked_shr(self, rhs)
    }
    fn read_kmer(seq: &PackedSeqVec, k: usize, pos: usize) -> Self {
        seq.read_kmer(k, pos)
    }
    fn read_revcomp_kmer(seq: &PackedSeqVec, k: usize, pos: usize) -> Self {
        seq.read_revcomp_kmer(k, pos)
    }
}

impl MssKey for u128 {
    const BITS: usize = 128;
    const ZERO: Self = 0;
    const ONE: Self = 1;
    const MAX: Self = u128::MAX;
    fn from_usize(value: usize) -> Self {
        value as u128
    }
    fn div_ceil(self, rhs: Self) -> Self {
        u128::div_ceil(self, rhs)
    }
    fn checked_shr(self, rhs: u32) -> Option<Self> {
        u128::checked_shr(self, rhs)
    }
    fn read_kmer(seq: &PackedSeqVec, k: usize, pos: usize) -> Self {
        seq.read_kmer_u128(k, pos)
    }
    fn read_revcomp_kmer(seq: &PackedSeqVec, k: usize, pos: usize) -> Self {
        seq.read_revcomp_kmer_u128(k, pos)
    }
}

#[derive(Clone, Copy)]
#[repr(C, packed(1))]
pub struct HeadOrTail<K: MssKey> {
    key: K,
    /// Twice the contig index, plus one for reverse complement; MAX means used.
    id: u32,
}

const _: () = {
    assert!(std::mem::align_of::<HeadOrTail<u64>>() == 1);
    assert!(std::mem::size_of::<HeadOrTail<u64>>() == 12);
    assert!(std::mem::align_of::<HeadOrTail<u128>>() == 1);
    assert!(std::mem::size_of::<HeadOrTail<u128>>() == 20);
};

impl<K: MssKey> PartialEq for HeadOrTail<K> {
    fn eq(&self, other: &Self) -> bool {
        let key = self.key;
        let other_key = other.key;
        key == other_key
    }
}

impl<K: MssKey> PartialOrd for HeadOrTail<K> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let key = self.key;
        let other_key = other.key;
        Some(key.cmp(&other_key))
    }
}

impl Radixable<u64> for HeadOrTail<u64> {
    type Key = u64;
    fn key(&self) -> Self::Key {
        self.key
    }
}

impl Radixable<u128> for HeadOrTail<u128> {
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
fn merge_head_part<K: MssKey>(
    input: &[HeadOrTail<K>],
    output: &mut [HeadOrTail<K>],
    mask: K,
    part: MergePart,
) {
    let mut cursors = part.start;
    for slot in output {
        let mut best_base = 0;
        let mut best_key = K::MAX;
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
fn resort_heads<K: MssKey>(
    heads: &mut Vec<HeadOrTail<K>>,
    scratch: &mut Vec<HeadOrTail<K>>,
    overlap: usize,
    threads: usize,
) {
    if heads.is_empty() {
        return;
    }
    let shift = 2 * overlap;
    let mask = (K::ONE << shift) - K::ONE;
    let mut bounds = [0; 5];
    for base in 1..4 {
        bounds[base] = heads.partition_point(|entry| entry.key >> shift < K::from_usize(base));
    }
    bounds[4] = heads.len();

    let chunks = threads.max(1).min(heads.len());
    let width = (K::ONE << shift).div_ceil(K::from_usize(chunks));
    let mut cursors = [bounds[0], bounds[1], bounds[2], bounds[3]];
    let mut parts = Vec::with_capacity(chunks);
    for chunk in 1..=chunks {
        let upper = width * K::from_usize(chunk);
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
        HeadOrTail {
            key: K::ZERO,
            id: 0,
        },
    );
    let input = heads.as_slice();
    rayon::scope(|scope| {
        let mut remaining = scratch.as_mut_slice();
        for part in parts {
            let (output, rest) = remaining.split_at_mut(part.len);
            remaining = rest;
            if !output.is_empty() {
                scope.spawn(move |_| merge_head_part(input, output, mask, part));
            }
        }
    });
    std::mem::swap(heads, scratch);
}

#[derive(Clone, Copy)]
struct Links {
    /// Neighbors at the forward head (0) and forward tail (1).
    nbs: [Option<Link>; 2],
}

#[derive(Clone, Copy)]
#[repr(C, packed(1))]
struct Link {
    /// Twice the contig index, plus one for reverse complement.
    id: u32,
    overlap: u8,
}

/// Input: unitigs, simplitigs=pathtigs, eulertigs (or (greedy) matchtigs).
///
/// Output: a superstring that contains all the input contigs, with newly
/// introduced "stitching" kmers masked out by lowercasing them.
///
/// If the input contains duplicate kmers, in case of (greedy) matchtigs, those
/// will be preserved in the output.
pub fn masked_superstring<K: MssKey>(
    k: usize,
    seq: PackedSeqVec,
    mut ranges: Vec<Range<usize>>,
) -> Vec<u8>
where
    HeadOrTail<K>: Radixable<K, Key = K>,
{
    info!("masked_superstring: k={}, contigs={}", k, ranges.len());
    assert!(k > 0 && k <= K::BITS / 2);

    info!("Collect tig ends");
    let mut total_merged = 0;
    let mut total_len = ranges.iter().map(|c| c.len()).sum::<usize>();
    let threads = rayon::current_num_threads();
    // Filter first so indices in both entry vectors match `ends` even when
    // the input contains contigs shorter than k.
    ranges.retain(|c| c.len() >= k);
    let num_contigs = ranges.len();
    // Leave room for the orientation bit and reserve u32::MAX as the sentinel.
    assert!(
        num_contigs <= (u32::MAX >> 1) as usize,
        "too many contigs for packed orientation IDs"
    );

    let empty_entry = HeadOrTail {
        key: K::ZERO,
        id: 0,
    };
    let mut heads = vec![empty_entry; 2 * num_contigs];
    let mut tails = vec![empty_entry; 2 * num_contigs];
    let first_overlap = k - 1;
    let first_head_mask = (K::ONE << (2 * first_overlap)) - K::ONE;
    ranges
        .into_par_iter()
        .zip(heads.par_chunks_mut(2))
        .zip(tails.par_chunks_mut(2))
        .enumerate()
        .for_each(|(index, ((r, head_entries), tail_entries))| {
            let head_fw = K::read_kmer(&seq, k, r.start);
            let head_rc = K::read_revcomp_kmer(&seq, k, r.start);
            let tail_fw = K::read_kmer(&seq, k, r.end - k);
            let tail_rc = K::read_revcomp_kmer(&seq, k, r.end - k);
            for reverse in [false, true] {
                let direction = usize::from(reverse);
                head_entries[direction] = HeadOrTail {
                    key: [head_fw, tail_rc][direction] & first_head_mask,
                    id: ((index as u32) << 1) | u32::from(reverse),
                };
                tail_entries[direction] = HeadOrTail {
                    key: [tail_fw, head_rc][direction],
                    id: ((index as u32) << 1) | u32::from(reverse),
                };
            }
        });
    drop(seq);

    let mut links = vec![Links { nbs: [None, None] }; num_contigs];

    // Sorting full tail k-mers also sorts every shorter prefix used below.
    info!("sorting {} heads..", heads.len());
    heads.voracious_mt_sort(threads);
    info!("sorting {} tails..", tails.len());
    tails.voracious_mt_sort(threads);
    info!("Reserve scratch");
    let mut scratch = vec![];

    let mut marked_heads = 0;
    let mut marked_tails = 0;

    for overlap in (0..=first_overlap).rev() {
        debug!("overlap: {}", overlap);
        let tail_shift = 2 * (k - overlap);
        if overlap < first_overlap {
            resort_heads(&mut heads, &mut scratch, overlap, threads);
        }

        debug!("merging..");
        let chunks = threads.max(1);
        let width = (K::ONE << (2 * overlap)).div_ceil(K::from_usize(chunks));
        let mut head_start = 0;
        let mut tail_start = 0;
        let mut ranges = Vec::with_capacity(chunks);
        for chunk in 1..=chunks {
            let upper = width * K::from_usize(chunk);
            let head_end = head_start
                + heads[head_start..].partition_point(|head| {
                    let key = head.key;
                    key < upper
                });
            let tail_end = tail_start
                + tails[tail_start..].partition_point(|tail| {
                    tail.key.checked_shr(tail_shift as u32).unwrap_or(K::ZERO) < upper
                });
            ranges.push((head_start..head_end, tail_start..tail_end));
            head_start = head_end;
            tail_start = tail_end;
        }
        debug_assert_eq!(head_start, heads.len());
        debug_assert_eq!(tail_start, tails.len());

        // Each mutex owns a disjoint mutable slice. Lock both slices in index
        // order when a link crosses shards, so workers cannot deadlock.
        let shard_len = links.len().div_ceil(chunks.saturating_mul(16)).max(1);
        let shards: Vec<_> = links.chunks_mut(shard_len).map(Mutex::new).collect();
        let (merged, newly_marked_heads, newly_marked_tails) = {
            let mut jobs = Vec::new();
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
                jobs.push((head_slice, tail_slice));
            }
            jobs.into_par_iter()
                .map(|(head_slice, tail_slice)| {
                    let mut merged = 0;
                    let mut marked_heads = 0;
                    let mut marked_tails = 0;
                    let mut i = 0;
                    let mut j = 0;
                    while i < head_slice.len() && j < tail_slice.len() {
                        let head = head_slice[i];
                        let tail = tail_slice[j];
                        if head.id == u32::MAX {
                            i += 1;
                            continue;
                        }
                        if tail.id == u32::MAX {
                            j += 1;
                            continue;
                        }
                        let tail_key = tail.key.checked_shr(tail_shift as u32).unwrap_or(K::ZERO);
                        let head_key = head.key;
                        if head_key < tail_key {
                            i += 1;
                            continue;
                        }
                        if head_key > tail_key {
                            j += 1;
                            continue;
                        }
                        if head.id == (tail.id ^ 1) {
                            j += 1;
                            continue;
                        }

                        let head_index = (head.id >> 1) as usize;
                        let tail_index = (tail.id >> 1) as usize;
                        let head_slot = (head.id & 1) as usize;
                        let tail_slot = ((tail.id ^ 1) & 1) as usize;
                        let head_shard = head_index / shard_len;
                        let tail_shard = tail_index / shard_len;
                        let (head_used, tail_used) = if head_shard == tail_shard {
                            let mut shard = shards[head_shard].lock().unwrap();
                            let head_local = head_index % shard_len;
                            let tail_local = tail_index % shard_len;
                            let head_used = shard[head_local].nbs[head_slot].is_some();
                            let tail_used = shard[tail_local].nbs[tail_slot].is_some();
                            if !head_used && !tail_used {
                                shard[head_local].nbs[head_slot] = Some(Link {
                                    id: tail.id,
                                    overlap: overlap as u8,
                                });
                                shard[tail_local].nbs[tail_slot] = Some(Link {
                                    id: head.id,
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
                                    &mut lower[head_index % shard_len],
                                    &mut upper[tail_index % shard_len],
                                )
                            } else {
                                (
                                    &mut upper[head_index % shard_len],
                                    &mut lower[tail_index % shard_len],
                                )
                            };
                            let head_used = head_end.nbs[head_slot].is_some();
                            let tail_used = tail_end.nbs[tail_slot].is_some();
                            if !head_used && !tail_used {
                                head_end.nbs[head_slot] = Some(Link {
                                    id: tail.id,
                                    overlap: overlap as u8,
                                });
                                tail_end.nbs[tail_slot] = Some(Link {
                                    id: head.id,
                                    overlap: overlap as u8,
                                });
                            }
                            (head_used, tail_used)
                        };
                        if !head_used && !tail_used {
                            head_slice[i].id = u32::MAX;
                            tail_slice[j].id = u32::MAX;
                            marked_heads += 1;
                            marked_tails += 1;
                            merged += 1;
                            i += 1;
                            j += 1;
                        } else {
                            if head_used {
                                head_slice[i].id = u32::MAX;
                                marked_heads += 1;
                                i += 1;
                            }
                            if tail_used {
                                tail_slice[j].id = u32::MAX;
                                marked_tails += 1;
                                j += 1;
                            }
                        }
                    }
                    (merged, marked_heads, marked_tails)
                })
                .reduce(
                    || (0, 0, 0),
                    |(m, h, t), (dm, dh, dt)| (m + dm, h + dh, t + dt),
                )
        };
        drop(shards);
        marked_heads += newly_marked_heads;
        marked_tails += newly_marked_tails;
        for (entries, marked) in [
            (&mut heads, &mut marked_heads),
            (&mut tails, &mut marked_tails),
        ] {
            // Only filter out dead entries if that's more than half of them.
            if *marked > entries.len() / 2 {
                info!("Drop {} of {} entries..", *marked, entries.len());
                let chunk_len = entries.len().div_ceil(chunks).max(1);
                let removed = entries
                    .par_chunks(chunk_len)
                    .map(|slice| slice.iter().filter(|entry| entry.id == u32::MAX).count())
                    .collect::<Vec<_>>();
                let removed_total = removed.iter().sum::<usize>();
                debug_assert_eq!(removed_total, *marked);
                let retained = entries.len() - removed_total;
                scratch.resize(
                    retained,
                    HeadOrTail {
                        key: K::ZERO,
                        id: 0,
                    },
                );
                rayon::scope(|scope| {
                    let mut output = scratch.as_mut_slice();
                    for (input, removed) in entries.chunks(chunk_len).zip(removed) {
                        let (part, rest) = output.split_at_mut(input.len() - removed);
                        output = rest;
                        scope.spawn(move |_| {
                            let mut write = 0;
                            for &entry in input {
                                if entry.id != u32::MAX {
                                    part[write] = entry;
                                    write += 1;
                                }
                            }
                            debug_assert_eq!(write, part.len());
                        });
                    }
                });
                std::mem::swap(entries, &mut scratch);
                info!("Final size: {} entries", entries.len());
                *marked = 0;
            }
        }
        total_merged += merged;
        total_len -= merged * overlap;
        info!(
            "overlap {overlap} merged {:>9} total merged {total_merged:>9} remaining {:>9} total len {total_len:>11}",
            merged,
            links.len() - total_merged
        );
    }
    assert_eq!(total_merged, links.len());
    // TODO break cycles

    unimplemented!("break cycles");

    let mut num_cycles = 0;
    let mut done = vec![false; links.len()];
    for i in 0..links.len() {
        if done[i] {
            continue;
        }
        let mut j = i;
        let mut reverse = false;
        loop {
            done[j] = true;
            let link = links[j].nbs[usize::from(!reverse)].as_ref().unwrap();
            j = (link.id >> 1) as usize;
            reverse = link.id & 1 != 0;
            if j == i {
                break;
            }
        }
        num_cycles += 1;
    }
    eprintln!("Number of cycles: {num_cycles}");

    todo!()
}
