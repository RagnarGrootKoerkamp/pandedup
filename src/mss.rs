//! Compute a masked superstring using kmercamel's bigreedy algorithm.
//!
//! See: https://github.com/OndrejSladky/kmercamel, and the paper:
//!
//! Ondřej Sladký, Pavel Veselý, and Karel Břinda:
//! Masked superstrings as a unified framework for textual k-mer set representations.
//! bioRxiv 2023.02.01.526717, 2023. https://doi.org/10.1101/2023.02.01.526717

use crate::{log_file_stats, timing::StageTiming};
use core::ops::{BitAnd, Mul, Range, Shl, Shr, Sub};
use packed_seq::{PackedSeqVec, SeqVec, complement_char};
use rayon::prelude::*;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use tracing::{debug, info, trace};
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

pub trait MssIndex: Copy + Eq + Send + Sync + Default {
    const MAX: Self;
    const MAX_CONTIGS: usize;
    fn from_parts(index: usize, reverse: bool) -> Self;
    fn contig_index(self) -> usize;
    fn reverse(self) -> bool;
    fn flip(self) -> Self;
}

impl MssIndex for u32 {
    const MAX: Self = u32::MAX;
    const MAX_CONTIGS: usize = (u32::MAX >> 1) as usize;
    fn from_parts(index: usize, reverse: bool) -> Self {
        ((index as u32) << 1) | u32::from(reverse)
    }
    fn contig_index(self) -> usize {
        (self >> 1) as usize
    }
    fn reverse(self) -> bool {
        self & 1 != 0
    }
    fn flip(self) -> Self {
        self ^ 1
    }
}

impl MssIndex for u64 {
    const MAX: Self = u64::MAX;
    const MAX_CONTIGS: usize = (u64::MAX >> 1) as usize;
    fn from_parts(index: usize, reverse: bool) -> Self {
        ((index as u64) << 1) | u64::from(reverse)
    }
    fn contig_index(self) -> usize {
        (self >> 1) as usize
    }
    fn reverse(self) -> bool {
        self & 1 != 0
    }
    fn flip(self) -> Self {
        self ^ 1
    }
}

#[derive(Clone, Copy)]
#[repr(C, packed(1))]
pub struct HeadOrTail<K: MssKey, I: MssIndex> {
    key: K,
    /// Twice the contig index, plus one for reverse complement; MAX means used.
    id: I,
}

const _: () = {
    assert!(std::mem::align_of::<HeadOrTail<u64, u32>>() == 1);
    assert!(std::mem::size_of::<HeadOrTail<u64, u32>>() == 12);
    assert!(std::mem::size_of::<HeadOrTail<u64, u64>>() == 16);
    assert!(std::mem::align_of::<HeadOrTail<u128, u64>>() == 1);
    assert!(std::mem::size_of::<HeadOrTail<u128, u32>>() == 20);
    assert!(std::mem::size_of::<HeadOrTail<u128, u64>>() == 24);
};

impl<K: MssKey, I: MssIndex> PartialEq for HeadOrTail<K, I> {
    fn eq(&self, other: &Self) -> bool {
        let key = self.key;
        let other_key = other.key;
        key == other_key
    }
}

impl<K: MssKey, I: MssIndex> PartialOrd for HeadOrTail<K, I> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let key = self.key;
        let other_key = other.key;
        Some(key.cmp(&other_key))
    }
}

impl<I: MssIndex> Radixable<u64> for HeadOrTail<u64, I> {
    type Key = u64;
    fn key(&self) -> Self::Key {
        self.key
    }
}

impl<I: MssIndex> Radixable<u128> for HeadOrTail<u128, I> {
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
fn merge_head_part<K: MssKey, I: MssIndex>(
    input: &[HeadOrTail<K, I>],
    output: &mut [HeadOrTail<K, I>],
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
fn resort_heads<K: MssKey, I: MssIndex>(
    heads: &mut Vec<HeadOrTail<K, I>>,
    scratch: &mut Vec<HeadOrTail<K, I>>,
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
            id: I::default(),
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
struct Links<I: MssIndex> {
    /// Neighbors at the forward head (0) and forward tail (1).
    nbs: [Option<Link<I>>; 2],
}

#[derive(Clone, Copy)]
#[repr(C, packed(1))]
struct Link<I: MssIndex> {
    /// Oriented contig to traverse after crossing this link.
    id: I,
    overlap: u8,
}

fn append_contig(
    output: &mut Vec<u8>,
    seq: &PackedSeqVec,
    range: &Range<usize>,
    reverse: bool,
    overlap: usize,
    k: usize,
) {
    // FIXME: Use packed_seq utils for this somehow
    let mut bases = seq.slice(range.clone()).unpack();
    if reverse {
        bases.reverse();
        for base in &mut bases {
            *base = complement_char(*base);
        }
    }
    assert!(overlap < k && overlap <= bases.len() && overlap <= output.len());
    if !output.is_empty() {
        let end = output.len();
        debug_assert!(
            output[end - overlap..]
                .iter()
                .zip(&bases[..overlap])
                .all(|(left, right)| left.to_ascii_uppercase() == *right)
        );
        // The first k-1-overlap appended bases end the new crossing k-mers.
        // Leave the existing output, including the cycle's first contig, intact.
        let mask_len = (k - 1 - overlap).min(bases.len() - overlap);
        for base in &mut bases[overlap..overlap + mask_len] {
            *base = base.to_ascii_lowercase();
        }
    }
    output.extend_from_slice(&bases[overlap..]);
}

fn reconstruct_cycles<I: MssIndex>(
    k: usize,
    seq: &PackedSeqVec,
    ranges: &[Range<usize>],
    short_ranges: &[Range<usize>],
    links: &[Links<I>],
    capacity: usize,
) -> Vec<u8> {
    let mut output = Vec::with_capacity(capacity);
    let mut done = vec![false; links.len()];
    let mut num_cycles = 0;
    for start in 0..links.len() {
        if done[start] {
            continue;
        }
        let mut id = I::from_parts(start, false);
        let mut overlap = 0;
        loop {
            let index = id.contig_index();
            assert!(!done[index], "cycle revisits contig {index} before closing");
            done[index] = true;
            append_contig(&mut output, seq, &ranges[index], id.reverse(), overlap, k);
            let link =
                links[index].nbs[usize::from(!id.reverse())].expect("missing link at contig tail");
            let next = link.id;
            if next.contig_index() == start {
                break;
            }
            id = next;
            overlap = link.overlap as usize;
        }
        num_cycles += 1;
    }
    // Short contigs contain no input k-mers, but retain their sequences.
    for range in short_ranges {
        if !range.is_empty() {
            append_contig(&mut output, seq, range, false, 0, k);
        }
    }
    info!(
        "Number of cycles: {num_cycles}; output length: {}",
        output.len()
    );
    output
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
    ranges: Vec<Range<usize>>,
) -> Vec<u8>
where
    HeadOrTail<K, u32>: Radixable<K, Key = K>,
    HeadOrTail<K, u64>: Radixable<K, Key = K>,
{
    info!("masked_superstring: k={}, contigs={}", k, ranges.len());
    assert!(k > 0 && k <= K::BITS / 2);
    let total_len = ranges.iter().map(|c| c.len()).sum::<usize>();
    // Filter first so indices in both entry vectors match `ends` even when
    // the input contains contigs shorter than k.
    let (ranges, short_ranges): (Vec<_>, Vec<_>) = ranges.into_iter().partition(|c| c.len() >= k);
    if ranges.len() < (1usize << 31) {
        masked_superstring_with_index::<K, u32>(k, seq, ranges, short_ranges, total_len)
    } else {
        masked_superstring_with_index::<K, u64>(k, seq, ranges, short_ranges, total_len)
    }
}

/// Read contigs, construct a masked superstring, and write it as FASTA.
pub fn run(input: &Path, output: Option<&Path>, k: usize, threads: Option<usize>) -> PathBuf {
    let mut pool_builder = rayon::ThreadPoolBuilder::new();
    if let Some(threads) = threads {
        pool_builder = pool_builder.num_threads(threads);
    }
    let pool = pool_builder.build().unwrap();
    let timing = StageTiming::start();
    info!("Reading input..");
    let (seq, ranges) = PackedSeqVec::from_fastx(input);
    let input_bases = ranges.iter().map(|range| range.len()).sum();
    log_file_stats("Read", input, Some((ranges.len(), input_bases))).unwrap();

    let superstring = pool.install(|| {
        if k <= 32 {
            masked_superstring::<u64>(k, seq, ranges)
        } else {
            masked_superstring::<u128>(k, seq, ranges)
        }
    });
    let output = output
        .map(Path::to_path_buf)
        .unwrap_or_else(|| input.with_extension("msfa"));
    let mut writer = BufWriter::new(std::fs::File::create(&output).unwrap());
    writeln!(writer, ">masked-superstring").unwrap();
    writer.write_all(&superstring).unwrap();
    writer.write_all(b"\n").unwrap();
    drop(writer);
    log_file_stats("Wrote", &output, Some((1, superstring.len()))).unwrap();
    info!(
        "mss: {}, input {} bytes, output {} bases, {} bytes ({})",
        timing.finish(),
        std::fs::metadata(input).unwrap().len(),
        superstring.len(),
        std::fs::metadata(&output).unwrap().len(),
        output.display()
    );
    output
}

fn masked_superstring_with_index<K: MssKey, I: MssIndex>(
    k: usize,
    seq: PackedSeqVec,
    ranges: Vec<Range<usize>>,
    short_ranges: Vec<Range<usize>>,
    mut total_len: usize,
) -> Vec<u8>
where
    HeadOrTail<K, I>: Radixable<K, Key = K>,
{
    info!("Collect tig ends");
    let mut total_merged = 0;
    let threads = rayon::current_num_threads();
    let num_contigs = ranges.len();
    // Leave room for the orientation bit and reserve I::MAX as the sentinel.
    assert!(
        num_contigs <= I::MAX_CONTIGS,
        "too many contigs for packed orientation IDs"
    );

    let empty_entry = HeadOrTail {
        key: K::ZERO,
        id: I::from_parts(0, false),
    };
    let mut heads = vec![empty_entry; 2 * num_contigs];
    let mut tails = vec![empty_entry; 2 * num_contigs];
    let first_overlap = k - 1;
    let first_head_mask = (K::ONE << (2 * first_overlap)) - K::ONE;
    ranges
        .par_iter()
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
                    id: I::from_parts(index, reverse),
                };
                tail_entries[direction] = HeadOrTail {
                    key: [tail_fw, head_rc][direction],
                    id: I::from_parts(index, reverse),
                };
            }
        });

    let mut links: Vec<Links<I>> = vec![Links { nbs: [None, None] }; num_contigs];

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
        trace!("overlap: {}", overlap);
        let tail_shift = 2 * (k - overlap);
        if overlap < first_overlap {
            resort_heads(&mut heads, &mut scratch, overlap, threads);
        }

        trace!("merging..");
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
                        let head_id = head.id;
                        let tail_id = tail.id;
                        if head_id == I::MAX {
                            i += 1;
                            continue;
                        }
                        if tail_id == I::MAX {
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
                        if head_id == tail_id.flip() {
                            j += 1;
                            continue;
                        }

                        let head_index = head_id.contig_index();
                        let tail_index = tail_id.contig_index();
                        let head_slot = usize::from(head_id.reverse());
                        let tail_slot = usize::from(!tail_id.reverse());
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
                                    id: tail_id.flip(),
                                    overlap: overlap as u8,
                                });
                                shard[tail_local].nbs[tail_slot] = Some(Link {
                                    id: head_id,
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
                                    id: tail_id.flip(),
                                    overlap: overlap as u8,
                                });
                                tail_end.nbs[tail_slot] = Some(Link {
                                    id: head_id,
                                    overlap: overlap as u8,
                                });
                            }
                            (head_used, tail_used)
                        };
                        if !head_used && !tail_used {
                            head_slice[i].id = I::MAX;
                            tail_slice[j].id = I::MAX;
                            marked_heads += 1;
                            marked_tails += 1;
                            merged += 1;
                            i += 1;
                            j += 1;
                        } else {
                            if head_used {
                                head_slice[i].id = I::MAX;
                                marked_heads += 1;
                                i += 1;
                            }
                            if tail_used {
                                tail_slice[j].id = I::MAX;
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
                debug!("Drop {} of {} entries..", *marked, entries.len());
                let chunk_len = entries.len().div_ceil(chunks).max(1);
                let removed = entries
                    .par_chunks(chunk_len)
                    .map(|slice| {
                        slice
                            .iter()
                            .filter(|entry| {
                                let id = entry.id;
                                id == I::MAX
                            })
                            .count()
                    })
                    .collect::<Vec<_>>();
                let removed_total = removed.iter().sum::<usize>();
                debug_assert_eq!(removed_total, *marked);
                let retained = entries.len() - removed_total;
                scratch.resize(retained, empty_entry);
                rayon::scope(|scope| {
                    let mut output = scratch.as_mut_slice();
                    for (input, removed) in entries.chunks(chunk_len).zip(removed) {
                        let (part, rest) = output.split_at_mut(input.len() - removed);
                        output = rest;
                        scope.spawn(move |_| {
                            let mut write = 0;
                            for &entry in input {
                                let id = entry.id;
                                if id != I::MAX {
                                    part[write] = entry;
                                    write += 1;
                                }
                            }
                            debug_assert_eq!(write, part.len());
                        });
                    }
                });
                std::mem::swap(entries, &mut scratch);
                *marked = 0;
            }
        }
        total_merged += merged;
        total_len -= merged * overlap;
        debug!(
            "overlap {overlap} merged {:>9} total merged {total_merged:>9} remaining {:>9} total len {total_len:>11}",
            merged,
            links.len() - total_merged
        );
    }
    assert_eq!(total_merged, links.len());

    info!("Reconstructing cycles..");
    reconstruct_cycles(k, &seq, &ranges, &short_ranges, &links, total_len)
}
