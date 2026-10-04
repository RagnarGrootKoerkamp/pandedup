//! Greedy matching of unitig ends through the compacted de Bruijn graph.
//!
//! Each unitig has two oriented copies. Their outgoing tails and incoming
//! heads are grouped by their `(k-1)`-mers to build tail-to-head graph edges.
//! Following one of these edges joins two unitigs with a direct `(k-1)`-base
//! overlap. Traversing an oriented unitig from head to tail costs its length
//! minus `k-1`: the number of new bases contributed by that unitig.
//!
//! First, greedily link free ends across direct graph edges (distance zero).
//! For longer overlaps, `DistanceField` stores the shortest distance from a
//! free outgoing tail to each oriented tail. A path to a free receiving head
//! can be read from the same field by reversing its orientation. Thus a path
//! meeting inside oriented unitig `v` has cost
//!
//! distance to v's head + traversal cost of v + distance from v's tail.
//!
//! Distances are stored only below `k/2`. Every path with positive overlap
//! has a meeting unitig for which both sides fit within that bound. The
//! resulting cost `d` gives an overlap of `k-1-d`; costs at least `k-1` give
//! no overlap and are left for the final arbitrary pairing.
//!
//! Each traversal with a possible positive-overlap bridge initially gets one
//! entry in the bucket for its cheapest cost. Buckets are processed in
//! increasing cost. When a bridge links two free ends, those ends cease to
//! be distance-zero sources.
//! `remove_source` repairs affected distances on both orientations. Deleting
//! sources can only increase bridge costs, so each queued traversal is
//! rechecked when popped and moved to a later bucket if necessary. A bridge
//! may return to its own physical end; this closes that end of an output path.
//! Finally, remaining ends are paired with zero overlap, and the links are
//! written as a masked superstring.

use crate::{log_file_stats, mss::MssKey, timing::StageTiming};
use packed_seq::{PackedSeqVec, SeqVec, complement_char};
use std::fmt::{self, Display, Formatter};
use std::io::{BufWriter, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};
use tracing::info;
use voracious_radix_sort::{RadixSort, Radixable};

const DEAD: u32 = u32::MAX;

/// Decimal SI units with three significant digits for counts of at least 1k.
struct Compact(u128);

fn compact(value: impl TryInto<u64>) -> Compact {
    Compact(value.try_into().ok().expect("count exceeds u64") as u128)
}

impl Display for Compact {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        const SUFFIXES: [&str; 7] = ["", "k", "M", "G", "T", "P", "E"];
        let value = self.0;
        let mut unit = 0;
        let mut divisor = 1u128;
        while unit + 1 < SUFFIXES.len() && value >= divisor * 1_000 {
            unit += 1;
            divisor *= 1_000;
        }
        if unit == 0 {
            return write!(f, "{value}");
        }
        loop {
            let whole = value / divisor;
            let mut decimals = if whole < 10 {
                2
            } else if whole < 100 {
                1
            } else {
                0
            };
            loop {
                let factor = 10u128.pow(decimals);
                let rounded = (value * factor + divisor / 2) / divisor;
                if rounded >= 1_000 {
                    if decimals > 0 {
                        decimals -= 1;
                        continue;
                    }
                    if unit + 1 < SUFFIXES.len() {
                        unit += 1;
                        divisor *= 1_000;
                        break;
                    }
                }
                return match decimals {
                    2 => write!(
                        f,
                        "{}.{:02}{}",
                        rounded / 100,
                        rounded % 100,
                        SUFFIXES[unit]
                    ),
                    1 => write!(f, "{}.{}{}", rounded / 10, rounded % 10, SUFFIXES[unit]),
                    _ => write!(f, "{}{}", rounded, SUFFIXES[unit]),
                };
            }
        }
    }
}

#[derive(Clone, Copy, Default)]
struct VisitStats {
    direct_overlap: u64,
    meet: u64,
    connect: u64,
    remove_source: u64,
}

#[derive(Clone, Copy)]
#[repr(C, packed(1))]
#[doc(hidden)]
pub struct End<K> {
    kmer: K,
    id: u32,
}

const _: () = {
    assert!(std::mem::size_of::<End<u64>>() == 12);
    assert!(std::mem::size_of::<End<u128>>() == 20);
};

impl Radixable<u64> for End<u64> {
    type Key = u64;
    fn key(&self) -> Self::Key {
        self.kmer
    }
}

impl Radixable<u128> for End<u128> {
    type Key = u128;
    fn key(&self) -> Self::Key {
        self.kmer
    }
}

impl<K: MssKey> PartialEq for End<K> {
    fn eq(&self, other: &Self) -> bool {
        let left = self.kmer;
        let right = other.kmer;
        left == right
    }
}

impl<K: MssKey> PartialOrd for End<K> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let left = self.kmer;
        let right = other.kmer;
        Some(left.cmp(&right))
    }
}

#[derive(Clone, Copy, PartialEq, PartialOrd)]
struct Edge {
    from: u32,
    to: u32,
}

impl Radixable<u64> for Edge {
    type Key = u64;
    fn key(&self) -> Self::Key {
        (u64::from(self.from) << 32) | u64::from(self.to)
    }
}

/// CSR adjacency indexed by an oriented outgoing tail. An edge's destination
/// is the oriented unitig entered at its head; its low bit is the orientation.
struct Graph {
    edges: Vec<u32>,
    offsets: Vec<u32>,
}

impl Graph {
    fn outgoing(&self, id: u32) -> &[u32] {
        &self.edges[self.offsets[id as usize] as usize..self.offsets[id as usize + 1] as usize]
    }

    /// Incoming neighbors, recovered from reverse-complement graph edges.
    fn incoming(&self, id: u32) -> impl ExactSizeIterator<Item = u32> + '_ {
        self.outgoing(id ^ 1).iter().map(|&reverse| reverse ^ 1)
    }
}

fn graph<K: MssKey>(k: usize, seq: &PackedSeqVec, ranges: &[Range<usize>]) -> Graph
where
    End<K>: Radixable<K, Key = K>,
{
    let overlap = k - 1;
    let mut heads = Vec::with_capacity(2 * ranges.len());
    let mut tails = Vec::with_capacity(2 * ranges.len());
    info!("Building graph of {} unitigs..", compact(ranges.len()));
    for (index, range) in ranges.iter().enumerate() {
        let id = (index as u32) * 2;
        let (head, tail, rc_head, rc_tail): (K, K, K, K) = (
            K::read_kmer(seq, overlap, range.start),
            K::read_kmer(seq, overlap, range.end - overlap),
            K::read_revcomp_kmer(seq, overlap, range.start),
            K::read_revcomp_kmer(seq, overlap, range.end - overlap),
        );
        heads.push(End { kmer: head, id });
        heads.push(End {
            kmer: rc_tail,
            id: id + 1,
        });
        tails.push(End { kmer: tail, id });
        tails.push(End {
            kmer: rc_head,
            id: id + 1,
        });
    }
    info!("Sorting heads");
    heads.voracious_mt_sort(rayon::current_num_threads());
    info!("Sorting tails");
    tails.voracious_mt_sort(rayon::current_num_threads());

    info!("Matching unitig ends");
    let mut pairs = Vec::new();
    let (mut h, mut t) = (0, 0);
    while h < heads.len() && t < tails.len() {
        // Copy keys out of packed entries before borrowing them for comparison.
        let head_key = heads[h].kmer;
        let tail_key = tails[t].kmer;
        match head_key.cmp(&tail_key) {
            std::cmp::Ordering::Less => h += 1,
            std::cmp::Ordering::Greater => t += 1,
            std::cmp::Ordering::Equal => {
                let key = heads[h].kmer;
                let mut he = h;
                let mut te = t;
                while he < heads.len() {
                    let next = heads[he].kmer;
                    if next != key {
                        break;
                    }
                    he += 1;
                }
                while te < tails.len() {
                    let next = tails[te].kmer;
                    if next != key {
                        break;
                    }
                    te += 1;
                }
                for tail in &tails[t..te] {
                    for head in &heads[h..he] {
                        // The reverse orientation of the same endpoint does not
                        // traverse a unitig or introduce a new graph edge.
                        if tail.id != (head.id ^ 1) {
                            pairs.push(Edge {
                                from: tail.id,
                                to: head.id,
                            });
                        }
                    }
                }
                h = he;
                t = te;
            }
        }
    }
    info!("Sorting {} edges", compact(pairs.len()));
    pairs.voracious_mt_sort(rayon::current_num_threads());
    // let old_len = pairs.len();
    // info!("Dedup edges");
    // pairs.dedup();
    // info!("Dedup edges => {} left", pairs.len());
    // assert_eq!(old_len, pairs.len(), "duplicate graph edges");
    assert!(pairs.len() < u32::MAX as usize, "too many graph edges");
    info!("Building CSR adjacency");
    let mut offsets = vec![0; 2 * ranges.len() + 1];
    info!(
        "Filling offsets: {}B",
        compact(std::mem::size_of_val(offsets.as_slice()))
    );
    for edge in &pairs {
        offsets[edge.from as usize + 1] += 1;
    }
    for i in 1..offsets.len() {
        offsets[i] += offsets[i - 1];
    }
    info!("Mapping pairs");
    Graph {
        edges: pairs.into_iter().map(|edge| edge.to).collect(),
        offsets,
    }
}

#[derive(Clone, Copy)]
#[repr(C, packed(1))]
struct Link {
    target: u32,
    overlap: u8,
}

const _: () = assert!(std::mem::size_of::<Link>() == 5);

impl Link {
    const EMPTY: Self = Self {
        target: DEAD,
        overlap: 0,
    };

    fn is_empty(self) -> bool {
        self.target == DEAD
    }
}

const UNREACHABLE: u8 = u8::MAX;

/// The cost of crossing unitig `node`.
fn traversal_weight(ranges: &[Range<usize>], node: u32, k: usize) -> usize {
    ranges[(node / 2) as usize].len() - (k - 1)
}

/// Distance to the nearest still-free outgoing end within a bounded radius.
struct DistanceField {
    /// Distances to all unitig ends.
    distance: Vec<u8>,
    /// Scratch buffer `remove_source`.
    worklist: Vec<(u32, u8)>,
}

impl DistanceField {
    fn new(nodes: usize) -> Self {
        Self {
            distance: vec![0; nodes],
            worklist: Vec::new(),
        }
    }

    /// Compute the distance to a node by iterating over its predecessors.
    fn replacement_distance(
        &self,
        node: u32,
        graph: &Graph,
        ranges: &[Range<usize>],
        k: usize,
        radius: usize,
        visited: &mut u64,
    ) -> u8 {
        let weight = traversal_weight(ranges, node, k);
        let predecessors = graph.incoming(node);
        *visited += 1 + predecessors.len() as u64;
        let mut best = radius;
        for predecessor in predecessors {
            best = best.min(self.distance[predecessor as usize] as usize + weight);
        }
        if best == radius {
            UNREACHABLE
        } else {
            best as u8
        }
    }

    /// Propagate increases after deleting a source. A neighbour only needs
    /// repair if its current label was tight through the changed node.
    fn remove_source(
        &mut self,
        source: u32,
        graph: &Graph,
        ranges: &[Range<usize>],
        k: usize,
        radius: usize,
        visited: &mut u64,
    ) {
        assert_eq!(
            self.distance[source as usize], 0,
            "active source has nonzero distance"
        );
        self.worklist.clear();
        self.distance[source as usize] =
            self.replacement_distance(source, graph, ranges, k, radius, visited);
        debug_assert!(self.distance[source as usize] > 0);
        self.worklist.push((source, 0));
        while let Some((node, old_distance)) = self.worklist.pop() {
            let next_nodes = graph.outgoing(node);
            *visited += 1 + next_nodes.len() as u64;
            for &next in next_nodes {
                let old_next = self.distance[next as usize];
                if old_next != UNREACHABLE
                    && old_distance as usize + traversal_weight(ranges, next, k)
                        == old_next as usize
                {
                    let new_next =
                        self.replacement_distance(next, graph, ranges, k, radius, visited);
                    debug_assert!(new_next >= old_next);
                    if new_next > old_next {
                        self.distance[next as usize] = new_next;
                        self.worklist.push((next, old_next));
                    }
                }
            }
        }
    }

    /// Follow tight incoming edges to a free source. Positive traversal
    /// weights make the distance decrease at every step.
    ///
    /// Note that now this is just a linear scan to follow a path backwards, not a DFS anymore.
    fn source(
        &self,
        mut node: u32,
        graph: &Graph,
        ranges: &[Range<usize>],
        k: usize,
        visited: &mut u64,
    ) -> u32 {
        loop {
            *visited += 1;
            let distance = self.distance[node as usize];
            assert_ne!(
                distance, UNREACHABLE,
                "source requested for unreachable node"
            );
            if distance == 0 {
                return node;
            }
            let weight = traversal_weight(ranges, node, k);
            assert!(weight <= distance as usize, "invalid distance label");
            let mut predecessors = graph.incoming(node);
            *visited += predecessors.len() as u64;
            node = predecessors
                .find(|&predecessor| {
                    let prior = self.distance[predecessor as usize];
                    prior != UNREACHABLE && prior as usize + weight == distance as usize
                })
                .expect("finite distance has no tight predecessor");
        }
    }
}

fn tail_slot(id: u32) -> usize {
    (id ^ 1) as usize
}
fn head_slot(id: u32) -> usize {
    id as usize
}

/// Connect unmatched unitigs ends `from` to `to`.
fn connect(
    links: &mut [Link],
    from: u32,
    to: u32,
    overlap: usize,
    total_bases: &mut usize,
    visited: &mut u64,
) {
    debug_assert!(links[tail_slot(from)].is_empty());
    debug_assert!(links[head_slot(to)].is_empty());
    debug_assert!(overlap <= u8::MAX as usize);
    links[tail_slot(from)] = Link {
        target: to,
        overlap: overlap as u8,
    };
    links[head_slot(to)] = Link {
        target: from ^ 1,
        overlap: overlap as u8,
    };
    *total_bases -= overlap;
    *visited += if tail_slot(from) == head_slot(to) {
        1
    } else {
        2
    };
}

/// A path through `node` costs the distance to its head, its own traversal
/// weight, and the distance from its tail to a free receiving head. The latter
/// is read from the same field through reverse-complement orientation.
fn meet_at(
    node: u32,
    wanted: usize,
    field: &DistanceField,
    graph: &Graph,
    ranges: &[Range<usize>],
    k: usize,
    visited: &mut u64,
) -> Option<(u32, u32)> {
    let weight = traversal_weight(ranges, node, k);
    if weight > wanted {
        return None;
    }
    let mut successor_at = [DEAD; 64];
    let successors = graph.outgoing(node);
    *visited += 1 + successors.len() as u64;
    for &successor in successors {
        let distance = field.distance[(successor ^ 1) as usize];
        if distance != UNREACHABLE && successor_at[distance as usize] == DEAD {
            successor_at[distance as usize] = successor;
        }
    }
    let predecessors = graph.incoming(node);
    *visited += 1 + predecessors.len() as u64;
    for predecessor in predecessors {
        let left = field.distance[predecessor as usize] as usize;
        if left > wanted - weight {
            continue;
        }
        let right = wanted - weight - left;
        if right >= successor_at.len() || successor_at[right] == DEAD {
            continue;
        }
        let source = field.source(predecessor, graph, ranges, k, visited);
        let receiver_end = field.source(successor_at[right] ^ 1, graph, ranges, k, visited);
        return Some((source, receiver_end ^ 1));
    }
    None
}

/// Cheapest bridge through one oriented unitig at or above `minimum`.
/// A bucket owns its entry until it is popped, so deletions only require
/// rechecking and possibly moving that entry to a later bucket.
fn bridge_cost(
    node: u32,
    minimum: usize,
    field: &DistanceField,
    graph: &Graph,
    ranges: &[Range<usize>],
    k: usize,
    visited: &mut u64,
) -> Option<usize> {
    let weight = traversal_weight(ranges, node, k);
    if weight >= k - 1 {
        return None;
    }
    let mut right_mask = 0u64;
    let successors = graph.outgoing(node);
    *visited += 1 + successors.len() as u64;
    for &successor in successors {
        let right = field.distance[(successor ^ 1) as usize];
        if right != UNREACHABLE {
            right_mask |= 1u64 << right;
        }
    }
    if right_mask == 0 {
        return None;
    }
    let predecessors = graph.incoming(node);
    *visited += 1 + predecessors.len() as u64;
    let mut best = k - 1;
    for predecessor in predecessors {
        let left = field.distance[predecessor as usize];
        if left == UNREACHABLE {
            continue;
        }
        let base = left as usize + weight;
        if base >= best {
            continue;
        }
        let required_right = minimum.saturating_sub(base);
        if required_right >= 64 {
            continue;
        }
        let eligible = right_mask >> required_right;
        if eligible != 0 {
            best = best.min(base + required_right + eligible.trailing_zeros() as usize);
        }
    }
    (best < k - 1).then_some(best)
}

/// Greedy distance-ordered matching through central unitig traversals.
fn match_ends(
    k: usize,
    ranges: &[Range<usize>],
    graph: &Graph,
    initial_bases: usize,
) -> (Vec<Link>, usize) {
    let nodes = 2 * ranges.len();
    let radius = k / 2;
    let mut links = vec![Link::EMPTY; nodes];
    let mut field = DistanceField::new(nodes);
    let mut free_ends = nodes;
    let mut visits = VisitStats::default();
    let mut estimated_bases = initial_bases;
    info!(
        "Matching {} unitigs, {} ends, k {}: two ends per unitig; overlap = k minus one minus distance; examined = bridge entries checked; links made = connections; a self-link occupies one end, other links occupy two",
        compact(ranges.len()),
        compact(nodes),
        compact(k),
    );
    info!(
        "Meeting search: source and receiver fronts each store distances below {}; a path through one unitig has cost left distance + unitig extensions + right distance",
        compact(radius),
    );
    info!(
        "Before distance {}: {} receiving ends examined, {} links made, {} receiving ends remain, {} free ends, {} estimated bases",
        compact(0u64),
        compact(0u64),
        compact(0u64),
        compact(nodes),
        compact(nodes),
        compact(initial_bases),
    );
    info!(
        "Node visits count repeated inspections: direct_overlap scans graph neighbors at distance zero; meet checks and schedules central unitig traversals and traces endpoints; connect touches link slots; remove_source repairs both bounded fronts"
    );

    // An edge in the unitig graph is a direct (k-1)-character overlap. Match
    // those edges without querying shortest-path labels.
    // TODO / HOT: This is slow
    let mut direct_examined = 0usize;
    let mut direct_made = 0usize;
    for target in 0..nodes as u32 {
        if k == 1 {
            break;
        }
        if !links[head_slot(target)].is_empty() {
            continue;
        }
        direct_examined += 1;
        for source in graph.incoming(target) {
            visits.direct_overlap += 1;
            if !links[tail_slot(source)].is_empty() {
                continue;
            }
            connect(
                &mut links,
                source,
                target,
                k - 1,
                &mut estimated_bases,
                &mut visits.connect,
            );
            free_ends -= if source == target ^ 1 { 1 } else { 2 };
            direct_made += 1;
            field.remove_source(source, graph, ranges, k, radius, &mut visits.remove_source);
            if source != target ^ 1 {
                field.remove_source(
                    target ^ 1,
                    graph,
                    ranges,
                    k,
                    radius,
                    &mut visits.remove_source,
                );
            }
            break;
        }
    }
    info!(
        "After distance {}: {} receiving ends examined, {} links made, {} receiving ends remain, {} free ends, {} estimated bases",
        compact(0u64),
        compact(direct_examined),
        compact(direct_made),
        compact(free_ends),
        compact(free_ends),
        compact(estimated_bases),
    );
    info!(
        "After distance {} node visits: direct_overlap {} (total {}), meet {} (total {}), connect {} (total {}), remove_source {} (total {})",
        compact(0u64),
        compact(visits.direct_overlap),
        compact(visits.direct_overlap),
        compact(0u64),
        compact(0u64),
        compact(visits.connect),
        compact(visits.connect),
        compact(visits.remove_source),
        compact(visits.remove_source),
    );

    // Every oriented end starts as a distance-zero source. After the direct
    // links, remove_source has already repaired the bounded distance field.
    // Seed each internal traversal once; later deletions only increase its
    // bridge cost, so we recheck and move its single entry when popped.
    let mut buckets = vec![Vec::<u32>::new(); k];
    for node in 0..nodes as u32 {
        if let Some(cost) = bridge_cost(node, 1, &field, graph, ranges, k, &mut visits.meet) {
            buckets[cost].push(node);
        }
    }
    for distance in 1..k - 1 {
        let before = visits;
        let mut examined = 0usize;
        let mut made = 0usize;
        while let Some(node) = buckets[distance].pop() {
            examined += 1;
            let Some(cost) =
                bridge_cost(node, distance, &field, graph, ranges, k, &mut visits.meet)
            else {
                continue;
            };
            if cost > distance {
                buckets[cost].push(node);
                continue;
            }
            let (source, target) =
                meet_at(node, distance, &field, graph, ranges, k, &mut visits.meet)
                    .expect("bridge cost has no meeting path");
            // The same traversal may connect another pair at this distance.
            buckets[distance].push(node);
            connect(
                &mut links,
                source,
                target,
                k - 1 - distance,
                &mut estimated_bases,
                &mut visits.connect,
            );
            free_ends -= if source == target ^ 1 { 1 } else { 2 };
            made += 1;
            field.remove_source(source, graph, ranges, k, radius, &mut visits.remove_source);
            if source != target ^ 1 {
                field.remove_source(
                    target ^ 1,
                    graph,
                    ranges,
                    k,
                    radius,
                    &mut visits.remove_source,
                );
            }
        }
        info!(
            "After distance {}: {} bridge entries examined, {} links made, {} receiving ends remain, {} free ends, {} estimated bases",
            compact(distance),
            compact(examined),
            compact(made),
            compact(free_ends),
            compact(free_ends),
            compact(estimated_bases),
        );
        info!(
            "After distance {} node visits: direct_overlap {} (total {}), meet {} (total {}), connect {} (total {}), remove_source {} (total {})",
            compact(distance),
            compact(0u64),
            compact(visits.direct_overlap),
            compact(visits.meet - before.meet),
            compact(visits.meet),
            compact(visits.connect - before.connect),
            compact(visits.connect),
            compact(visits.remove_source - before.remove_source),
            compact(visits.remove_source),
        );
    }

    // Link up (concatenate) remaining seqs.
    let outgoing: Vec<_> = (0..nodes as u32)
        .filter(|&id| links[tail_slot(id)].is_empty())
        .collect();
    let mut pairs = outgoing.chunks_exact(2);
    for pair in &mut pairs {
        let target = pair[1] ^ 1;
        connect(
            &mut links,
            pair[0],
            target,
            0,
            &mut estimated_bases,
            &mut visits.connect,
        );
    }
    // A self-link without overlap is fine.
    if let [source] = pairs.remainder() {
        connect(
            &mut links,
            *source,
            *source ^ 1,
            0,
            &mut estimated_bases,
            &mut visits.connect,
        );
    }
    info!(
        "Final node visits (including closure links): direct_overlap {}, meet {}, connect {}, remove_source {}",
        compact(visits.direct_overlap),
        compact(visits.meet),
        compact(visits.connect),
        compact(visits.remove_source),
    );
    (links, estimated_bases)
}

fn append(
    output: &mut Vec<u8>,
    seq: &PackedSeqVec,
    range: &Range<usize>,
    reverse: bool,
    overlap: usize,
    k: usize,
) {
    let mut bases = seq.slice(range.clone()).unpack();
    if reverse {
        bases.reverse();
        for base in &mut bases {
            *base = complement_char(*base);
        }
    }
    if !output.is_empty() {
        let mask = (k - 1 - overlap).min(bases.len() - overlap);
        for base in &mut bases[overlap..overlap + mask] {
            *base = base.to_ascii_lowercase();
        }
    }
    output.extend_from_slice(&bases[overlap..]);
}

/// Build a masked superstring from unitigs. Newly created crossing k-mers are
/// lowercase, as in `mss::masked_superstring`.
pub fn masked_superstring<K: MssKey>(
    k: usize,
    seq: PackedSeqVec,
    ranges: Vec<Range<usize>>,
) -> Vec<u8>
where
    End<K>: Radixable<K, Key = K>,
{
    assert!(k > 0 && k <= K::BITS / 2);
    let initial_bases: usize = ranges.iter().map(Range::len).sum();
    let (ranges, short): (Vec<_>, Vec<_>) = ranges.into_iter().partition(|r| r.len() >= k);
    assert!(ranges.len() <= (u32::MAX as usize / 2), "too many unitigs");
    let graph = graph::<K>(k, &seq, &ranges);

    let (links, estimated_bases) = match_ends(k, &ranges, &graph, initial_bases);

    info!("Reconstruct output");
    let mut output = Vec::new();
    let mut done = vec![false; ranges.len()];
    // Start open paths at their self-linked head, oriented away from it.
    // Remaining components have no such endpoint and can be traversed as cycles.
    for pass in 0..2 {
        for start in 0..ranges.len() {
            if done[start] {
                continue;
            }
            let forward = (start as u32) * 2;
            let reverse = forward ^ 1;
            let mut id = if links[head_slot(forward)].target == forward {
                forward
            } else if links[head_slot(reverse)].target == reverse {
                reverse
            } else if pass == 0 {
                continue;
            } else {
                forward
            };
            let mut overlap = 0;
            loop {
                let index = (id / 2) as usize;
                assert!(!done[index], "link cycle revisits unitig before closing");
                done[index] = true;
                append(&mut output, &seq, &ranges[index], id & 1 != 0, overlap, k);
                let link = links[tail_slot(id)];
                assert!(!link.is_empty(), "unmatched unitig end");
                if link.target == id ^ 1 || (link.target / 2) as usize == start {
                    break;
                }
                id = link.target;
                overlap = link.overlap as usize;
            }
        }
    }
    for range in &short {
        if !range.is_empty() {
            append(&mut output, &seq, range, false, 0, k);
        }
    }
    info!(
        "Length summary: {} estimated bases after links, {} output bases",
        compact(estimated_bases),
        compact(output.len()),
    );
    output
}

/// Read unitigs and write one masked superstring in FASTA format.
pub fn run(input: &Path, output: Option<&Path>, k: usize) -> PathBuf {
    let timing = StageTiming::start();
    info!("Reading unitigs..");
    let (seq, ranges) = PackedSeqVec::from_fastx(input);
    let input_bases = ranges.iter().map(Range::len).sum();
    log_file_stats("Read", input, Some((ranges.len(), input_bases))).unwrap();
    let superstring = if k <= 32 {
        masked_superstring::<u64>(k, seq, ranges)
    } else {
        masked_superstring::<u128>(k, seq, ranges)
    };
    let output = output
        .map(Path::to_path_buf)
        .unwrap_or_else(|| input.with_extension("matchtigs.msfa"));
    let mut writer = BufWriter::new(std::fs::File::create(&output).unwrap());
    writeln!(writer, ">matchtigs-masked-superstring").unwrap();
    writer.write_all(&superstring).unwrap();
    writer.write_all(b"\n").unwrap();
    drop(writer);
    log_file_stats("Wrote", &output, Some((1, superstring.len()))).unwrap();
    info!(
        "matchtigs: {}, output {} bases ({})",
        timing.finish(),
        compact(superstring.len()),
        output.display()
    );
    output
}
