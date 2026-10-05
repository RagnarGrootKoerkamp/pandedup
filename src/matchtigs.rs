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
//! free outgoing tail to each oriented head, before traversing that unitig.
//! A path from an oriented tail to a free receiving head can be read from the
//! same field by reversing its orientation. Thus a path meeting inside `v`
//! has cost
//!
//! distance[v] + traversal cost of v + distance[v^1].
//!
//! Distances are grown one layer at a time, up to `k/2`. At half-distance `d`,
//! scan the nodes at distance `d` and relax their outgoing edges. When both
//! orientations of a unitig have reached labels, queue its bridge by total
//! cost. Buckets at costs `2d-1` and `2d` are then matched; earlier bridges
//! wait in those buckets. A bridge of cost `c` gives
//! an overlap of `k-1-c`; costs at least `k-1` have no overlap.
//!
//! After a link consumes its free ends, `remove_source` repairs every affected
//! finite label, including tentative labels in later layers. The same
//! traversal can connect several pairs at one cost, so it is rechecked after
//! a link. Finally, remaining ends are paired with zero overlap and written
//! as a masked superstring.

use crate::{default_msfa_output, log_file_stats, mss::MssKey, timing::StageTiming};
use packed_seq::{PackedSeqVec, SeqVec};
use rayon::prelude::*;
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
    build_distances: u64,
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

/// Settled distances through the current layer, plus valid tentative labels.
struct DistanceField {
    /// Distance to each oriented head, excluding that unitig's traversal.
    distance: Vec<u8>,
    /// Scratch buffer `remove_source`.
    worklist: Vec<(u32, u8)>,
}

impl DistanceField {
    /// Heads adjacent to a free outgoing tail have distance zero.
    fn new(links: &[Link], graph: &Graph, visited: &mut u64) -> Self {
        let mut field = Self {
            distance: vec![UNREACHABLE; links.len()],
            worklist: Vec::new(),
        };
        *visited += links.len() as u64;
        for source in 0..links.len() as u32 {
            if links[tail_slot(source)].is_empty() {
                let successors = graph.outgoing(source);
                *visited += successors.len() as u64;
                for &next in successors {
                    field.distance[next as usize] = 0;
                }
            }
        }
        field
    }

    /// Relax the nodes currently labeled with this distance.
    fn relax_at(
        &mut self,
        distance: usize,
        graph: &Graph,
        ranges: &[Range<usize>],
        k: usize,
        radius: usize,
        bridges: &mut BridgeQueue,
        first_unprocessed_cost: usize,
        visited: &mut u64,
        meet_visited: &mut u64,
    ) -> usize {
        let mut settled = 0;
        *visited += self.distance.len() as u64;
        for node in 0..self.distance.len() as u32 {
            if self.distance[node as usize] as usize != distance {
                continue;
            }
            settled += 1;
            self.relax_from(node, distance, graph, ranges, k, radius, visited);
            bridges.schedule(
                node & !1,
                self,
                ranges,
                k,
                distance,
                first_unprocessed_cost,
                meet_visited,
            );
        }
        settled
    }

    fn relax_from(
        &mut self,
        node: u32,
        distance: usize,
        graph: &Graph,
        ranges: &[Range<usize>],
        k: usize,
        radius: usize,
        visited: &mut u64,
    ) {
        let weight = traversal_weight(ranges, node, k);
        if weight >= radius {
            return;
        }
        let next_distance = distance + weight;
        if next_distance >= radius {
            return;
        }
        let successors = graph.outgoing(node);
        *visited += 1 + successors.len() as u64;
        for &next in successors {
            if next_distance < self.distance[next as usize] as usize {
                self.distance[next as usize] = next_distance as u8;
            }
        }
    }

    /// Find a replacement through a free tail or a predecessor in a reached
    /// layer. Later tentative predecessors have not been expanded yet.
    fn replacement_distance(
        &self,
        node: u32,
        links: &[Link],
        graph: &Graph,
        ranges: &[Range<usize>],
        k: usize,
        radius: usize,
        reached_layer: usize,
        visited: &mut u64,
    ) -> u8 {
        let predecessors = graph.incoming(node);
        *visited += 1 + predecessors.len() as u64;
        let mut best = radius;
        for predecessor in predecessors {
            if links[tail_slot(predecessor)].is_empty() {
                return 0;
            }
            let prior = self.distance[predecessor as usize] as usize;
            if prior <= reached_layer {
                best = best.min(prior + traversal_weight(ranges, predecessor, k));
            }
        }
        if best == radius {
            UNREACHABLE
        } else {
            best as u8
        }
    }

    /// Repair all finite labels affected by deleting a source, including
    /// tentative labels for later layers.
    fn remove_source(
        &mut self,
        source: u32,
        links: &[Link],
        graph: &Graph,
        ranges: &[Range<usize>],
        k: usize,
        radius: usize,
        reached_layer: usize,
        first_unprocessed_cost: usize,
        bridges: &mut BridgeQueue,
        visited: &mut u64,
        meet_visited: &mut u64,
    ) {
        debug_assert!(!links[tail_slot(source)].is_empty());
        self.worklist.clear();
        // Removing a free tail removes its zero-cost edges to adjacent heads.
        let successors = graph.outgoing(source);
        *visited += 1 + successors.len() as u64;
        for &next in successors {
            if self.distance[next as usize] != 0 {
                continue;
            }
            let replacement = self.replacement_distance(
                next,
                links,
                graph,
                ranges,
                k,
                radius,
                reached_layer,
                visited,
            );
            if replacement > 0 {
                self.distance[next as usize] = replacement;
                bridges.schedule(
                    next & !1,
                    self,
                    ranges,
                    k,
                    reached_layer,
                    first_unprocessed_cost,
                    meet_visited,
                );
                self.worklist.push((next, 0));
            }
        }
        while let Some((node, old_distance)) = self.worklist.pop() {
            let next_nodes = graph.outgoing(node);
            *visited += 1 + next_nodes.len() as u64;
            for &next in next_nodes {
                let old_next = self.distance[next as usize];
                if old_next != UNREACHABLE
                    && old_distance as usize + traversal_weight(ranges, node, k)
                        == old_next as usize
                {
                    let new_next = self.replacement_distance(
                        next,
                        links,
                        graph,
                        ranges,
                        k,
                        radius,
                        reached_layer,
                        visited,
                    );
                    debug_assert!(new_next >= old_next);
                    if new_next > old_next {
                        self.distance[next as usize] = new_next;
                        bridges.schedule(
                            next & !1,
                            self,
                            ranges,
                            k,
                            reached_layer,
                            first_unprocessed_cost,
                            meet_visited,
                        );
                        self.worklist.push((next, old_next));
                    }
                }
            }
        }
    }

    /// Follow tight incoming edges to a zero-distance head, then return one
    /// of its free incoming tails. Each step crosses the predecessor unitig.
    fn source(
        &self,
        mut node: u32,
        links: &[Link],
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
            let mut predecessors = graph.incoming(node);
            *visited += predecessors.len() as u64;
            if distance == 0 {
                return predecessors
                    .find(|&predecessor| links[tail_slot(predecessor)].is_empty())
                    .expect("zero distance has no free incoming tail");
            }
            node = predecessors
                .find(|&predecessor| {
                    let prior = self.distance[predecessor as usize];
                    prior != UNREACHABLE
                        && prior as usize + traversal_weight(ranges, predecessor, k)
                            == distance as usize
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

/// Read the cost of a bridge from the distances to its two oriented heads.
fn bridge_candidate(
    node: u32,
    field: &DistanceField,
    ranges: &[Range<usize>],
    k: usize,
    settled_limit: usize,
    visited: &mut u64,
) -> Option<usize> {
    *visited += 1;
    let forward = field.distance[node as usize];
    let backward = field.distance[(node ^ 1) as usize];
    (forward as usize <= settled_limit && backward as usize <= settled_limit)
        .then(|| forward as usize + traversal_weight(ranges, node, k) + backward as usize)
}

/// One pending bridge per physical unitig. A bridge is queued only after both
/// head labels have reached a layer, so later relaxation cannot lower its cost.
struct BridgeQueue {
    buckets: Vec<Vec<u32>>,
    // queued: Vec<u64>,
}

impl BridgeQueue {
    fn new(k: usize) -> Self {
        Self {
            buckets: vec![Vec::new(); k],
        }
    }

    fn schedule(
        &mut self,
        node: u32,
        field: &DistanceField,
        ranges: &[Range<usize>],
        k: usize,
        reached_layer: usize,
        first_unprocessed_cost: usize,
        visited: &mut u64,
    ) {
        debug_assert_eq!(node & 1, 0);
        let cost = bridge_candidate(node, field, ranges, k, reached_layer, visited);
        self.schedule_known(node, cost, k, first_unprocessed_cost);
    }

    fn schedule_known(
        &mut self,
        node: u32,
        cost: Option<usize>,
        k: usize,
        first_unprocessed_cost: usize,
    ) {
        let Some(cost) = cost else { return };
        if cost < first_unprocessed_cost || cost >= k - 1 {
            return;
        }
        self.buckets[cost].push(node);
    }

    fn pop(&mut self, cost: usize) -> Option<u32> {
        let node = self.buckets[cost].pop()?;
        Some(node)
    }
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
        "Meeting search: grow source and receiver distances below {} one layer at a time; a path through one unitig has cost left distance + unitig extensions + right distance",
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
        "Node visits count repeated inspections: direct_overlap scans graph neighbors at distance zero; build_distances scans and relaxes distance layers; meet checks and queues bridges and traces endpoints; connect touches link slots; remove_source repairs distance labels"
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
            break;
        }
    }
    let mut field = DistanceField::new(&links, graph, &mut visits.build_distances);
    let mut bridges = BridgeQueue::new(k);
    let initial_frontier = field.relax_at(
        0,
        graph,
        ranges,
        k,
        radius,
        &mut bridges,
        1,
        &mut visits.build_distances,
        &mut visits.meet,
    );
    info!(
        "Half-distance 0: {} nodes settled",
        compact(initial_frontier)
    );
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
        "After distance {} node visits: direct_overlap {} (total {}), build_distances {}, meet {} (total {}), connect {} (total {}), remove_source {} (total {})",
        compact(0u64),
        compact(visits.direct_overlap),
        compact(visits.direct_overlap),
        compact(visits.build_distances),
        compact(visits.meet),
        compact(visits.meet),
        compact(visits.connect),
        compact(visits.connect),
        compact(visits.remove_source),
        compact(visits.remove_source),
    );

    // Reached nodes schedule bridges once both sides are known. Process the
    // cost buckets in increasing order, rechecking entries after links.
    for half_distance in 1..=radius {
        if free_ends == 0 {
            break;
        }
        let mut before = visits;
        let frontier = if half_distance < radius {
            field.relax_at(
                half_distance,
                graph,
                ranges,
                k,
                radius,
                &mut bridges,
                2 * half_distance - 1,
                &mut visits.build_distances,
                &mut visits.meet,
            )
        } else {
            0
        };
        info!(
            "Half-distance {}: {} nodes settled",
            compact(half_distance),
            compact(frontier),
        );
        for distance in [2 * half_distance - 1, 2 * half_distance] {
            if distance >= k - 1 || free_ends == 0 {
                break;
            }
            let found = bridges.buckets[distance].len();
            let mut examined = 0usize;
            let mut made = 0usize;
            while let Some(node) = bridges.pop(distance) {
                examined += 1;
                let cost =
                    bridge_candidate(node, &field, ranges, k, half_distance, &mut visits.meet);
                if cost != Some(distance) {
                    bridges.schedule_known(node, cost, k, distance);
                    continue;
                }
                let source = field.source(node, &links, graph, ranges, k, &mut visits.meet);
                let receiver_end =
                    field.source(node ^ 1, &links, graph, ranges, k, &mut visits.meet);
                let target = receiver_end ^ 1;
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
                field.remove_source(
                    source,
                    &links,
                    graph,
                    ranges,
                    k,
                    radius,
                    half_distance,
                    distance,
                    &mut bridges,
                    &mut visits.remove_source,
                    &mut visits.meet,
                );
                if source != target ^ 1 {
                    field.remove_source(
                        target ^ 1,
                        &links,
                        graph,
                        ranges,
                        k,
                        radius,
                        half_distance,
                        distance,
                        &mut bridges,
                        &mut visits.remove_source,
                        &mut visits.meet,
                    );
                }
                // This traversal may connect another pair at the same or a
                // later cost after its sources have been removed.
                bridges.schedule(
                    node,
                    &field,
                    ranges,
                    k,
                    half_distance,
                    distance,
                    &mut visits.meet,
                );
            }
            info!(
                "After distance {}: {} bridge entries queued, {} entries examined, {} links made, {} receiving ends remain, {} free ends, {} estimated bases",
                compact(distance),
                compact(found),
                compact(examined),
                compact(made),
                compact(free_ends),
                compact(free_ends),
                compact(estimated_bases),
            );
            info!(
                "After distance {} node visits: direct_overlap {} (total {}), build_distances {} (total {}), meet {} (total {}), connect {} (total {}), remove_source {} (total {})",
                compact(distance),
                compact(0u64),
                compact(visits.direct_overlap),
                compact(visits.build_distances - before.build_distances),
                compact(visits.build_distances),
                compact(visits.meet - before.meet),
                compact(visits.meet),
                compact(visits.connect - before.connect),
                compact(visits.connect),
                compact(visits.remove_source - before.remove_source),
                compact(visits.remove_source),
            );
            before = visits;
        }
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
        "Final node visits (including closure links): direct_overlap {}, build_distances {}, meet {}, connect {}, remove_source {}",
        compact(visits.direct_overlap),
        compact(visits.build_distances),
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
    let start = output.len();
    if !reverse {
        seq.slice(range.start + overlap..range.end)
            .unpack_into(output);
    } else {
        // FIXME TEST THIS
        seq.slice(range.start..range.end - overlap)
            .unpack_rc_into(output);
    }

    // FIXME TEST THIS
    let num_lowercase = (k - 1 - overlap).min(range.len() - overlap);
    for base in &mut output[start..start + num_lowercase] {
        *base = base.to_ascii_lowercase();
    }
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

    let output = reconstruct_output::<K>(k, seq, ranges, short, links);
    info!(
        "Length summary: {} estimated bases after links, {} output bases",
        compact(estimated_bases),
        compact(output.len()),
    );
    output
}

fn reconstruct_output<K: MssKey>(
    k: usize,
    seq: packed_seq::private::PackedSeqVecBase<2>,
    ranges: Vec<Range<usize>>,
    short: Vec<Range<usize>>,
    links: Vec<Link>,
) -> Vec<u8> {
    info!("Reconstruct output");
    let mut done = vec![false; ranges.len()];
    // make output seqs of length at most 32Mbp, to keep threads busy.
    const PART_BASES: usize = 32 * 1024 * 1024;
    let mut parts = Vec::<Vec<u32>>::new();
    let mut first_overlaps = Vec::<u8>::new();
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
            let mut part = Vec::new();
            let mut part_bases = 0usize;
            let mut first_overlap = 0u8;
            let mut overlap = 0u8;
            let first_part = parts.len();
            loop {
                let index = (id / 2) as usize;
                assert!(!done[index], "link cycle revisits unitig before closing");
                done[index] = true;
                let added_bases = ranges[index].len() - overlap as usize;
                // Keep unitigs intact; an individual contribution may exceed
                // the target part size.
                if !part.is_empty() && part_bases + added_bases > PART_BASES {
                    parts.push(std::mem::take(&mut part));
                    first_overlaps.push(first_overlap);
                    part_bases = 0;
                    first_overlap = overlap;
                }
                part.push(id);
                part_bases += added_bases;
                let link = links[tail_slot(id)];
                assert!(!link.is_empty(), "unmatched unitig end");
                if link.target == id ^ 1 || (link.target / 2) as usize == start {
                    break;
                }
                id = link.target;
                overlap = link.overlap;
            }
            if pass == 1 && parts.len() == first_part {
                // The whole cycle fits in one part. Omit its shortest overlap
                // by opening the cycle at the unitig immediately after it.
                let cut = (0..part.len())
                    .min_by_key(|&i| links[head_slot(part[i])].overlap)
                    .unwrap();
                part.rotate_left(cut);
            }
            parts.push(part);
            first_overlaps.push(first_overlap);
        }
    }
    info!("Reconstructing {} parts in parallel", compact(parts.len()));
    let strings: Vec<Vec<u8>> = parts
        .into_par_iter()
        .zip(first_overlaps)
        .with_max_len(1)
        .map(|(part, first_overlap)| {
            let mut string = Vec::new();
            for (position, id) in part.into_iter().enumerate() {
                let overlap = if position == 0 {
                    first_overlap
                } else {
                    links[head_slot(id)].overlap
                };
                append(
                    &mut string,
                    &seq,
                    &ranges[(id / 2) as usize],
                    id & 1 != 0,
                    overlap as usize,
                    k,
                );
            }
            string
        })
        .collect();
    info!("Concatenating strings..");
    let mut output = Vec::with_capacity(strings.iter().map(Vec::len).sum());
    for string in strings {
        output.extend_from_slice(&string);
    }
    // FIXME: These should be separate contigs, or we should insert padding characters.
    for range in &short {
        if !range.is_empty() {
            append(&mut output, &seq, range, false, 0, k);
        }
    }
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
        .unwrap_or_else(|| default_msfa_output(input, "greedytigs"));
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
