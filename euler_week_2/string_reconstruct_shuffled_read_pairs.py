"""
StringReconstructionFromReadPairs
=================================
Reconstructs a genome from an unordered collection of paired k-mers
(read pairs) separated by a known gap d, assuming perfect coverage.

Pipeline:
- Each read pair (Left|Right) becomes an edge in a paired De Bruijn graph.
- The source node of an edge is (prefix(Left), prefix(Right)), and the
  destination node is (suffix(Left), suffix(Right)).
- Because the read pairs are given in arbitrary order, we first recover
  their genomic order by finding an Eulerian path through the graph.
- To find the Eulerian path, temporarily add one artificial edge from the
  end node to the start node, turning the graph into an Eulerian cycle.
- Find the Eulerian cycle using Hierholzer's algorithm, then remove the
  artificial edge to recover the original Eulerian path.
- Because the graph may contain repeated (k-1)-mers, there can be more
  than one topologically valid Eulerian path, but only one corresponds to
  the true genome (the one consistent with the gap d). We therefore try
  multiple randomized Eulerian paths until we find one that is consistent.
- The ordered path of paired (k-1)-mers is then treated exactly like the
  ordered gapped patterns from the string_reconstruct_from_ordered_read_pairs.py problem.
- Glue together the left (k-1)-mers to obtain PrefixString and the right
  (k-1)-mers to obtain SuffixString.
- PrefixString and SuffixString describe the same genome, shifted by
  (k + d) positions, so their overlapping characters must match.
- If they are consistent, append the last (k+d) characters of
  SuffixString to PrefixString to obtain the reconstructed genome.
"""

import random
from collections import defaultdict

def read_gapped_patterns(file):
    """
    Reads input:
        line 1: "k d"
        line 2: space-separated (k,d)-mers in arbitrary order (unordered) formatted as "aaaa|bbbb"
    Returns k, d, and a list of (prefix_kmer, suffix_kmer) tuples.
    """
    with open(file, "r") as f:
        lines = [line.rstrip() for line in f if line.rstrip() != ""]

    k, d = map(int, lines[0].split()) #map() applies int() to both of the elements of line[0] -> k, d

    gapped_patterns = []
    for token in lines[1].split():
        a, b = token.split("|")       # split the first and second kmer of a read pair
        gapped_patterns.append((a, b)) # append tuples

    return k, d, gapped_patterns


def build_paired_debruijn_graph(gapped_patterns):

    """
    Builds the paired De Bruijn graph.
    """
    graph = defaultdict(list)

    for left, right in gapped_patterns:
        prefix = (left[:-1], right[:-1])
        suffix = (left[1:], right[1:])
        graph[prefix].append(suffix)

    return graph


def find_start_and_end(graph):
    """
    Determines the start and end nodes of the Eulerian PATH (not cycle)
    by finding degree imbalances:
        out_degree - in_degree == +1  -> this node is the path's START
        out_degree - in_degree == -1  -> this node is the path's END
    All other nodes must be balanced (in_degree == out_degree).
    """
    # Gather every node mentioned anywhere: as a key, or as someone's neighbor
    all_nodes = set(graph.keys())
    for neighbors in graph.values():
        all_nodes.update(neighbors)

    out_degree = {node: len(graph.get(node, [])) for node in all_nodes}
    in_degree = {node: 0 for node in all_nodes}
    for node, neighbors in graph.items():
        for n in neighbors:
            in_degree[n] += 1

    start = end = None
    for node in all_nodes:
        diff = out_degree[node] - in_degree[node]
        if diff == 1:
            start = node
        elif diff == -1:
            end = node

    return start, end


def eulerian_cycle(graph, start):
    """
    Hierholzer's algorithm.
    """

    remaining = {node: list(neighbors) for node, neighbors in graph.items()}

    stack = [start]
    cycle = []

    while stack:
        current = stack[-1]

        if current in remaining and remaining[current]:
            next_node = remaining[current].pop()
            stack.append(next_node)
        else:
            cycle.append(stack.pop())

    cycle.reverse()
    return cycle


def eulerian_path(graph):
    """
    Turns the graph into an Eulerian cycle by adding one artificial edge,
    finds the cycle, then cuts the cycle open to recover the path.
    """

    start, end = find_start_and_end(graph)

    # graph already has an Eulerian cycle
    if start is None:
        start = end = next(iter(graph))

    # add artificial edge end -> start
    augmented = {node: list(neighbors) for node, neighbors in graph.items()}
    augmented.setdefault(end, []).append(start)

    cycle = eulerian_cycle(augmented, start)

    # Cut the cycle at the artificial edge to recover the real path
    for i in range(len(cycle) - 1):
        if cycle[i] == end and cycle[i + 1] == start:
            path = cycle[i + 1:] + cycle[1:i + 1]
            return path

    return cycle  # if graph was already a true cycle


def eulerian_path_consistent(graph, k, d, max_tries=2000):
    """
    The paired De Bruijn graph may admit multiple topologically valid
    Eulerian paths (when it contains repeated (k-1)-mer nodes). Only one such path is consistent with the
    gap d and actually spells out the true genome. This function retries
    with randomly shuffled adjacency-list orderings until it finds an
    Eulerian path whose PrefixString/SuffixString are consistent.
    """
    for _ in range(max_tries):
        shuffled = {node: random.sample(neighbors, len(neighbors)) #Hierholzer's algorithm just returns some Eulerian path — determined by whatever order edges happen to sit in the adjacency lists (since remaining[current].pop() just pops from the end of a Python list, which is essentially an arbitrary order
                    for node, neighbors in graph.items()} #randomly shuffle the order of neighbors to generate a different Eulerian path, if multiple exists
        path = eulerian_path(shuffled)
        genome = string_spelled_by_gapped_patterns(path, k, d)
        if genome:
            return genome
    return None #if after max_tries no Eulerian path produced a consistent prefix string and suffix string


def path_to_genome(patterns):
    """
    Glues overlapping k-1mer patterns into one string.
    """
    genome = patterns[0]

    for pattern in patterns[1:]:
        genome += pattern[-1]

    return genome


def string_spelled_by_gapped_patterns(path, k, d):
    """
        Returns PrefixString + last (k+d) chars of SuffixString (reconstructing genome)
    """


    first_patterns = [a for a, b in path]  
    second_patterns = [b for a, b in path]  

    # glue each list into one long string 
    prefix_string = path_to_genome(first_patterns)    
    suffix_string = path_to_genome(second_patterns)   

    # consistency check (same genome, just a k+d shift) 
    # PrefixString's character at position i = SuffixString's character at position (i - k - d),
    for i in range(k + d, len(prefix_string)):
        if prefix_string[i] != suffix_string[i - k - d]:
            return None  #there is no string spelled by the gapped patterns

    # build the final genome
    return prefix_string + suffix_string[-(k + d):] #last k+d char of suffix string



if __name__ == "__main__":

    k, d, gapped_patterns = read_gapped_patterns("euler_week_2/gapped.txt")

    graph = build_paired_debruijn_graph(gapped_patterns)

    genome = eulerian_path_consistent(graph, k, d)

    print(genome if genome else "No genome can be reconstructed.")