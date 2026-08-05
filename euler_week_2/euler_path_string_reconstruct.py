"""
String Reconstruction Problem
==============================
Given a collection of k-mers (Patterns) (not in the order they found in genome), reconstruct a string Text (linear genome) 
whose k-mer composition equals Patterns.

Pipeline:
    1. DeBruijn(Patterns)   -> build the de Bruijn graph from the k-mers
    2. EulerianPath(dB)     -> find an Eulerian path through that graph
    3. PathToGenome(path)   -> spell out the linear genome string from that path

"""

def read_patterns(file):
    """
    Reads the input file.
    First line: integer k (the k-mer length).
    Second line: the k-mers themselves (Patterns), separated by " ".
    """
    with open(file, "r") as f:
        lines = [line.rstrip() for line in f if line.rstrip() != ""]
    k = int(lines[0])
    patterns = lines[1].split()
    return k, patterns


def de_bruijn_graph(patterns):
    """
    Builds the de Bruijn graph from a list of k-mers.

    - Each k-mer becomes a single DIRECTED EDGE, not a node.
    - The edge goes FROM the k-mer's PREFIX (first k-1 chars)
                    TO   the k-mer's SUFFIX (last k-1 chars).
    - Nodes are therefore (k-1)-mers, and edges are the original k-mers.

    If the same (prefix, suffix) pair occurs from multiple k-mers
    (duplicate k-mers), we must add a SEPARATE edge each time --
    this is exactly why we use a list (not a set) for neighbors.
    """
    graph = {}
    for pattern in patterns:
        prefix = pattern[:-1]   # first k-1 characters
        suffix = pattern[1:]    # last k-1 characters
        graph.setdefault(prefix, []).append(suffix) #if prefix isn't already a key in graph, this creates it with an empty list as its value and returns that list; if it is already a key, it just returns the existing list. Either way, .append(suffix) then adds suffix onto that list
    return graph


def eulerian_cycle(graph, start):
    """
    Hierholzer's algorithm:
    Consumes edges out of a MUTABLE COPY (`remaining`) as it walks, backtracking (popping onto `cycle`) whenever it dead-ends.
    """
    remaining = {node: list(neighbors) for node, neighbors in graph.items()} # Make a mutable copy of the adjacency lists so we can remove edges without changing the original graph -> "augmented" in this case

    stack = [start]  #for Eulerian path, we can't start anywhere — start at the node with the extra outgoing edge (passed into function as argument)
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


def eulerian_path(graph):
    """
    Finds an Eulerian PATH by temporarily converting the graph into an
    Eulerian CYCLE (adding one artificial edge end -> start), running the
    ordinary cycle algorithm, then "cutting" the cycle open at that
    artificial edge to unroll it into the Eulerian path.
    """
    start, end = find_start_and_end(graph)

    if start is None:
        #This branch only matters for graphs that are already perfectly balanced, have an Eulerian cycle rather than a strict path — in that case any node can serve as both start and end.
        start = end = next(iter(graph))

    # Build augmented graph with the artificial edge end -> start added
    augmented = {node: list(neighbors) for node, neighbors in graph.items()}
    augmented.setdefault(end, []).append(start)

    cycle = eulerian_cycle(augmented, start)

    # Cut the cycle at the artificial edge to recover the real path
    for i in range(len(cycle) - 1):
        if cycle[i] == end and cycle[i + 1] == start:
            path = cycle[i + 1:] + cycle[1:i + 1]
            return path

    return cycle  # if graph was already a true cycle


def path_to_genome(path):
    """
    String Spelled by a Genome Path Problem.

    `path` is a list of (k-1)-mers, e.g. ["GGC", "GCT", "CTT", "TTA", "TAC", "ACC", "CCA"]
    Consecutive nodes overlap in k-2 characters (since each edge/k-mer glued
    prefix->suffix by shifting one character over).

    To reconstruct the genome:
    - Start with the FULL first node (all of path[0]).
    - For every subsequent node, only its LAST character is genuinely "new"
      information -- everything before that was already included as the
      overlap with the previous node.
    """
    genome = path[0]                 # start with the entire first (k-1)-mer
    for node in path[1:]:
        genome += node[-1]           # append only the last character of each next node
    return genome


def string_reconstruction(k, patterns):
    """
    Ties the whole pipeline together
    """
    dB = de_bruijn_graph(patterns)
    path = eulerian_path(dB)
    text = path_to_genome(path)
    return text


if __name__ == "__main__":
    k, patterns = read_patterns("euler_week_2/reconstruction.txt")
    text = string_reconstruction(k, patterns)
    print(text)