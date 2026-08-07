"""
ContigGeneration
=================
Assembles contigs from a collection of k-mers with imperfect coverage,
where a full Eulerian path across the whole genome isn't possible.

Pipeline:
- Build a de Bruijn graph from the k-mers/reads: each k-mer becomes a directed
  edge from its prefix (k-1)-mer to its suffix (k-1)-mer.
- Break the graph into maximal non-branching paths (plus isolated
  cycles) rather than one single Eulerian path, since gaps in coverage
  mean the graph isn't guaranteed to have one.
- Each maximal non-branching path corresponds to a contig: a stretch of
  genome we can reconstruct with confidence, even if we can't stitch
  every contig together into the full genome (we also don't know the order of contigs in genome and to which strand any belongs).
- Glue each path's (k-1)-mers into a contig by taking the first node in
  full, then appending only the last character of every subsequent node.
"""

def read_patterns(file):
    """
    Reads the input file.
    First line: the k-mers themselves (Patterns), separated by " ", which have imperfect coverage (so we won't search for Eulerian path for whole genome).
    """
    with open(file, "r") as f:
        for line in f: 
            line.rstrip()
    
    patterns = line.split()
    return patterns

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


def maximal_non_branching_paths(adj):
    """
    adj: {node: [list of outgoing neighbors]} 
    returns: list of paths, covering all maximal non-branching paths plus isolated cycles.
    """
    # Collect all nodes (including ones that only appear as targets)
    nodes = set(adj.keys())
    for ws in adj.values():
        nodes.update(ws)

    outdeg = {v: len(adj.get(v, [])) for v in nodes}
    indeg = {v: 0 for v in nodes}
    for v, ws in adj.items():
        for w in ws:
            indeg[w] += 1

    def is_1in1out(v): #nested function, so it can access indeg and outdeg
        return indeg.get(v, 0) == 1 and outdeg.get(v, 0) == 1 #return True or False

    # Mutable copy of adjacency lists so we can "consume" edges
    adj_copy = {v: list(ws) for v, ws in adj.items()}

    paths = []
    visited = set()  # 1-in-1-out nodes already used inside some path

    # 1. Paths starting at non-(1-in-1-out) nodes
    for v in nodes:
        if not is_1in1out(v) and outdeg.get(v, 0) > 0:
            for w in list(adj_copy.get(v, [])): # list() creates one new copy list containing the same elements, so when we remove w, iteration does not get affected
                adj_copy[v].remove(w)  # consume this edge instance
                path = [v, w]
                cur = w
                while is_1in1out(cur):
                    visited.add(cur)
                    nxt = adj_copy[cur].pop()
                    path.append(nxt)
                    cur = nxt
                paths.append(path)

    # 2. Isolated cycles: all-1-in-1-out components not yet visited
    for v in nodes:
        if is_1in1out(v) and v not in visited and adj_copy.get(v):
            cycle = [v]
            visited.add(v)
            cur = v
            while True:
                nxt = adj_copy[cur].pop()
                cycle.append(nxt)
                if nxt == v:
                    break # assumes we will eventually loop back to v, making a cycle
                visited.add(nxt)
                cur = nxt
            paths.append(cycle)

    return paths

def path_to_contig(path):
    """
    `path` is a list of (k-1)-mers, e.g. ["GGC", "GCT", "CTT", "TTA", "TAC", "ACC", "CCA"]
    Consecutive nodes overlap in k-2 characters (since each edge/k-mer glued
    prefix->suffix by shifting one character over).

    To reconstruct the contig (a long contigous segment of genome):
    - Start with the FULL first node
    - For every subsequent node, only its LAST character is genuinely "new"
      information -- everything before that was already included as the
      overlap with the previous node.
    """
    contig = path[0]                 # start with the entire first (k-1)-mer
    for node in path[1:]:
        contig += node[-1]           # append only the last character of each next node
    return contig

if __name__ == "__main__":
    patterns=read_patterns("euler_week_2/sample_input.txt")
    graph=de_bruijn_graph(patterns)
    paths = maximal_non_branching_paths(graph)
    contigs=[]
    for path in paths: 
        contigs.append(path_to_contig(path)) #find contig from every maximal nonbranching path one at a time
    print(*contigs) #unpack the list of contigs
    