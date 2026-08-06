"""
MaximalNonBranchingPaths
=========================
Breaks a directed graph into all maximal non-branching paths.

A node is "1-in-1-out" if it has exactly one incoming and one outgoing
edge. A maximal non-branching path starts/ends at nodes that are not
1-in-1-out, with only 1-in-1-out nodes in between.

Pipeline:
- Compute in-degree and out-degree for every node.
- From each non-1-in-1-out node (with outdeg>0), follow each outgoing edge forward
  through consecutive 1-in-1-out nodes until hitting another branch
  point or dead end (non-1-in-1-out).
- Separately catch isolated cycles: components where every node is
  1-in-1-out, which the step above never starts a path from.
"""

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


def read_adjacency_list_from_file(filepath):
    adj = {}
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            left, right = line.split(":")
            v = int(left.strip())
            ws = [int(x) for x in right.strip().split()]
            adj[v] = ws
    return adj


if __name__ == "__main__":
    adj = read_adjacency_list_from_file("euler_week_2/sample_input.txt")
    result = maximal_non_branching_paths(adj)
    for path in result:
        print(" ".join(map(str, path)))