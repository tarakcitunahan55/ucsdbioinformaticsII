def read_graph(file):
    with open(file, "r") as f:
        graph = {}
        for line in f:
            line = line.rstrip()
            if line == "":
                continue
            node, neighbors = line.split(": ")
            graph[node] = neighbors.split()  # returns a list of neighbors
    return graph


def eulerian_cycle(graph, start):
    # Make a mutable copy of the adjacency lists so we can remove edges without changing the original graph -> "augmented" in this case
    remaining = {node: list(neighbors) for node, neighbors in graph.items()}

    stack = [start] #for Eulerian path, we can't start anywhere — start at the node with the extra outgoing edge (passed into function as argument)
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
    """figures out where the Eulerian path must begin and end"""
    # Collect every node that appears anywhere (as a key OR as someone's neighbor)
    all_nodes = set(graph.keys()) #graph keys already do not have duplicates, but we are using set for removing duplicates when we update
    for neighbors in graph.values():
        all_nodes.update(neighbors) #combined collection will have overlaps, so we need something that automatically discards duplicates when merged. That's exactly what a set gives you via .update() -> unpacks the neighbors list and adds to set (.add() cannot unpack the list, so it would raise error)

    out_degree = {node: len(graph.get(node, [])) for node in all_nodes} #For each node in all_nodes, look up its neighbor list in graph (using .get(node, []) so that if the node isn't a key at all, we get an empty list (len=0) instead of a KeyError, and take its length. len(graph[node]) would crash on any node that's never a key in graph
    in_degree = {node: 0 for node in all_nodes}
    for node, neighbors in graph.items():
        for n in neighbors:
            in_degree[n] += 1 #increment in_degree[n] by 1, because n is being "arrived at."

    #Find the node whose out-degree exceeds in-degree by 1 (start), and vice versa (end)
    start = end = None
    for node in all_nodes:
        diff = out_degree[node] - in_degree[node]
        if diff == 1:
            start = node
        elif diff == -1:
            end = node

    return start, end


def eulerian_path(graph):
    start, end = find_start_and_end(graph)

    if start is None:
        #This branch only matters for graphs that are already perfectly balanced, have an Eulerian cycle rather than a strict path — in that case any node can serve as both start and end.
        start = end = next(iter(graph))

    # Add the artificial edge end -> start so every node becomes balanced (Eulerian cycle)
    augmented = {node: list(neighbors) for node, neighbors in graph.items()} #makes a copy of graph
    augmented.setdefault(end, []).append(start) #augmented.setdefault(end, []) — since end = "4" and "4" isn't currently a key in augmented at all, this creates the entry augmented["4"] = [] and returns that new empty list. .append(start) — appends "6" onto that list. Net effect: augmented["4"] = ["6"]. We've added a brand-new edge 4 → 6 that didn't exist in the original graph.

    # Run the ordinary Eulerian cycle algorithm, but force it to begin at "start"
    cycle = eulerian_cycle(augmented, start) #closed walk

    # Find the artificial edge (end -> start) inside the cycle and cut cycle there to turn it into the Eulerian path.
    for i in range(len(cycle) - 1):
        if cycle[i] == end and cycle[i + 1] == start: #Loops through consecutive pairs in cycle
            # The path is: everything after the artificial edge, followed by everything up to and including the artificial edge's tail.
            path = cycle[i + 1:] + cycle[1:i + 1]
            return path

    # If we never find the artificial edge separately (e.g. graph was already a true cycle and end == start), the cycle itself IS a valid answer.
    return cycle


if __name__ == "__main__":
    graph = read_graph("euler_week_2/euler.txt")
    path = eulerian_path(graph)
    print(" ".join(path))