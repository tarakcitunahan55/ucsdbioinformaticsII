def read_graph(file):
    with open(file, "r") as f:
        graph = {}
        for line in f:
            line = line.rstrip()
            if line == "":
                continue
            node, neighbors = line.split(": ")
            graph[node] = neighbors.split() #split() returns a list of neighbors -> it maps node: [list of neighbor nodes]. Each entry means "there's a directed edge from node to each neighbor in the list
    return graph


def eulerian_cycle(graph):
    """reads a directed graph from a file (given as adjacency lists) and finds an Eulerian cycle — a closed walk that uses every edge exactly once 
    — using Hierholzer's algorithm.
    Assumes the graph actually has an Eulerian cycle. It doesn't check the standard prerequisites (every node's in-degree equals its out-degree (balanced), and the graph is strongly connected)"""

    # Make a mutable copy of the adjacency lists so we can remove edges with .pop as we use them without changing original graph
    remaining = {node: list(neighbors) for node, neighbors in graph.items()} #list(neighbors) copies every list, so "remaining" is a new dict with new lists

    start = next(iter(graph)) #"iter(graph)" creates an iterator over the graph's keys and "next" returns the first key/node
    stack = [start] #current path being explored (a list acting as a stack)
    cycle = [] #accumulates the final Eulerian cycle, built in reverse order

    while stack: #continues until the stack is empty — meaning every edge has been consumed and fully backtracked through
        current = stack[-1] #Look at the node on top of the stack (current)
        if current in remaining and remaining[current]: #If current still has unused outgoing edges. Note: "current in remaining" guards against nodes that appear only as destinations and never as keys in graph, i.e., nodes with no outgoing edges
            next_node = remaining[current].pop() #.pop() removes and returns the last item in a list (Pop one neighbor off that list (this "consumes" that edge)) Note: pop() → removes and returns the last element. pop(i) → removes and returns the element at index i.
            stack.append(next_node) #Push that neighbor onto the stack — we now walk into it
        else: #no more unused edges from current
            cycle.append(stack.pop()) #Pop current off the stack and append it to cycle. Classic Hierholzer move: when you get stuck, that node is finalized as part of the cycle, and you backtrack to the previous node on the stack to look for other unused edges from it.

    cycle.reverse() #Since nodes were appended to cycle in the order they were "finished" (dead-ended), which is the reverse of the actual traversal order, we reverse the list to get the correct Eulerian cycle from start to finish.
    return cycle


graph = read_graph("euler_week_2/euler.txt")
cycle = eulerian_cycle(graph)
print(" ".join(cycle))
