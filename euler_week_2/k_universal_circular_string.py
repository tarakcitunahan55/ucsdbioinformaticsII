"""Find a circular binary string of length 2^k that contains every possible binary string of length k 
exactly once as a substring"""
from itertools import product
k=8

def all_kmers(k):
    kmers = []
    for p in product('01', repeat=k): #generate all possible kmers
        kmers.append(''.join(p))
    return kmers

def de_bruijn_graph(patterns):
    graph = {}
    for pattern in patterns:
        prefix = pattern[:-1]   # first k-1 characters
        suffix = pattern[1:]    # last k-1 characters
        graph.setdefault(prefix, []).append(suffix) #if prefix isn't already a key in graph, this creates it with an empty list as its value and returns that list; if it is already a key, it just returns the existing list. Either way, .append(suffix) then adds suffix onto that list
    return graph

def eulerian_cycle(graph):
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

def cycle_to_genome(cycle):
    """
    - Start with the FULL first node 
    - For every subsequent node, only its LAST character is genuinely "new"
      information -- everything before that was already included as the
      overlap with the previous node.
    """
    genome = cycle[0]                 # start with the entire first (k-1)-mer
    for node in cycle[1:]:
        genome += node[-1]           # append only the last character of each next node

    return genome[:-k+1] # -k+1 = -(k-1)) cut off the last k-1 characters from the end of the string
    #because the cycle repeats its start node (k-1 length) at the very end, this linear spelling redundantly includes it twice: once at the beginning (the starting node itself) and once again at the very end

kmers=all_kmers(k)
graph=de_bruijn_graph(kmers)
cycle = eulerian_cycle(graph)
genome=cycle_to_genome(cycle)
print(genome)

