def read_kmers(file):
    with open(file, "r") as f:
        line = f.readline().rstrip()
        kmers = line.split()
    return kmers


kmers = read_kmers("genomeassembly_week_1/graph.txt") #kmers/reads given as inputs in any order
for kmer in kmers:
    k = len(kmer)

def prefix():
    kmer_and_prefixes={}
    for kmer in kmers:
        kmer_and_prefixes[kmer] = kmer[:k-1]
    return kmer_and_prefixes

def suffix():
    kmer_and_suffixes={}
    for kmer in kmers:
        kmer_and_suffixes[kmer] = kmer[1:]
    return kmer_and_suffixes

def compare():
    """Find a kmer whose suffix matches another kmer's prefix"""
    result = {}
    kmer_and_prefixes = prefix()
    kmer_and_suffixes = suffix()
    for kmer1, s in kmer_and_suffixes.items():
        for kmer2, p in kmer_and_prefixes.items():
            if s == p and kmer1 != kmer2: 
                if kmer1 in result:
                    result[kmer1].append(kmer2)
                else:
                    result[kmer1] = [kmer2]
    return result

overlap_graph = compare()
for kmer, neighbors in overlap_graph.items():
    print(f"{kmer}: {' '.join(neighbors)}") #further, you can connect the overlapping kmers in a Hamiltonian path (every vertex/node/kmer once) to reconstruct genome (there can be multiple Hamiltonian paths due to repeated kmers, meaning different ways to reconstruct genome)
#The key challenge is that while we are guided by Euler’s Theorem in solving the Eulerian Cycle Problem, an analogous simple condition for the Hamiltonian Cycle Problem remains unknown.
#Use debruijn graph and Euler's theorem instead to find Eulerian path(s)

""" Alternative: CORRECT, BUT SLOWER 
kmers= "example"
def overlap_graph():
    kmer_list = kmers.split() #split returns a list of strings - kmers
    holder={}
    for kmer1 in kmer_list:
        holder[kmer1]=""
        for kmer2 in kmer_list:
            k=len(kmer2)
            if kmer1!=kmer2 and kmer1[1:]==kmer2[:k-1]:
                holder[kmer1]=holder[kmer1]+f" {kmer2}"
    return holder

holder=overlap_graph()

for key,value in holder.items():
    if value: #only print those with a value
        print(f"{key}:{value}")"""



