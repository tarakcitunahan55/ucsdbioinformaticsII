def read_kmers(file):
    with open(file, "r") as f:
        line = f.readline().rstrip()
        kmers = line.split()
    return kmers

kmers = read_kmers("genomeassembly_week_1/debrujin_from_kmers.txt") #kmers/reads given as inputs in any order
for kmer in kmers:
    k = len(kmer)


def prefix(kmer):
    return kmer[:k-1]

def suffix(kmer):
    return kmer[1:]

def build():
    """Get prefix and suffix of each kmer"""
    holder={}
    for kmer in kmers:
        pre=prefix(kmer)
        suf=suffix(kmer)
        if pre not in holder:
            holder[pre]=[suf]
        else:
            holder[pre].append(suf) #if prefix is already present with a suffix, then append the new suffix (don't try to create a new element in holder and don't change the value of the existing one)

    return holder
 #further you can connect k-1mers (pre and suf) as a Eulerian path to reconstruct genome

debrujin_graph = build()
for prefix, suffix in debrujin_graph.items():
    print(f"{prefix}: {' '.join(suffix)}")



