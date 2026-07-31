def read_data(file):
    with open(file, "r") as f:
        line = f.readline().rstrip() #all inputs are one one line (k and whole genome sequence)
        args=line.split()
        k=int(args[0])
        sequence=args[1]
    return k, sequence

k, sequence = read_data("genomeassembly_week_1/debrujin.txt")

kmers=[sequence[i:i+k] for i in range(len(sequence)-k+1)] #get all possible kmers/reads from the sequence 

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


debrujin_graph = build()
for prefix, suffix in debrujin_graph.items():
    print(f"{prefix}: {' '.join(suffix)}")


