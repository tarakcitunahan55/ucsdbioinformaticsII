genome="TAATGCCATGGGATGTT"
k=3
d=2

def paired_comp():
    holder=[]
    for i in range(len(genome)-(2*k+d)+1):
        holder.append(f"({genome[i:i+k]}|{genome[i+k+d:i+2*k+d]})") #read pair kmers separated by length d
    holder.sort()
    return holder

print (*paired_comp())

