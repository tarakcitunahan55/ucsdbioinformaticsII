"""
StringSpelledByGappedPatterns
==============================
Reconstructs a string from a sequence of (k,d)-mers - pairs of k-mers that
are known to be separated by a gap of d unknown characters in the original
genome (paired-end reads) in (location) order.

Pipeline:
- FirstPatterns  = just "left" k-mer of every (k,d)-mer, in order
- SecondPatterns = just "right" k-mer of every (k,d)-mer, in order
- Since consecutive FirstPatterns overlap by k-1 characters (same as in
  ordinary genome-path reconstruction), we can glue them all together into
  one PrefixString using the same logic as path_to_genome. Same for
  SecondPatterns -> SuffixString.
- PrefixString and SuffixString describe the SAME underlying genome, just
  shifted relative to one another by (k + d) positions (distance from the start of the left k-mer to the start of the right
  k-mer in a read pair).
- So we can check consistency: characters of PrefixString from position
  (k+d) onward MUST match the corresponding earlier-positioned characters
  of SuffixString. If they don't line up, no valid string exists.
- If everything matches, the final genome is all of PrefixString, plus
  whatever new tail information SuffixString has that PrefixString
  didn't reach - LAST (k+d) characters of SuffixString.
"""


def read_gapped_patterns(file):
    """
    Reads input:
        line 1: "k d"
        line 2: space-separated (k,d)-mers in (genomic) order formatted as "aaaa|bbbb"
    Returns k, d, and a list of (prefix_kmer, suffix_kmer) tuples.
    """
    with open(file, "r") as f:
        lines = [line.rstrip() for line in f if line.rstrip() != ""]

    k, d = map(int, lines[0].split()) #map() applies int() to both of the elements of line[0] -> k, d

    gapped_patterns = []
    for token in lines[1].split():
        a, b = token.split("|")       # split the first and second kmer of a read pair
        gapped_patterns.append((a, b)) # append tuples

    return k, d, gapped_patterns


def path_to_genome(patterns):
    """
    Glues a list of overlapping k-mers into one string: start with the
    first pattern in full, then append only the last character of each
    subsequent pattern 
    """
    genome = patterns[0]          # start with whole first k-mer
    for p in patterns[1:]:
        genome += p[-1]           # append only the last character each time
    return genome


def string_spelled_by_gapped_patterns(gapped_patterns, k, d):
    """
        Returns PrefixString + last (k+d) chars of SuffixString
    """

    # split the list of (a, b) pairs into two separate lists 
    # gapped_patterns looks like: [("GACC","GCGC"), ("ACCG","CGCC"), ...]
    first_patterns = [a for a, b in gapped_patterns]  
    second_patterns = [b for a, b in gapped_patterns]  

    # glue each list into one long string 
    prefix_string = path_to_genome(first_patterns)    
    suffix_string = path_to_genome(second_patterns)   

    # consistency check (same genome, just a k+d shift) 
    # PrefixString's character at position i = SuffixString's character at position (i - k - d),
    for i in range(k + d, len(prefix_string)):
        if prefix_string[i] != suffix_string[i - k - d]:
            return None  #there is no string spelled by the gapped patterns

    # build the final genome
    return prefix_string + suffix_string[-(k + d):] #last k+d char of suffix string


if __name__ == "__main__":
    k, d, gapped_patterns = read_gapped_patterns("euler_week_2/gapped.txt")
    result = string_spelled_by_gapped_patterns(gapped_patterns, k, d)
    print(result if result is not None else "there is no string spelled by the gapped patterns")