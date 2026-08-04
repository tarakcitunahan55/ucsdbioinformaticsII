"""
StringSpelledByGappedPatterns
==============================
Reconstructs a string from a sequence of (k,d)-mers -- pairs of k-mers that
are known to be separated by a gap of d unknown characters in the original
genome (this happens with paired-end reads).

Idea (matches the given pseudocode exactly):
- FirstPatterns  = just the "left" k-mer of every (k,d)-mer, in order
- SecondPatterns = just the "right" k-mer of every (k,d)-mer, in order
- Since consecutive FirstPatterns overlap by k-1 characters (same as in
  ordinary genome-path reconstruction), we can glue them all together into
  one PrefixString using the same logic as path_to_genome. Same for
  SecondPatterns -> SuffixString.
- PrefixString and SuffixString describe the SAME underlying genome, just
  shifted relative to one another by (k + d) positions (because that's the
  distance from the start of the left k-mer to the start of the right
  k-mer in each pair).
- So we can check consistency: characters of PrefixString from position
  (k+d) onward MUST match the corresponding earlier-positioned characters
  of SuffixString. If they don't line up, no valid string exists.
- If everything matches, the final genome is: all of PrefixString, plus
  whatever "new" tail information SuffixString has that PrefixString
  didn't reach -- that's exactly the LAST (k+d) characters of SuffixString.
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

    k, d = map(int, lines[0].split())

    gapped_patterns = []
    for token in lines[1].split():
        a, b = token.split("|")       # split "GACC|GCGC" into "GACC", "GCGC"
        gapped_patterns.append((a, b))

    return k, d, gapped_patterns


def path_to_genome(patterns):
    """
    StringSpelledByPatterns (reused from the earlier genome-path code).
    Glues a list of overlapping k-mers into one string: start with the
    first pattern in full, then append only the last character of each
    subsequent pattern (since the rest is guaranteed overlap, because
    consecutive k-mers in the list share a (k-1)-character overlap).
    """
    genome = patterns[0]          # start with the FULL first k-mer
    for p in patterns[1:]:
        genome += p[-1]           # append only the NEW last character each time
    return genome


def string_spelled_by_gapped_patterns(gapped_patterns, k, d):
    """
    Implements the pseudocode directly:

        FirstPatterns  <- initial k-mers from GappedPatterns
        SecondPatterns <- terminal k-mers from GappedPatterns
        PrefixString   <- StringSpelledByPatterns(FirstPatterns, k)
        SuffixString   <- StringSpelledByPatterns(SecondPatterns, k)
        check overlap consistency
        return PrefixString + last (k+d) chars of SuffixString
    """

    # ---- Step 1: split the list of (a, b) pairs into two separate lists ----
    # gapped_patterns looks like: [("GACC","GCGC"), ("ACCG","CGCC"), ...]
    # We only want the "a" halves in order, and separately the "b" halves in order.
    first_patterns = [a for a, b in gapped_patterns]   # e.g. ["GACC","ACCG","CCGA","CGAG","GAGC"]
    second_patterns = [b for a, b in gapped_patterns]  # e.g. ["GCGC","CGCC","GCCG","CCGG","CGGA"]

    # ---- Step 2: glue each list into one long string ----
    # This works because consecutive FirstPatterns entries overlap by k-1
    # characters with each other (same guarantee as ordinary genome-path
    # reconstruction) -- and likewise for SecondPatterns.
    prefix_string = path_to_genome(first_patterns)    # e.g. "GACCGAGC"
    suffix_string = path_to_genome(second_patterns)   # e.g. "GCGCCGGA"

    # ---- Step 3: consistency check ----
    # PrefixString and SuffixString both describe the SAME genome, but
    # SuffixString is "shifted forward" by (k+d) positions relative to
    # PrefixString (because that's exactly the gap between where a left
    # k-mer starts and where its paired right k-mer starts).
    #
    # So: PrefixString's character at position i should equal
    #     SuffixString's character at position (i - k - d),
    # for every i where both strings actually have a character there.
    #
    # (This loop is a direct 0-indexed translation of the pseudocode's
    # 1-indexed "for i = k+d+1 to |PrefixString|" line.)
    for i in range(k + d, len(prefix_string)):
        if prefix_string[i] != suffix_string[i - k - d]:
            return None  # "there is no string spelled by the gapped patterns"

    # ---- Step 4: build the final genome ----
    # PrefixString already correctly covers the whole genome EXCEPT for its
    # very last (k+d) characters, which only SuffixString has "seen" (since
    # SuffixString extends (k+d) positions further to the right).
    # So: take all of PrefixString, then tack on just that missing tail
    # from the end of SuffixString.
    return prefix_string + suffix_string[-(k + d):]


if __name__ == "__main__":
    k, d, gapped_patterns = read_gapped_patterns("euler_week_2/gapped.txt")
    result = string_spelled_by_gapped_patterns(gapped_patterns, k, d)
    print(result if result is not None else "there is no string spelled by the gapped patterns")