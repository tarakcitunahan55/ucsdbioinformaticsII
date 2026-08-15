"""
CyclopeptideSequencingWithScoring
==================================
Implements the Leaderboard Cyclopeptide Sequencing algorithm for
identifying an unknown cyclic peptide from its experimental mass
spectrum.

A peptide is represented as a list of
integer amino acid masses rather than letters.

The algorithm grows candidate peptides one amino acid at a time
(branch-and-bound):

1. Start with a single empty candidate peptide.
2. Expand: extend every candidate by each of the 18 possible
   amino acid masses.
3. Bound: discard any candidate whose total mass exceeds the
   parent mass (the full peptide's mass, taken from the experimental
   spectrum) — it can no longer be part of a valid answer.
4. Score: for candidates whose mass exactly matches the parent
   mass, compare their theoretical cyclic spectrum against the
   experimental spectrum and keep track of the best-scoring
   leader peptide seen so far.
5. Trim: cut the candidate list down to the top N scorers (by
   linear, not cyclic, spectrum score) to keep the search tractable,
   since the full branching tree grows exponentially.
6. Repeat steps 2-5 until no candidates remain.

The result is the best scoring peptide found, expressed as a list of
integer amino acid masses (since several amino acids share the
same integer mass, the exact letter sequence can't always be
recovered uniquely). Usage of best-scoring instead of looking for direct spectra match
is more realistic, accounting for experimental spectra errors. Here, assume that parent mass is found exactly in
experimental spectrum (which is easy to find in real practice).
"""
from collections import Counter

amino_acid_mass = [
    57, 71, 87, 97, 99, 101, 103, 113, 114,
    115, 128, 129, 131, 137, 147, 156, 163, 186
]


def peptide_mass(peptide):
    """Total mass of a peptide."""
    return sum(peptide)


def linear_spectrum(peptide):
    """List of all contiguous (non-wrapping) linear subpeptide masses, plus 0."""
    n = len(peptide)
    prefix_mass = [0] * (n + 1)
    for i in range(1, n + 1):
        prefix_mass[i] = prefix_mass[i - 1] + peptide[i - 1]

    spectrum = [0]
    for i in range(n):
        for j in range(i + 1, n + 1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])
    return spectrum


def cyclic_spectrum(peptide):
    """List of all subpeptide masses including wrap-around ones, plus 0 and total mass."""
    n = len(peptide)
    prefix_mass = [0] * (n + 1)
    for i in range(1, n + 1):
        prefix_mass[i] = prefix_mass[i - 1] + peptide[i - 1]
    total_mass = prefix_mass[n]

    spectrum = [0]
    for i in range(n):
        for j in range(i + 1, n + 1):
            sub_mass = prefix_mass[j] - prefix_mass[i]
            spectrum.append(sub_mass)
            if i > 0 and j < n:
                spectrum.append(total_mass - sub_mass)
    return spectrum


def expand(peptides):
    """Branching step: extend every candidate peptide by each possible amino acid mass."""
    return [peptide + [a] for peptide in peptides for a in amino_acid_mass]


def read_experimental_spectrum(file):
    """Reads a single line of space separated integer masses."""
    with open(file, "r") as f:
        return [int(x) for x in f.readline().split()]


def score(theoretical_spectrum, experimental_spectrum):
    """Count the number of matching masses. Higher score means higher resemblance."""
    match=0

    for mass1, count1 in Counter(theoretical_spectrum).items():
        for mass2, count2 in Counter(experimental_spectrum).items():
            if mass1==mass2:
                match+=min(count1,count2) #if masses are the same in both spectra, take the lowest count and add to score

    return match


def trim(N, experimental_spectrum, leaderboard):
    """(Linear) Score every peptide on the leaderboard and keep only the top N (with ties)."""
    linear_scores={}

    for peptide in leaderboard:
        theoretical_spectrum = linear_spectrum(peptide)
        linear_scores[tuple(peptide)]=score(theoretical_spectrum,experimental_spectrum) #since peptide is a list, it can't be a key in dict -> make tuple

    sorted_scores = sorted(linear_scores.items(), key=lambda x: x[1], reverse=True) #gives a sorted list of tuples based on highest to lowest score -> [((113,128), 32), ((97,71), 15), ((186,), 8)]

    if N >= len(sorted_scores): #keep everyone if leaderboard is smaller than N
        #trimming is naturally skipped only in the very early rounds — mainly the first one or two iterations, when peptides are still length 0 or 1. Since expand multiplies the leaderboard size by 18 each round, it takes almost no iterations before the board size blows way past N (e.g. even N=236 gets exceeded after just 2 expansions: 1 → 18 → 324). So this guard clause is really only a startup safety net — after that, trim is doing real work every round
        return leaderboard

    (x,y)=sorted_scores[N-1] # N-1 since list indexing starts counting from 0
    trimmed_leaderboard=[list(a) for (a,b) in sorted_scores if b>=y]
    leaderboard=trimmed_leaderboard
    return leaderboard


def leaderboard_cyclopeptide_sequencing(spectrum, N):
    """Returns leader peptide whose cyclic spectrum best matches the experimental spectrum."""
    parent_mass = max(spectrum)

    leaderboard = [[]]
    leaderpeptide = []

    while leaderboard:
        leaderboard = expand(leaderboard)
        # filter out overweight peptides without mutating while iterating (do not use ".remove()", it mutates the list while iterating over the for loop, altering indexes and causing error) (.pop() and del also cause the same error)
        leaderboard = [p for p in leaderboard if peptide_mass(p) <= parent_mass]

        for peptide in leaderboard:
            if peptide_mass(peptide) == parent_mass:
                if score(cyclic_spectrum(peptide), spectrum) > score(cyclic_spectrum(leaderpeptide), spectrum):
                    leaderpeptide = peptide

        leaderboard = trim(N, spectrum, leaderboard)

    return leaderpeptide


if __name__ == "__main__":
    experimental_spectrum = read_experimental_spectrum("antibiotics_week_3/spectrum.txt")
    leaderpeptide = leaderboard_cyclopeptide_sequencing(experimental_spectrum, 161)
    print("-".join(map(str,leaderpeptide)))