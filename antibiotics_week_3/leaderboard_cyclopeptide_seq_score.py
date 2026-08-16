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

The result is the best scoring peptides found, expressed as lists of
integer amino acid masses (since several amino acids share the
same integer mass, the exact letter sequence can't always be
recovered uniquely). Usage of best-scoring instead of looking for exact spectra match
is more realistic, accounting for experimental spectra errors. Here, assume that parent mass is found exactly in
experimental spectrum (which is easy to find in real practice).

Caveat: as the experimental spectrum gets noisy (more false/missing masses), this heuristic
algorithm can return incorrect peptides.

Trimming is irreversible: once a peptide is cut from the leaderboard, it can never come back.
False masses (noise peaks) can make an incorrect peptide score higher than it should, while
missing masses (true peaks lost to noise) can make the true peptide score lower than it should.
If this pushes the true peptide out of the top N at any intermediate round, it is discarded
for good — even though it may have gone on to be the correct final answer.

Even if the true peptide survives every trimming round, noise can still cause an incorrect
peptide to tie or outscore it at the final comparison, causing the wrong peptide to be
returned as the leader.
"""
from collections import Counter

#18 unique amino acid residue integer masses (unique residue masses of 20 proteinogenic amino acids, excluding pyrrolysine and selenocysteine)
#I, L (113 Da); K, Q (128 Da) have same integer masses
amino_acid_mass = [
    57, 71, 87, 97, 99, 101, 103, 113, 114,
    115, 128, 129, 131, 137, 147, 156, 163, 186
]

"""
# Extended amino acid masses of all proteinogenic and non-proteinogenic amino acids since non-ribosomal peptides/NRPs (our antibiotic cyclopeptide in this case)
# can include non-proteinogenic amino acids as NRPs are not synthesized in ribosomes.
# However, using this algorithm with extended mass table will generate lots of incorrect candidate peptides.
# Instead, we must determine the amino acid composition of a peptide from its spectrum so that we may run LeaderboardCyclopeptideSequencing on this smaller alphabet of amino acids.
amino_acid_mass = list(range(57, 201))"""

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
    trimmed_leaderboard=[list(a) for (a,b) in sorted_scores if b>=y] #"list(a)" is here because of the earlier step we had to convert peptides into tuples to use them as dictionary keys — and now we convert them back into lists, since that's the data type the rest of algorithm expects.
    leaderboard=trimmed_leaderboard
    return leaderboard


def leaderboard_cyclopeptide_sequencing(spectrum, N):
    """Returns all leader peptides whose cyclic spectrum best matches the experimental spectrum."""
    parent_mass = max(spectrum)

    leaderboard = [[]]
    leaderpeptides = []
    leader_score = 0

    while leaderboard:
        leaderboard = expand(leaderboard)
        # filter out overweight peptides without mutating while iterating (do not use ".remove()", it mutates the list while iterating over the for loop, altering indexes and causing error) (.pop() and del also cause the same error)
        leaderboard = [p for p in leaderboard if peptide_mass(p) <= parent_mass]

        for peptide in leaderboard:
            if peptide_mass(peptide) == parent_mass:
                peptide_score = score(cyclic_spectrum(peptide), spectrum)
                if peptide_score > leader_score:
                    leaderpeptides = [peptide]
                    leader_score = peptide_score
                elif peptide_score == leader_score:
                    leaderpeptides.append(peptide)

        leaderboard = trim(N, spectrum, leaderboard)

    return leaderpeptides #if there are multiple highest scoring peptides, all of them are getting returned


if __name__ == "__main__":
    experimental_spectrum = read_experimental_spectrum("antibiotics_week_3/spectrum.txt")
    leaderpeptides = leaderboard_cyclopeptide_sequencing(experimental_spectrum, 1000)
    for peptide in leaderpeptides: #each one represents a cyclic peptide, written linearly in answers, so we can get multiple linear answers which represent the same cyclic peptide like -> [1,2,3] and [2,3,1] and [3,1,2]
        print("-".join(map(str,peptide)),end=" ")
