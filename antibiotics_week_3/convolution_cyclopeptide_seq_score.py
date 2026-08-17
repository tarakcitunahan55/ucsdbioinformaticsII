"""
ConvolutionCyclopeptideSequencing
==================================

Extends Leaderboard Cyclopeptide Sequencing by first inferring the
amino acid alphabet directly from the experimental spectrum, rather
than assuming the 18 standard integer masses (or >100 masses including all proteinogenic and non-proteinogenic aa.s).

Since non-ribosomal peptides (cyclopeptide antibiotics like Tyrocidines in our case) can contain
non-proteinogenic amino acids, the 18-mass table used before isn't
the right alphabet. Since extended table has 144 aa.s, spectral convolution
recovers plausible amino acid masses empirically rather than using all 144 and generating lots of incorrect peptides:

1. Compute the convolution: every positive pairwise difference
   between masses in the spectrum. Differences corresponding to
   real amino acid masses should occur often, since each amino
   acid's mass shows up repeatedly as the gap between many pairs
   of subpeptide masses.
2. Keep only differences in the range 57-200, the plausible mass
   range for a single amino acid (proteinogenic and non-proteinogenic).
3. Keep the M most frequent differences (with ties), producing a
   custom candidate alphabet.

This alphabet is then fed into the same branch-and-bound search:
expand candidate peptides by each mass in the alphabet (now hugely restricted by convolution),
discard any that exceed the parent mass, score full-length
candidates by their cyclic spectrum, and trim intermediate
candidates each round by their linear spectrum score, keeping the
top N (with ties) to keep the search tractable.
As before, the result may include multiple linear peptides that
represent the same underlying cyclic peptide (rotations of one
another)

Caveat: Although we use a very restricted set of amino acids from the given exp. spectrum and convolution makes the code more robust, 
the search is still heuristic and a noisy spectrum can still cause the true peptide to be trimmed early or lose a final tie.
"""
from collections import Counter

def spectral_convolution(spectrum, M):
    n=len(spectrum)
    convolution =[]
    for x in range(0,n-1):
        for y in range(x+1,n):
            diff=spectrum[y]-spectrum[x] #since the spectrum is given in increasing order all diff are >=0
            if 57<=diff<=200: #amino acids can be any integer mass 57-200 since NRPs can have non-proteinogenic aa.s as well as proteinogenic aa.s
                convolution.append(diff)
    aa_freq=Counter(convolution) #frequency dictionary
    sorted_freqs = sorted(aa_freq.items(), key=lambda x: x[1], reverse=True) #gives a sorted list of tuples based on highest to lowest freq -> [(113, 32), (97, 15), (186, 8)]
    
    if M >= len(sorted_freqs): #keep everyone if smaller than M
        return convolution
    
    (x,y)=sorted_freqs[M-1] # M-1 since list indexing starts counting from 0
    m_freq_convolution=[a for (a,b) in sorted_freqs if b>=y] 

    return m_freq_convolution

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


def expand(peptides,convolution):
    """Branching step: extend every candidate peptide by each possible amino acid mass."""
    return [peptide + [a] for peptide in peptides for a in convolution]


def read_experimental_spectrum(file):
    """Reads M (integer for most frequent M aa masses in convolution with ties), N (integer for highest scoring N peptides with ties) 
    and a single line of space separated integer masses in a total of three lines."""
    with open(file, "r") as f:
        M=int(f.readline().strip())
        N=int(f.readline().strip())
        spectrum=[int(x) for x in f.readline().split()] #list(map(int,f.readline().split())) 
        return M, N, spectrum


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


def leaderboard_cyclopeptide_sequencing(spectrum, N, convolution):
    """Returns all leader peptides whose cyclic spectrum best matches the experimental spectrum."""
    parent_mass = max(spectrum)

    leaderboard = [[]]
    leaderpeptides = []
    leader_score = 0

    while leaderboard:
        leaderboard = expand(leaderboard,convolution)
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
    M, N, experimental_spectrum = read_experimental_spectrum("antibiotics_week_3/spectrum.txt")
    convolution=spectral_convolution(experimental_spectrum, M)
    leaderpeptides = leaderboard_cyclopeptide_sequencing(experimental_spectrum, N, convolution)
    for peptide in leaderpeptides: #each one represents a cyclic peptide, written linearly in answers, so we can get multiple linear answers which represent the same cyclic peptide like -> [1,2,3] and [2,3,1] and [3,1,2]
        print("-".join(map(str,peptide)),end=" ")
