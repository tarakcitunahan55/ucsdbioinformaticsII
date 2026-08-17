"""
LinearPeptideScoring
========================
Linear_score will score the growing linear peptides from Leaderboard (here we have just one as an example) 
against the experimental spectrum.
"""
from collections import Counter

# Standard integer mass table for 20 amino acids.
amino_acid_mass = {
    'G': 57,  'A': 71,  'S': 87,  'P': 97,  'V': 99,
    'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
    'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
    'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186,
}


def linear_spectrum(peptide):
    """List of all contiguous (non-wrapping) linear subpeptide masses, plus 0."""
    n = len(peptide)
    prefix_mass = [0] * (n + 1)
    for i in range(1, n + 1):
        amino_acid = peptide[i - 1]
        prefix_mass[i] = prefix_mass[i - 1] + amino_acid_mass[amino_acid]

    spectrum = [0]
    for i in range(n):
        for j in range(i + 1, n + 1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])
    return spectrum


theoretical_spectrum = linear_spectrum("PEEP")

def read_experimental_spectrum(file):
    """Reads a single line of space separated integer masses."""
    with open(file, "r") as f:
        return [int(x) for x in f.readline().split()]

experimental_spectrum=read_experimental_spectrum("antibiotics_week_3/spectrum.txt")

def linear_score():
    """Count the number of matching masses. Higher score means higher resemblance."""
    match=0

    for mass1, count1 in Counter(theoretical_spectrum).items():
        for mass2, count2 in Counter(experimental_spectrum).items():
            if mass1==mass2:
                match+=min(count1,count2) #if masses are the same in both spectra, take the lowest count and add to score

    return match

print(linear_score())