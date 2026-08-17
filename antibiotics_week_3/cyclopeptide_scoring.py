"""
CyclopeptideScoring
========================
Compare an experimental spectrum with the theoretical cyclic spectrum of the peptide. Count the number of matching masses.
This is important when the experimental spectrum is not ideal and have erroneous or missing masses.
"""
from collections import Counter

# Standard integer mass table for 20 amino acids.
# Some amino acids share the same integer mass, so the Cyclic Spectrum can't always uniquely determine the peptide.

amino_acid_mass = {
    'G': 57,  'A': 71,  'S': 87,  'P': 97,  'V': 99,
    'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
    'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
    'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186,
}


def cyclic_spectrum(peptide):
    """
    Returns the theoretical (cyclic) spectrum of peptide (including 0 and total mass).
    """
    n = len(peptide)

    # PrefixMass[i] = total mass of the first i amino acids (PrefixMass[0] = 0)
    prefix_mass = [0] * (n + 1)
    for i in range(1, n + 1):
        amino_acid = peptide[i - 1]
        prefix_mass[i] = prefix_mass[i - 1] + amino_acid_mass[amino_acid]

    total_mass = prefix_mass[n]  # mass of the entire peptide

    spectrum = [0]  # the empty subpeptide always has mass 0

    for i in range(n):
        for j in range(i + 1, n + 1):
            subpeptide_mass = prefix_mass[j] - prefix_mass[i]
            spectrum.append(subpeptide_mass)

            # only "interior" subpeptides (not already touching either
            # end) have a genuine wrap-around counterpart -- a
            # subpeptide touching position 0 or n is already the full
            # remaining stretch, so its "complement" would be empty/full
            # and is already covered elsewhere.
            if i > 0 and j < n:
                spectrum.append(total_mass - subpeptide_mass)

    return spectrum


theoretical_spectrum = cyclic_spectrum("MAMA")

def read_experimental_spectrum(file):
    """Reads a single line of space separated integer masses."""
    with open(file, "r") as f:
        return [int(x) for x in f.readline().split()]

experimental_spectrum=read_experimental_spectrum("antibiotics_week_3/spectrum.txt")

def score():
    """Count the number of matching masses. Higher score means higher resemblance."""
    match=0

    for mass1, count1 in Counter(theoretical_spectrum).items():
        for mass2, count2 in Counter(experimental_spectrum).items():
            if mass1==mass2:
                match+=min(count1,count2) #if masses are the same in both spectra, take the lowest count and add to score

    return match

print(score())