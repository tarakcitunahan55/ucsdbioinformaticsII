"""
LinearSpectrum
===============
Computes the theoretical linear spectrum of a peptide: the sorted list
of masses of every possible contiguous subpeptide (including the empty
subpeptide of mass 0 and the full peptide itself).

Pipeline:
- Build a table of the 20 amino acid masses.
- Compute PrefixMass[i] = total mass of the first i amino acids of the
  peptide, using a running sum (PrefixMass[0] = 0 by definition).
- Every contiguous subpeptide from position i to j has mass
  PrefixMass[j] - PrefixMass[i] -- so we just need to take the
  difference between every pair of prefix masses (i < j) rather than
  re-summing amino acids over and over for each subpeptide.
- Collect all these differences, plus the mass 0 for the empty
  subpeptide, and sort the result.
- One Simplification: The amino acid masses correspond to residues forming peptide bonds (aa in a protein), so the corresponding
N-terminus is NH, while C-terminus is CO (H2O lost in peptide bond formation). The linear peptide has 
1 free N-termnius (NH3+) and 1 free C-terminus (COO-). Actually, we need to add 1 mol H20 mass (18 Da) to subpeptides.
"""

# Standard integer mass table for 20 amino acids.
# Some amino acids share the same integer mass, so the LinearSpectrum can't always uniquely determine the peptide.
AMINO_ACID_MASS = {
    'G': 57,  'A': 71,  'S': 87,  'P': 97,  'V': 99,
    'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
    'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
    'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186,
}


def linear_spectrum(peptide, amino_acid_mass=AMINO_ACID_MASS):
    """
    Returns the sorted linear spectrum (list of subpeptide masses) of peptide.
    """
    n = len(peptide)

    # PrefixMass[i] = total mass of the first i amino acids (PrefixMass[0] = 0)
    prefix_mass = [0] * (n + 1)
    for i in range(1, n + 1):
        amino_acid = peptide[i - 1]
        prefix_mass[i] = prefix_mass[i - 1] + amino_acid_mass[amino_acid]

    spectrum = [0]  # the empty subpeptide always has mass 0

    # every subpeptide from position i to j (i < j) has mass PrefixMass[j] - PrefixMass[i]
    for i in range(n):
        for j in range(i + 1, n + 1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])

    return sorted(spectrum) # sorts by numeric int value, smallest to biggest


def read_peptide(file):
    """Reads a single peptide string."""
    with open(file, "r") as f:
        return f.readline().strip()


if __name__ == "__main__":
    peptide = read_peptide("antibiotics_week_3/spectrum.txt")
    spectrum = linear_spectrum(peptide)
    print(*spectrum)