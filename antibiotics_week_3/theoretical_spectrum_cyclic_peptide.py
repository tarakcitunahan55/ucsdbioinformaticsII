"""
CyclicSpectrum
===============
Computes the theoretical spectrum of a cyclic peptide: the sorted list
of masses of every contiguous subpeptide, including subpeptides that
wrap around the end of the string back to the start (since the peptide
is a loop, not a straight chain), plus mass 0 and the mass of the
whole peptide.

Pipeline:
- Build PrefixMass exactly as in LinearSpectrum: PrefixMass[i] = total
  mass of the first i amino acids, PrefixMass[0] = 0.
- Any non-wrapping subpeptide from i to j has mass
  PrefixMass[j] - PrefixMass[i], same as the linear spectrum.
- For cyclic peptides: for every interior subpeptide (one that
  doesn't already touch position 0 or the very end), there is also a
  wrap-around subpeptide covering everything else -- and because the
  peptide is a cycle, that leftover stretch is itself a valid
  contiguous subpeptide (it just continues past the end back to the
  start). Its mass: total_mass - (PrefixMass[j] - PrefixMass[i])
  This lets us account for every wrap-around subpeptide with one
  extra subtraction, without ever physically wrapping the string.
- Note: Cyclic peptides do not have a free N and C terminus end. Do not add H20 (18 Da) mass
to the resulting subpeptides, not as a simplification like in linear spectrum, but as correct biochemistry. 
"""

# Standard integer mass table for 20 amino acids.
# Some amino acids share the same integer mass, so the LinearSpectrum can't always uniquely determine the peptide.
AMINO_ACID_MASS = {
    'G': 57,  'A': 71,  'S': 87,  'P': 97,  'V': 99,
    'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114,
    'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131,
    'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186,
}


def cyclic_spectrum(peptide, amino_acid_mass=AMINO_ACID_MASS):
    """
    Returns the sorted theoretical (cyclic) spectrum of peptide.
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

    return sorted(spectrum)


def read_peptide(file):
    """Reads a single peptide string from the first line of the file."""
    with open(file, "r") as f:
        return f.readline().strip()


if __name__ == "__main__":
    peptide = read_peptide("antibiotics_week_3/spectrum.txt")
    spectrum = cyclic_spectrum(peptide)
    print(*spectrum)