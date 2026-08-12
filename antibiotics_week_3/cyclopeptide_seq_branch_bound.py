"""
CyclopeptideSequencing
========================
Reconstructs a cyclic peptide from its theoretical (cyclic) spectrum
using branch-and-bound (brute force algorithms that enumerate all candidate solutions but discard large subsets of hopeless candidates).

Peptides here are represented as lists of integer masses (not amino
acid letters), since several amino acids share the same mass and are
therefore indistinguishable from the spectrum alone.

Branch-and-bound loop:
- Branch: extend every candidate peptide by each of the 18 possible
  amino acid masses (Expand) (one caveat: we extend the peptides by all amino acids, not just the ones present in spectrum
  -> unneccsarily longer list to trim (bound) later)
- Bound: for each extended candidate,
    - if its mass equals the spectrum's parent mass, check whether its
      full cyclic spectrum matches the given spectrum exactly -- if so
      it's a valid answer; either way it's removed from further growth
      (it can't get longer).
    - otherwise, discard it unless its linear spectrum is still a
      sub-multiset of the given spectrum ("consistent") -- inconsistent
      candidates can never grow into a valid answer, so cutting them
      early is what keeps this tractable
"""

from collections import Counter

AMINO_ACID_MASSES = [
    57, 71, 87, 97, 99, 101, 103, 113, 114,
    115, 128, 129, 131, 137, 147, 156, 163, 186
]


def peptide_mass(peptide):
    """Total mass of a peptide."""
    return sum(peptide)


def linear_spectrum(peptide):
    """Sorted list of all contiguous (non-wrapping) linear subpeptide masses, plus 0."""
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
    """Sorted list of all subpeptide masses including wrap-around ones, plus 0 and total mass."""
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
            if i > 0 and j < n:  # interior subpeptide -> has a wrap-around complement
                spectrum.append(total_mass - sub_mass)
    return spectrum


def expand(peptides):
    """Branching step: extend every candidate peptide by each possible 18 amino acid mass."""
    return [peptide + [a] for peptide in peptides for a in AMINO_ACID_MASSES] #+ on two lists means concatenation — it joins them end-to-end into a brand new list


def is_consistent(peptide, spectrum_counts):
    """
    Bounding check: peptide survives only if every mass in its linear spectrum occurs no more often than in the target spectrum.
    """
    for mass, count in Counter(linear_spectrum(peptide)).items():
        if count > spectrum_counts.get(mass, 0):
            return False
    return True


def cyclopeptide_sequencing(spectrum):
    """Returns all peptides (as lists of masses) whose cyclic spectrum matches given spectrum."""
    spectrum_counts = Counter(spectrum)
    parent_mass = max(spectrum)

    candidate_peptides = [[]]     # start from the empty peptide
    final_peptides = []

    while candidate_peptides: #[[]] is True since an empty list is the 0th element of the outer list
        candidate_peptides = expand(candidate_peptides)
        survivors = []

        for peptide in candidate_peptides:
            if peptide_mass(peptide) == parent_mass:
                if Counter(cyclic_spectrum(peptide)) == spectrum_counts:
                    if peptide not in final_peptides:
                        final_peptides.append(peptide)

            elif is_consistent(peptide, spectrum_counts):
                survivors.append(peptide)
            # else: inconsistent -> not appended to survivors, so discarded (bounding step)
 #peptide whose mass equals the parent mass is never appended to survivors (cannot be extended anymore)
        candidate_peptides = survivors #new candidates to be extended

    return final_peptides


def read_spectrum(file):
    """Reads a single line of space separated integer masses."""
    with open(file, "r") as f:
        return [int(x) for x in f.readline().split()]


if __name__ == "__main__":
    spectrum = read_spectrum("antibiotics_week_3/spectrum.txt")
    peptides = cyclopeptide_sequencing(spectrum)
    print(*("-".join(map(str, p)) for p in peptides)) #make each element of peptide list p a string since ".join" requires strings