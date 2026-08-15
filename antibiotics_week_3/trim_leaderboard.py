"""
TrimLeaderboard
========================
Implement linear peptide scoring to a collection of growing peptides (Leaderboard) and update the leaderboard
with just highest scoring N peptides with respect to experimental spectrum
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


def read_experimental_spectrum(file):
    """Reads a single line of space separated integer masses."""
    with open(file, "r") as f:
        return [int(x) for x in f.readline().split()]


def linear_score(theoretical_spectrum,experimental_spectrum):
    """Count the number of matching masses. Higher score means higher resemblance."""
    match=0

    for mass1, count1 in Counter(theoretical_spectrum).items():
        for mass2, count2 in Counter(experimental_spectrum).items():
            if mass1==mass2:
                match+=min(count1,count2) #if masses are the same in both spectra, take the lowest count and add to score

    return match

def trim(N):
    experimental_spectrum=read_experimental_spectrum("antibiotics_week_3/spectrum.txt")
    linear_scores={}
    peptides="MMDRSSVCGDSMCIEWKLIIPQANHLFPWDHYQNSCTQYICGISNTDN WYKMSGYSMWCVDNVSIMMYWLIKNTAKQYQHADETWCIQWRWGNKAC SKTFWFLQYHFCVAYSMFQRGHETHIATKWAATKCMHDIAWVDMTTFF TDMKHRYGAQEAVQCDRSCLEYVPGTVFLTLRHTIKRVRCDCYICGSF EAELHVSNGKIQARNDDDMTCRAFAPQAQESQCRLDVLGTYWHQQINM HQWWESRYIENSQYLNGVKQWWPSKILASYFRSNESWNYYELYRPRVY QSMLDGGGRFKLNLCTWCGADRRHMISNAPWNNTDMRHVVFLTLKMNQ MLFCADDGYIDVIAGFVADLNIWNSKIDEKSLECIYVDTFHNNISWKR LIWVPEDKLYWVPFVSTVSTCPVRYNDECWIINPCQYLINIKDCLIEY DQRNSNEKRPSFCWQCDMDIDRETKKPTIEKYKRFDWKGHDELSPEMY CAKLNDKGRMYMHLNPCQDTWPLCFVEEAIIMPEEEEYRYLWMWPIFA SESRASYAVLQKGANMNRYFCDFWAYKCWIICTLFPFYDKVWPWCHPK TNKMHVYNTMNDACPIWQFLTPSCKSDTMFPCQAHFRTSFGRQAAWRF IYPTELLEYFIAVQRTYGHLNYPVYAYWPGATGYEYNGQWRQRGNLAA KNLWIILWDHYPDTMWERMRVCFHCPYQCTKDWQSAGRMHNNMEVYWG QLGEDEPIVRAHKMPVFLGLYAPSTMINLSLCYHLEVYQMFTRENHCE ILRVMLRFAFWHHNGPCGARYNLDTPHKSMPAWQPIKDPGGGAHPVEK EMFDNTYMFGQINKTCVDWPIEGHRCLVNMEQTETWCFTSKKLFDTYA ETRCECVIKANEAENYRKEKIDIEYLIRISAELKIPFVYQLYTMDEDN WPDEEIVGHPDDHRWQQHKWHDETWQTWSPRRPNRTVCNPFNFLITWF"
    leaderboard=peptides.split() # a list of space separated peptides

    for peptide in leaderboard:
        theoretical_spectrum = linear_spectrum(peptide)
        linear_scores[peptide]=linear_score(theoretical_spectrum,experimental_spectrum)

    sorted_scores = sorted(linear_scores.items(), key=lambda x: x[1], reverse=True) #gives a sorted list of tuples based on highest to lowest score -> [("PEPTIDE2", 32), ("PEPTIDE1", 15), ("PEPTIDE3", 8)]
    (x,y)=sorted_scores[N-1] # N-1 since list indexing starts counting from 0
    trimmed_leaderboard=[a for (a,b) in sorted_scores if b>=y]
    leaderboard=trimmed_leaderboard
    return leaderboard

print(*trim(6))