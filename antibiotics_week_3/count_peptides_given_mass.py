AMINO_ACID_MASSES = [

        57, 71, 87, 97, 99, 101, 103, 113, 114,

        115, 128, 129, 131, 137, 147, 156, 163, 186

    ] # 18 amino acid massses (I,L and K,Q have the same masses) So, we would count AIKD and ALQD as just one peptide, not two.

 

def count_peptides_with_mass(target_mass):
#dynamic programming
    dp = [0] * (target_mass + 1)

    dp[0] = 1  # base case: one way to make mass 0 (empty peptide)

    for mass in range(1, target_mass + 1):

        for a in AMINO_ACID_MASSES:

            if mass - a >= 0:

                dp[mass] += dp[mass - a]

 

    return dp[target_mass]

 

print(count_peptides_with_mass(1265)) # pass the given/target mass
#The number of linear peptides having integer mass (m)