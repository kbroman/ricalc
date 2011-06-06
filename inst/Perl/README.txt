
Perl code for calculating the three-point probabilities in 4-way RILs
by sibling mating, for the autosome or the X chromosome.

gtypes.pl              Determines the set of all parental types and
                       produces a lookup table giving all possible
                       parental types and the corresponding equivalence
                       class for the reduced set after accounting for
                       symmetries.

calc_abs_probs_ver1.pl Actually calculates the 3-point probabilities,
                       by brute-force; requires the results of
                       gtypes.pl.  This version does the calculation
                       for each recombination fraction one at a time
                       (sequentially).

calc_abs_probs_ver2.pl Actually calculates the 3-point probabilities,
                       by brute-force; requires the results of
                       gtypes.pl.  This version does the calculation
                       for all recombination fractions at once, which
                       can be faster, but requires more memory.
