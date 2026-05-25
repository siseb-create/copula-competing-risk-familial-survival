############################################################
## Kinship matrix and kinship extraction functions
##
## This file defines the family kinship matrix (2*kinship) and provides
## functions to extract the kinship coefficients required for
## non-proband and bivariate likelihood contributions.
############################################################

kinship <- matrix(
  c(1,1/2,1/8,1/8,1/8,1/8,
    1/2,1,1/8,1/8,1/8,1/8,
    1/8,1/8,1,1/2,1/8,1/8,
    1/8,1/8,1/2,1,1/8,1/8,
    1/8,1/8,1/8,1/8,1,1/2,
    1/8,1/8,1/8,1/8,1/2,1),
  ncol=6)
kin <- kinship[1,2:6]

kinJK <- cbind(kinship[1, combn(2:6, 2)[1, ]],
               kinship[1, combn(2:6, 2)[2, ]],
               mapply(function(j, k) kinship[j, k], combn(2:6, 2)[1, ], combn(2:6, 2)[2, ])
)
