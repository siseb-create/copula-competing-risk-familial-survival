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
