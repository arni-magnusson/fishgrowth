set.seed(1)

vb <- function(t, L1, L2, k, t1, t2)
  L1 + (L2-L1) * (1-exp(-k*(t-t1))) / (1-exp(-k*(t2-t1)))
inv <- function(L, L1, L2, k, t1, t2)
  t1 - 1/k * log(1 - (L-L1) * (1-exp(-k*(t2-t1))) / (L2-L1))

# True parameters
L1 <- 30
L2 <- 80
k <- 0.34
sigma_1 <- 1.4
sigma_2 <- 1.8
t1 <- 0
t2 <- 4

# Simulate otoliths
Noto <- 32
Aoto <- exp(rnorm(Noto, m=-1, s=0.3))
Loto <- vb(Aoto, L1, L2, k, t1, t2)
Loto.obs <- Loto + rnorm(Noto, m=0, s=1.5)

# Simulate tags
Ntag <- 487
liberty <- exp(rnorm(Ntag, m=-1, s=0.5))
lenRel <- sqrt(rnorm(Ntag, m=2500, s=500))
ageRel <- inv(lenRel, L1, L2, k, t1, t2)
ageRec <- ageRel + liberty
lenRec <- vb(ageRec, L1, L2, k, t1, t2)
lenRel.obs <- lenRel + rnorm(Ntag, m=0, s=1.7)
lenRec.obs <- lenRec + rnorm(Ntag, m=0, s=1.7)

# Simulated observations
otoliths <- data.frame(age=round(Aoto, 2), len=round(Loto.obs, 1))
tags <- data.frame(lenRel=round(lenRel.obs, 1), lenRec=round(lenRec.obs, 1),
                   liberty=round(liberty, 2))
