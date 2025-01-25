set.seed(1)

vb <- function(t, L1, L2, k, t1, t2)
  L1 + (L2-L1) * (1-exp(-k*(t-t1))) / (1-exp(-k*(t2-t1)))
inv_vb <- function(L, L1, L2, k, t1, t2)
  t1 - 1/k * log(1 - (L-L1) * (1-exp(-k*(t2-t1))) / (L2-L1))

# True parameter values
L1 <- 30
L2 <- 80
k <- 0.34
sigma_min <- 1.4
sigma_max <- 1.8
t1 <- 0
t2 <- 4

# Simulate otoliths
Noto <- 32
Aoto <- exp(rnorm(Noto, m=-1, s=0.3))
Loto <- vb(Aoto, L1, L2, k, t1, t2)
Loto_obs <- Loto + rnorm(Noto, m=0, s=1.5)

# Simulate tags
Ntag <- 487
liberty <- exp(rnorm(Ntag, m=-1, s=0.5))
lenRel <- sqrt(rnorm(Ntag, m=2500, s=500))
ageRel <- inv_vb(lenRel, L1, L2, k, t1, t2)
ageRec <- ageRel + liberty
lenRec <- vb(ageRec, L1, L2, k, t1, t2)
L_min <- min(Loto, lenRel, lenRec) - 2*sigma_min
L_max <- max(Loto, lenRel, lenRec) + 2*sigma_max
sigma_slope <- (sigma_max - sigma_min) / (L_max - L_min)  # s <- a + b*age
sigma_intercept <- sigma_min - L_min * sigma_slope
sigma_lenRel <- sigma_intercept + sigma_slope*lenRel
sigma_lenRec <- sigma_intercept + sigma_slope*lenRec
lenRel_obs <- lenRel + rnorm(Ntag, m=0, s=sigma_lenRel)
lenRec_obs <- lenRec + rnorm(Ntag, m=0, s=sigma_lenRec)

# Simulated observations
otoliths <- data.frame(age=round(Aoto, 2), len=round(Loto_obs, 1))
tags <- data.frame(lenRel=round(lenRel_obs, 1), lenRec=round(lenRec_obs, 1),
                   liberty=round(liberty, 2))

write.table(otoliths, "otoliths_skj.tab", quote=FALSE, sep="\t",
            row.names=FALSE)
write.table(tags, "tags_skj.tab", quote=FALSE, sep="\t", row.names=FALSE)
