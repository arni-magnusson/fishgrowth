set.seed(1)

vb <- function(t, L1, L2, k, t1, t2)
  L1 + (L2-L1) * (1-exp(-k*(t-t1))) / (1-exp(-k*(t2-t1)))
inv_vb <- function(L, L1, L2, k, t1, t2)
  t1 - 1/k * log(1 - (L-L1) * (1-exp(-k*(t2-t1))) / (L2-L1))

# True parameter values
L1 <- 30
L2 <- 80
k <- 0.34
sigma_1 <- 1.4
sigma_2 <- 1.8
t1 <- 0
t2 <- 4
L_short <- 30
L_long <- 60

# Simulate otoliths
Noto <- 32
Aoto <- exp(rnorm(Noto, m=-1, s=0.3))
Loto <- vb(Aoto, L1, L2, k, t1, t2)
Loto_obs <- Loto + rnorm(Noto, m=0, s=1.5)

# Simulate tags
# (generate 2*Ntag datapoints initially, select final Ntag later)
Ntag <- 487
liberty <- exp(rnorm(2*Ntag, m=-1, s=0.5))
lenRel <- sqrt(rnorm(2*Ntag, m=2500, s=500))
ageRel <- inv_vb(lenRel, L1, L2, k, t1, t2)
ageRec <- ageRel + liberty
lenRec <- vb(ageRec, L1, L2, k, t1, t2)
sigma_slope <- (sigma_2 - sigma_1) / (L_long - L_short)  # s <- a + b*age
sigma_intercept <- sigma_1 - L_short * sigma_slope
sigma_lenRel <- sigma_intercept + sigma_slope*lenRel
sigma_lenRec <- sigma_intercept + sigma_slope*lenRec
lenRel_obs <- lenRel + rnorm(2*Ntag, m=0, s=sigma_lenRel)
lenRec_obs <- lenRec + rnorm(2*Ntag, m=0, s=sigma_lenRec)

# Simulated observations
otoliths <- data.frame(age=round(Aoto, 2), len=round(Loto_obs, 1))
tags <- data.frame(lenRel=round(lenRel_obs, 1), lenRec=round(lenRec_obs, 1),
                   liberty=round(liberty, 2))

# Apply data selection
dailyGrowth_obs <- (lenRec_obs - lenRel_obs) / 365
tags <- subset(tags, liberty >= 30/365      &  # min 1 month at liberty
                     liberty <= 3           &  # max 3 years at liberty
                     dailyGrowth_obs >= 0   &  # zero or positive growth
                     dailyGrowth_obs <= 0.2 &  # max 2 mm daily growth
                     lenRel >= 30)             # min 30 cm at release
tags <- tags[1:Ntag,]

write.table(otoliths, "otoliths_sim.tab", quote=FALSE, sep="\t",
            row.names=FALSE)
write.table(tags, "tags_sim.tab", quote=FALSE, sep="\t", row.names=FALSE)
