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

# Simulate otoliths
Noto <- 32
Aoto <- exp(rnorm(Noto, m=-1, s=0.3))
Loto <- vb(Aoto, L1, L2, k, t1, t2)
Loto.obs <- Loto + rnorm(Noto, m=0, s=1.5)

# Simulate tags
# (generate 2*Ntag datapoints initially, select final Ntag later)
Ntag <- 487
liberty <- exp(rnorm(2*Ntag, m=-1, s=0.5))
lenRel <- sqrt(rnorm(2*Ntag, m=2500, s=500))
ageRel <- inv_vb(lenRel, L1, L2, k, t1, t2)
ageRec <- ageRel + liberty
lenRec <- vb(ageRec, L1, L2, k, t1, t2)
lenRel.obs <- lenRel + rnorm(2*Ntag, m=0, s=1.7)
lenRec.obs <- lenRec + rnorm(2*Ntag, m=0, s=1.7)

# Simulated observations
otoliths <- data.frame(age=round(Aoto, 2), len=round(Loto.obs, 1))
tags <- data.frame(lenRel=round(lenRel.obs, 1), lenRec=round(lenRec.obs, 1),
                   liberty=round(liberty, 2))

# Apply data selection
dailyGrowth.obs <- (lenRec.obs - lenRel.obs) / 365
tags <- subset(tags, liberty >= 30/365      &  # min 1 month at liberty
                     liberty <= 3           &  # max 3 years at liberty
                     dailyGrowth.obs >= 0   &  # zero or positive growth
                     dailyGrowth.obs <= 0.2 &  # max 2 mm daily growth
                     lenRel >= 30)             # min 30 cm at release
tags <- tags[1:Ntag,]

write.table(otoliths, "otoliths.skj.tab", quote=FALSE, sep="\t",
            row.names=FALSE)
write.table(tags, "tags.skj.tab", quote=FALSE, sep="\t", row.names=FALSE)
