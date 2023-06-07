# testing using datasets from metadat
library("metadat")

# load packages
data(dat.anand1999,
     dat.axfors2021,
     dat.cannon2006,
     dat.collins1985a,
     dat.collins1985b,
     dat.damico2009,
     dat.colditz1994,
     dat.hartmannboyce2018,
     dat.yusuf1985)

## Introducing continuity correction before and during the workflow
fit <- rareIV(ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.anand1999,
              measure="logOR", method="REML", cc="constant")

