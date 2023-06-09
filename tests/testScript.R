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

###
#1. INTRODUCTING CONTINUITY CORRECTIN AT DIFFERENT STAGES OF THE WORKFLOW
###

# we compare different ways of fitting the same model on the same dataset and compare results

## random effects meta analysis with REML-heterogeneity estimator and
## constant continuity correction (adding 0.5) to all 0 studies

# 1: using the rareIV function directly (continuity correction happens inside)
fit1 <- rareIV(ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.anand1999,
              measure="logOR", cc = "constant", ccval = 0.5, ccto = "only0", method="FE")

# 2: using rareIV on previously continuity corrected data
d2 <- dat.anand1999[,c("ai", "n1i", "ci", "n2i")]

# dropping 00-studies
dzstudies <- which((d2$ai == 0 & d2$ci == 0 ) | (d2$ci == d2$n2i & d2$ai == d2$n1i))
d2 <- d2[-dzstudies,]

# continuity correction
ccstudies <- which(d2$ai == 0 | d2$ai == d2$n1i | d2$ci == 0 | d2$ci == d2$n2i)
d2[ccstudies,c("ai", "ci")] = d2[ccstudies, c("ai", "ci")] + 0.5
d2[ccstudies,c("n1i", "n2i")] = d2[ccstudies, c("n1i", "n2i")] + 1


fit2 <- rareIV(ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=d2,
               measure="logOR", method="FE")

# 3: using rareIV on pre-processed data (rareData object) with usage of rareCC by hand
x3 <- rareDescribe(ai=ai, n1i=n1i, ci=ci, n2i=n2i, data=dat.anand1999)
x3 <- rareCC(x3, cc = "constant", ccval = 0.5, ccto = "only0")

fit3 <- rareIV(x3, measure = "logOR", method = "FE")

# comparison

test_that("rareIV works the same with the application of the continuity correction at different stages", {
  expect_identical(
    fit1[c("model","b","se","zval","pval","ci.lb","ci.ub","vb","tau2","se.tau2")],
    fit2[c("model","b","se","zval","pval","ci.lb","ci.ub","vb","tau2","se.tau2")]
  )

  expect_identical(
    fit2[c("model","b","se","zval","pval","ci.lb","ci.ub","vb","tau2","se.tau2")],
    fit3[c("model","b","se","zval","pval","ci.lb","ci.ub","vb","tau2","se.tau2")]
 )
})



###
#3. TESTING FOR THE SIGN ERROR IN 'rarePeto'
###

library(metafor)

## 3.1
d <- dat.axfors2021
d <- data.frame(ai = d$hcq_arm_event, n1i = d$hcq_arm_total, ci = d$control_arm_event, n2i = d$control_arm_total)

fit1 <- rarePeto(ai = hcq_arm_event, n1i = hcq_arm_total, ci = control_arm_event, n2i = control_arm_total, data = dat.axfors2021)

fit2 <- rma.peto(ai = hcq_arm_event, n1i = hcq_arm_total, ci = control_arm_event, n2i = control_arm_total, data = dat.axfors2021)

# comparison

test_that("rarePeto and rma.peto calculate the same values", {
  expect_identical(
    fit1[c("b","se","zval","pval","ci.lb","ci.ub")],
    fit2[c("b","se","zval","pval","ci.lb","ci.ub")]
  )
})

# introduced as.numeric to prevent integer overflow?

## 3.2 same analysis, different dataset
d <- dat.cannon2006
d <- data.frame(ai = d$ep1t, n1i = d$nt, ci = d$ep1c, n2i = d$nc)

fit1 <- rarePeto(ai = ep1t, n1i = nt, ci = ep1c, n2i = nc, data = dat.cannon2006)

fit2 <- rma.peto(ai = ep1t, n1i = nt, ci = ep1c, n2i = nc, data = dat.cannon2006, to = "none")

# comparison

test_that("rarePeto and rma.peto calculate the same values", {
  expect_identical(
    fit1[c("b","se","zval","pval","ci.lb","ci.ub")],
    fit2[c("b","se","zval","pval","ci.lb","ci.ub")]
  )
})

# has different values because of different way of calculating the products
# ?solution: copy Viechtbauers representation?



