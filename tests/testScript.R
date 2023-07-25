# testing how raremeta fairs with NA entries in studies
data <- data.frame(
  ai = c(0, 3, 2, 0),
  bi = c(20, 18, 15, 19),
  ci = c(1, NA, 0, 0),
  di = c(19, NA, 16, 20),
  n1i = c(20, 21, 17, 19),
  n2i = c(20, NA, 16, 20)
)

x <- rareDescribe(
  ai = ai,
  bi = bi,
  ci = ci,
  di = di,
  n1i = n1i,
  n2i = n2i,
  data = data
)

x_CC <- rareCC(x)
x_CC

x_ES <- rareES(x, measure = "logOR", method = "FE", cc = "constant")
x_ES$yi


#rareCC








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

# Collins1985 example:

dat <- data.frame(
  ai = dat.collins1985a$b.xti,
  n1i = dat.collins1985a$nti,
  ci = dat.collins1985a$b.xci,
  n2i = dat.collins1985a$nci
)

rareIV(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat,
       measure = "logOR", method = "REML")
# missing values warning has to be reformatted

rareIV(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat,
       measure = "logOR", method = "REML",
       cc = "constant", ccto = "only0")
# missing values are not handled well

dat <- dat[-which(is.na(dat$ai)),]
dat

fit1 <- rareIV(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat,
       measure = "logOR", method = "REML",
       cc = "constant", ccto = "only0")

fit2 <- rma(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat,
    measure = "OR", method = "REML", drop00 = TRUE)

fit1$ai
fit1$ai.cc

dat_c <- dat

dat_c[which(dat$ai == 0),] <- dat_c[which(dat$ai == 0),] + 0.1
fit1  <- rareIV(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat_c,
               measure = "logOR", method = "REML",
               cc = "none")
fit1
fit1$ai.cc
fit1$ai

fit2 <- rma(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat_c,
            measure = "OR", method = "REML", drop00 = TRUE)
fit2

fit3 <- rareIV(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat,
               measure = "logOR", method = "REML",
               cc = "constant", cccval = 2, tccval = 3, ccto = "all")
fit3$ai.cc
fit3$ai
fit3$ci.cc
fit3$ci

fit4 <- rareIV(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat,
               measure = "logOR", method = "REML",
               cc = "constant", cccval = -5, tccval = -20, ccto = "all")

dat_cc <- rareCC(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat,
                 measure = "logOR", cc = "empirical", method = "DL")
dat_cc
# Ist der Output fuer dat_cc so passend? => ggf. eher sinnvoll,
# hier Info zur Kontinuitaetskorrektur zu geben
# (welche; wie viele Studien wurden korrigiert)

dat_nocc <-  rareDescribe(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat)

rarePeto(dat_cc)
rarePeto(dat_nocc)

rma.peto(ai = ai, n1i = n1i, ci = ci, n2i = n2i,
         data = dat)

rareMH(dat_cc, measure = "logOR")
rareMH(dat_nocc, measure = "logOR")

rma.mh(ai = ai, n1i = n1i, ci = ci, n2i = n2i,
       data = dat, drop00 = FALSE)

# rarePeto und rareMH nutzen die nicht-korrigierten ais, cis usw. (?)
# => ist i.O., sollte nur irgendwo deutlich gemacht werden

dat0 <- dat
dat0$ci <- 0

rareMH(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat0,
       measure = "logOR")
# Fehlermeldung etwas kryptisch - es wird nicht klar, warum die Daten nicht geeignet sind
# Kontinuitaetskorrekturen sollten auch in rareMH verfuegbar sein (aber nicht der Default)
# (das wird manchmal gemacht)

rarePeto(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat00)

dat00 <- dat0
dat00$ai <- 0

rarePeto(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat00)
# Fehlermeldung tritt erst auf, wenn das Outcome in beiden Gruppen in allen Studien 0 ist
# (sollte auch so sein (?).
# Vorschlag: Fehlermeldung umformulieren zu "All studies are double-zero studies -
# no estimation possible."


# Alternative models -----------------------------------------------------------

# Betabin:
rareBetabin(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat,
            measure = "logOR")

# GLMM, fixed intercepts:
rareGLMM(ai = ai, n1i = n1i, ci = ci, n2i = n2i,
         data = dat, measure = "logOR")

rma.glmm(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat,
         measure = "OR")


# GLMM, random intercept:
rareGLMM(ai = ai, n1i = n1i, ci = ci, n2i = n2i,
         data = dat, measure = "logOR", intercept = "random",
         slope = "random")

rma.glmm(ai = ai, n1i = n1i, ci = ci, n2i = n2i, data = dat,
         measure = "OR", model = "UM.RS")
# results differ - find out why!

rareGLMM(ai = ai, n1i = n1i, ci = ci, n2i = n2i,
         data = dat, measure = "logOR", intercept = "random",
         slope = "random", drop00 = TRUE)
# apparently the results differ because rma.glmm per default drops dz studies
# NB: we currently use a different type of LRT compared to rma.glmm -
# which one makes more sense?

# Hypergeometric-normal model:
rareGLMM(ai = ai, n1i = n1i, ci = ci, n2i = n2i,
         data = dat, measure = "logOR", intercept = "random",
         slope = "random", conditional = TRUE)
# Q: is p calculated correctly (currently, tau2 is not taken into account)
# apparently se is not correct - find out why

rma_nchg <- rma.glmm(ai = ai, n1i = n1i, ci = ci, n2i = n2i,
         data = dat, measure = "OR", model = "CM.EL")
