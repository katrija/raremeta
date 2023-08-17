# code to prepare `Yang 2019` dataset goes here

dat.yang2019 <- data.frame(
  study = 1:16,
  ai = c(1,6,1,14,53,3,2,13,6,0,3,0,5,24,1,4),
  ci = c(0,11,1,7,18,4,2,5,2,0,4,0,5,3,0,0),
  n1i = c(6,36,12,47,94,21,19,95,21,17,78,20,36,95,10,32),
  n2i = c(6,31,23,26,34,21,17,43,23,17,39,57,14,30,9,23)
)


usethis::use_data(dat.yang2019, overwrite = TRUE)
