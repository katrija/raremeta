# code to prepare `Okumura 2014` dataset goes here

dat.okumura2014 <- data.frame(
  study = 1:27,
  ai = c(7,1,0,4,4,2,3,0,20,1,1,15,11,1,32,45,9,0,0,12,3,4,43,3,29,4,0),
  ci = c(1,2,0,4,20,2,2,4,19,1,0,5,9,0,37,57,0,0,2,10,2,3,42,2,6,11,8),
  n1i = c(68,15,15,52,53,29,64,44,61,19,13,35,72,41,46,62,38,31,24,23,12,23,99,14,47,167,48),
  n2i = c(42,16,15,23,40,26,66,40,83,18,12,29,79,39,46,74,36,31,26,22,18,29,100,14,28,170,48)
)


usethis::use_data(dat.okumura2014, overwrite = TRUE)
