# code to prepare `Hoppen 2023` dataset goes here

dat.hoppen2023 <- data.frame(
  study = 1:7,
  ai = c(12, 5, 12, 5, 11, 0, 10),
  n1i = c(49, 21, 39, 25, 55, 12, 28),
  ci = c(16, 2, 5, 4, 8, 0, 0),
  n2i = c(49, 21, 29, 24, 47, 11, 29)
)


usethis::use_data(dat.hoppen2023, overwrite = TRUE)
