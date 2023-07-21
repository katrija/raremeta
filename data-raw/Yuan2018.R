# code to prepare `Yuan 2018` dataset goes here

dat.yuan2018 <- data.frame(
  study = c("Ackerson et al. 1998", "Cobham 2012", "Jacob and De Guzman 2016",
            "Lyneham and Rapee 2006", "Rapee et al. 2006", "Rohde et al. 2015",
            "Stice et al. 2010", "Thirlwall et al. 2013"),
  ai = c(3,0,0,9,29,6,4,29),
  ci = c(5,0,0,1,12,8,1,6),
  n1i = c(15,20,15,78,90,128,80,125),
  n2i = c(15,12,15,22,87,124,84,69)
)


usethis::use_data(dat.yuan2018, overwrite = TRUE)
