usethis::use_data_raw("Nissen2007", overwrite = TRUE)

## code to prepare `Nissen2007` dataset goes here

## COMMENTS
## Inconsistent   : Study 49653/330  1172 -> 1181 in nRosiglitazone; 377 -> 382 in nControl
##                in the first Table of the Paper (as opposed to third one)
##
## Method of coding
## "SulfonylureaX": 1-either glyburide or gliclazide, 2-glyburide, gliclazide or glipizide
##                  3-glyburide, gliclazide, chlorporopamide, glimepiride, or tolbatamide
##                  4-uncpecified,
##                  5-3-glyburide, gliclazide, chlorporopamide, glimepiride, or tolbatamide

nissen2007 <- data.frame(
  study = c(
    "49653/011", "49653/020", "49653/024", "49653/093", "49653/094",
    "100684", "49653/143", "49653/211", "49653/284", "712753/008",
    "AVM100264", "BRL 49653C/185", "BRL 49653/334", "BRL 49653/347",
    "49653/015", "49653/079", "49653/080", "49653/082", "49653/085",
    "49653/095", "49653/097", "49653/125", "49653/127", "49653/128",
    "49653/134", "49653/135", "49653/136", "49653/145", "49653/147",
    "49653/162", "49653/234", "49653/330", "49653/331", "49653/137",
    "SB-712753/002", "SB-712753/003", "SB-712753/007", "SB-712753/009",
    "49653/132", "AVA100193", "DREAM", "ADOPT"
  ),
  nRosiglitazone = c(
    357, 391, 774, 213, 232, 43, 121, 110, 382, 284, 294,
    563, 278, 418, 395, 203, 104, 212, 138, 196, 122, 175,
    56, 39, 561, 116, 148, 231, 89, 168, 116, 1172, 706, 204,
    288, 254, 314, 162, 442, 394, 2635, 1456
  ),
  miRosiglitazone = c(
    2, 2, 1, 0, 1, 0, 1, 5, 1, 1, 0, 2, 2, 2, 2, 1, 1,
    2, 3, 0, 0, 0, 1, 1, 0, 2, 1, 1, 1, 1, 0, 1, 0, 1,
    1, 1, 1, 0, 1, 1, 15, 27
  ),
  cvRosiglitazone = c(
    1, 0, 0, 0, 1, 0, 0, 3, 0, 0, 2, 0, 0, 0, 2, 1, 0, 1, 1, 1,
    0, 0, 0, 0, 1, 2, 2, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 12, 2
  ),
  nControl = c(
    176, 207, 185, 109, 116, 47, 124, 114, 384, 135, 302, 142, 279,
    212, 198, 106, 99, 107, 139, 96, 120, 173, 58, 38, 276, 111, 143,
    242, 88, 172, 61, 377, 325, 185, 280, 272, 154, 160, 112, 124, 2634, 2895
  ),
  miControl = c(
    0, 1, 1, 1, 0, 1, 0, 2, 0, 0, 1, 0, 1, 0, 1, 1, 2, 0, 1, 0, 1, 1, 0, 0, 2, 3,
    0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 9, 41
  ),
  cvControl = c(
    0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 10, 5
  ),
  methodControl = c("Placebo", "Glyburide", "Placebo", "Metformin",
                        "Metformin", "Glyburide", "Glyburide", "Usual care",
                        "Metformin", "Metformin", "Metformin and sulfonylurea1",
                        "Usual care with or without metformin", "Placebo", "Insulin",
                        "Sulfonylurea2", "Glyburide", "Glyburide", "Insulin", "Insulin", "Insulin",
                        "Glyburide", "Sulfonylurea3", "Glyburide", "Placebo", "Placebo", "Glipizide",
                        "Placebo", "Gliclazide", "Sulfonylurea4", "Glyburide", "Glimepiride", "Placebo", "Placebo",
                        "Glyburide and metformin", "Metformin", "Metformin", "Metformin",
                        "Insulin", "Sulfonylurea5", "Placebo", "Placebo", "Metformin and glyburide")
)
