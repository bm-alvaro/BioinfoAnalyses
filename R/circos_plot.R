# Alvaro Benitez, 2023

library(readxl)
library(canvasXpress)

data_ <- read_excel("Freqs/compilation.xlsx")
data_ <- as.data.frame(data_)
rownames(data_) <- c('Naturales', 'Artificiales', 'EMS')
data_ <- data_[,-1]

datx_ <- read_excel("Freqs/compilation.xlsx", sheet = "datx_")
datx_ <- as.data.frame(datx_)
datx_[1,] <- substr(datx_[1,], 1, nchar(datx_[1,])-1)

add_space_after_third_char <- function(x) {
  if (nchar(x) >= 3) {
    paste0(substr(x, 1, 3), " ", substr(x, 4, nchar(x)))
  } else {
    x
  }
}

datx_[1, ] <- lapply(datx_[1, ], add_space_after_third_char)

rownames(datx_) <- c('Chromosomes')

vars = c("Naturales",
         "Artificiales",
         "EMS")
valz = c(1,
         2,
         3)
smpz = c("Ring")
datz_ = as.data.frame(matrix(valz, nrow = 3, ncol = 1, byrow = FALSE, dimnames = list(vars, smpz)))
canvasXpress(
  data = data_,
  smpAnnot = datx_,
  varAnnot = datz_,
  "graphType" = "Circular",
  "legendKeyBackgroundBorderColor" = "rgba(255,255,255,0)",
  "legendKeyBackgroundColor" = "rgba(255,255,255,0)",
  "ringGraphType" = list("bar","bar","bar"),
  "ringGraphWeight" = list(25,25,25),
  "segregateSamplesBy" = list("Chromosomes"),
  "segregateVariablesBy" = list("Ring"),
  "showTransition" = FALSE,
  "smpOverlays" = list("Chromosomes"),
  "title" = "Circos Plot TILLING",
  "transitionStep" = 50,
  "transitionTime" = 1500,
  "graphOrientation" = "horizontal",
  "arcSegmentsSeparation" = 1,
  "showSampleNames" = FALSE,
  "circularCenterProportion" = 0.3,
  "ringSeparation" = 20,
  "axisAlgorithm" = "heckbert"
)

