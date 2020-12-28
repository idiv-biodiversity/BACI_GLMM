# Script to download the data from zenodo.org 

dir.create(path = "./Data")

download.file(url = "https://zenodo.org/record/4384785/files/baci_predation_data.txt?download=1",
              destfile = "./Data/baci_predation_data.txt")

download.file(url = "https://zenodo.org/record/4384785/files/metadata.txt?download=1",
              destfile = "./Data/metadata.txt")


# Rename columns and the data file to the names that were used when we published
# the R code for the analysis in the "Appendix S1: R code used in statistical
# analysis" at
# https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Frec.12678&file=rec12678-sup-0001-AppendixS1.pdf

BACI.dt <- read.table(file = "./Data/baci_predation_data.txt", header = TRUE, sep = "\t")

names(BACI.dt) <- c("obs", "Site_F", "Year_F", "Period.BA_F", "SiteClass.CI_F", "Method_F", "CL", "NOTAB")

write.csv(x = BACI.dt, file = "./Data/BACI.dt.csv", row.names = FALSE)
