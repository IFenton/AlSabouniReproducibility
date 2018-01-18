# IF reanalysis of Nadia's repeatability
# IF reanalysis
# Isabel Fenton
# Date created: 18 / 1 / 2017
# Date last edited: 18 / 1 / 2017
# 
# Re running the analyses done by Nadia in the repeatability analysis. Lab book repeatability
# 
# Previous file: N/A
# Next file:
#   

rm(list = ls())
# Inputs
# Outputs
# Source files / libraries
library("readxl")

# 1. Load in the data -----------------------------------------------------

# 1a. Load the datasets ---------------------------------------------------
# the original IDs
slide125 <- as.data.frame(read_excel("Data/PersonIDs.xlsx", sheet = "125slide"), stringsAsFactors = FALSE)
slide150 <- as.data.frame(read_excel("Data/PersonIDs.xlsx", sheet = "150slide"), stringsAsFactors = FALSE)
digital125 <- as.data.frame(read_excel("Data/PersonIDs.xlsx", sheet = "125digital"), stringsAsFactors = FALSE)
digital150 <- as.data.frame(read_excel("Data/PersonIDs.xlsx", sheet = "150digital"), stringsAsFactors = FALSE)

datasets <- list(slide125, slide150, digital125, digital150)

# the species abbreviations
sp.abb <- as.data.frame(read_excel("Data/PersonIDs.xlsx", sheet = "Abbreviations"))

# people metadata
people <- as.data.frame(read_excel("Data/PeopleMetadata.xlsx"))

# specimen size
size125 <- as.data.frame(read_excel("Data/SpecimenSize.xlsx", sheet = "Size125"))
size150 <- as.data.frame(read_excel("Data/SpecimenSize.xlsx", sheet = "Size150"))

# 1b. Sanity check --------------------------------------------------------
str(slide125)

# check that all the levels are correct and match to those in the abbreviations
for (i in datasets) {
  tmp <- sort(unique(unlist(sapply(i[, 2:ncol(i)], unique))))
  print(tmp %in% sp.abb$Abbreviation) 
}
# all of those are true, so all the abbreviations match to real species

# get the digital pair and the slide pair are the same size
dim(digital125) == dim(digital150)
dim(slide125) == dim(slide150)

# 2. Calculate the consensus ----------------------------------------------

# 2a. Consensus 50 --------------------------------------------------------
# calculate the cutoff for consensus 50
c50_cutoff_slide <- round(ncol(slide125[, nchar(names(slide125)) < 5])/2, 0)
c50_cutoff_digital <- round(ncol(digital125[, nchar(names(digital125)) < 5])/2, 0)

# calculate the consensus values
slide125$IFc50 <- apply(slide125[, nchar(names(slide125)) < 5], 1, function (x) ifelse(length(which(table(x) > c50_cutoff_slide)) > 0, names(table(x))[which(table(x) > c50_cutoff_slide)], "nc")) 
# if one name has the majority (i.e. the length of the table > 8 is 1), then return that name
# if the there are no names that have more than 8 (given there are 17 IDs, then consensus-50 needs at least 9 to match), then return "nc"
slide150$IFc50 <- apply(slide150[, nchar(names(slide150)) < 5], 1, function (x) ifelse(length(which(table(x) > c50_cutoff_slide)) > 0, names(table(x))[which(table(x) > c50_cutoff_slide)], "nc")) 
digital125$IFc50 <- apply(digital125[, nchar(names(digital125)) < 5], 1, function (x) ifelse(length(which(table(x) > c50_cutoff_digital)) > 0, names(table(x))[which(table(x) > c50_cutoff_digital)], "nc")) 
digital150$IFc50 <- apply(digital150[, nchar(names(digital150)) < 5], 1, function (x) ifelse(length(which(table(x) > c50_cutoff_digital)) > 0, names(table(x))[which(table(x) > c50_cutoff_digital)], "nc")) 


# see if this matches Nadia's results
write.csv(slide125[which(slide125$consensus50 != slide125$IFc50), ], "Outputs/Slide125Cmismatch.csv")
slide150[which(slide150$consensus50 != slide150$IFc50), ]
digital125[which(digital125$consensus50 != digital125$IFc50), ] # 2 errors in Nadia's results
digital150[which(digital150$consensus50 != digital150$IFc50), ] # matches for everything

# does removing in NAs help?
tmp <- apply(slide125[, nchar(names(slide125)) < 5], 1, function (x) ifelse(length(which(table(x) > (c50_cutoff_slide - sum(x == "na")))) > 0, names(table(x))[which(table(x) > (c50_cutoff_slide - sum(x == "na")))], "nc")) 
write.csv(slide125[which(slide125$consensus50 != tmp), ], "Outputs/Slide125CmismatchNA.csv")


# 2b. Consensus 20 --------------------------------------------------------

# Supp Table 1

# Supp Table 3

# 3. Agreement between workers (pairwise comparisons) ---------------------


# 3a. C20_score_participant -----------------------------------------------


# 3b. C20_score_group -----------------------------------------------------

# 3c. Average pairwise agreement scores -----------------------------------
# Figure 3
# Table 4

# 3d. Lumped pairwise agreement scores ------------------------------------


# 3d. Radial plots --------------------------------------------------------
# Figure 4

# 3e. Confusion matrix ----------------------------------------------------





# 4. NMDS -----------------------------------------------------------------
# Figure 2


# 5. Repeated analysis by workers -----------------------------------------
# ex figure 6

# 6. Digital vs. slides ---------------------------------------------------
# Table 6

# 7. SST ------------------------------------------------------------------
# Figure 6

# 8. Diversity ------------------------------------------------------------
# Table 5

# Figure 7

# Figure 8


# 9. Outliers -------------------------------------------------------------
# Table 7
# Table 8 
# 10. Size vs. maximum agreement ------------------------------------------
# Figure 5
# ex. figure 7

# 11. Comparison of different tests
# ex. Figure 12
