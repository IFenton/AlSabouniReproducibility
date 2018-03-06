# Code for the repeatability analysis 
# Al-Sabouni, Fenton, Telford & Kucera
# Isabel Fenton
# Project: IF reanalysis
# Date created: 5 / 3 / 2018
# Date last edited: 5 / 3 / 2018
# 
# All the code needed to run the Al-Sabouni et al repeatability paper. Lab book repeatability
# 
# Previous file: Reanalysis_NA_IF.R
# Next file:
# 

# ========================================
# replace "ASFigures" with "Figures"

rm(list = ls())

# Inputs ------------------------------------------------------------------
# these files are assumed to be in a folder called Data
# AlSabouni_PersonIDs.xlsx - the raw data and the species abbreviations
# AlSabouni_PeopleMetadata.xlsx - the person metadata for the participants
# AlSabouni_SpecimenSize.xlsx - the size data for the specimens
# AlSabouni_DiversityTemp.xlsx - the temperature data for the sites
# Siccha_ForCenS.xlsx - the ForCenS data for comparison of diversity ranges

# Outputs -----------------------------------------------------------------
## Figures (assuming a folder called "/Figures")
# CombCon_agreement_fullID.png / SepCon_agreement_fullID.png - the percentage agreement with the different consensus values
# Confusion matrices: CombCon_conf_*.png for the combined consensus and SepCon_conf_*.png for the separate consensus values

# Source files / libraries ------------------------------------------------
library(readxl) # reading in xlsx files
library(caret) # for the confusion matrix
library(colorRamps) # colours
library(RColorBrewer) # colours
library(stringr) # for confusion matrix axes
library(vegan) # for the diversity metrics
# library(cluster) # for the dendrogram

source("Code/Confusion_matrix.R")


# 1. Load in the data -----------------------------------------------------

# 1a. Load the datasets ---------------------------------------------------
# the original IDs
slide125 <- as.data.frame(read_excel("Data/AlSabouni_PersonIDs.xlsx", sheet = "125slide"), stringsAsFactors = FALSE)
slide150 <- as.data.frame(read_excel("Data/AlSabouni_PersonIDs.xlsx", sheet = "150slide"), stringsAsFactors = FALSE)
digital125 <- as.data.frame(read_excel("Data/AlSabouni_PersonIDs.xlsx", sheet = "125digital"), stringsAsFactors = FALSE)
digital150 <- as.data.frame(read_excel("Data/AlSabouni_PersonIDs.xlsx", sheet = "150digital"), stringsAsFactors = FALSE)

full.125 <- merge(slide125, digital125, by = "Specimen")
full.150 <- merge(slide150, digital150, by = "Specimen")

# the species abbreviations
sp.abb <- as.data.frame(read_excel("Data/AlSabouni_PersonIDs.xlsx", sheet = "Abbreviations"))

# people metadata
people <- as.data.frame(read_excel("Data/AlSabouni_PeopleMetadata.xlsx", na = "NA"))

# specimen size
size125 <- as.data.frame(read_excel("Data/AlSabouni_SpecimenSize.xlsx", sheet = "Size125"))
size150 <- as.data.frame(read_excel("Data/AlSabouni_SpecimenSize.xlsx", sheet = "Size150"))

# diversity / temperature 
divTemp <- as.data.frame(read_excel("Data/AlSabouni_DiversityTemp.xlsx", na = "NA"))

# 1b. Sanity check --------------------------------------------------------
str(slide125)

# check that all the levels are correct and match to those in the abbreviations
datasets <- list(slide125, slide150, digital125, digital150)

for (i in datasets) {
  tmp <- sort(unique(unlist(sapply(i[, 2:ncol(i)], unique))))
  print(tmp %in% sp.abb$Abbreviation) 
}
rm(i, datasets, tmp)
# all of those are true, so all the abbreviations match to real species

# check the digital pair and the slide pair are the same size
dim(digital125) == dim(digital150)
dim(slide125) == dim(slide150)

# create a list of the columns containing data
col.nam <- list(c125 = which(names(full.125) != "Specimen"), c150 = which(names(full.125) != "Specimen"), s125 = which(names(full.125) %in% names(slide125) & names(full.125) != "Specimen"), s150 = which(names(full.150) %in% names(slide150) & names(full.150) != "Specimen"), d125 = which(names(full.125) %in% names(digital125) & names(full.125) != "Specimen"), d150 = which(names(full.150) %in% names(digital150) & names(full.150) != "Specimen"))

rm(slide125, digital125, slide150, digital150)

# 2. Calculate the consensus ----------------------------------------------

# 2a. Strict consensus --------------------------------------------------------
# calculate the cutoff for the strict consensus
# for the odd numbers, rounding down gives the cutoff as anything greater than this value.
c50_cutoff <- list()
c50_cutoff$slide <- round(length(col.nam$s125)/2, 0)
c50_cutoff$digital <- round(length(col.nam$d125)/2, 0)
c50_cutoff$full <- round(length(col.nam$c125)/2, 0)

# calculate the strict consensus values using the full datasets (i.e. digital / slide combined)
full.125$cSC50 <- apply(full.125[, col.nam$c125], 1, function (x) ifelse(length(which(table(x) > c50_cutoff$full)) > 0, names(table(x))[which(table(x) > c50_cutoff$full)], "nc"))
# if one name has the majority (i.e. the length of the table > c50_cutoff$full is 1), then return that name
# if the there are no names that have more than c50_cutoff$full (given there are 17 IDs, then consensus-50 needs at least 9 to match), then return "nc"
full.150$cSC50 <- apply(full.150[, col.nam$c150], 1, function (x) ifelse(length(which(table(x) > c50_cutoff$full)) > 0, names(table(x))[which(table(x) > c50_cutoff$full)], "nc"))

# calculate the strict consensus values for the smaller datasets (i.e. slide and digital separated)
full.125$sSC50 <- apply(full.125[, col.nam$s125], 1, function (x) ifelse(length(which(table(x) > c50_cutoff$slide)) > 0, names(table(x))[which(table(x) > c50_cutoff$slide)], "nc"))
full.150$sSC50 <- apply(full.150[, col.nam$s150], 1, function (x) ifelse(length(which(table(x) > c50_cutoff$slide)) > 0, names(table(x))[which(table(x) > c50_cutoff$slide)], "nc"))
full.125$dSC50 <- apply(full.125[, col.nam$d125], 1, function (x) ifelse(length(which(table(x) > c50_cutoff$digital)) > 0, names(table(x))[which(table(x) > c50_cutoff$digital)], "nc"))
full.150$dSC50 <- apply(full.150[, col.nam$d150], 1, function (x) ifelse(length(which(table(x) > c50_cutoff$digital)) > 0, names(table(x))[which(table(x) > c50_cutoff$digital)], "nc"))
 
# 2b. Consensus ID --------------------------------------------------------
# for consensus minimum
# check what the maximum sum is for the consensus (so what is the highest count per specimen) ignoring na's
full.125$cMaxCon <- apply(full.125[, col.nam$c125], 1, function (x)  max(table(x[x != 'na'])))
full.150$cMaxCon <- apply(full.150[, col.nam$c150], 1, function (x) max(table(x[x != 'na'])))

# and split by smaller datasets
full.125$sMaxCon <- apply(full.125[, col.nam$s125], 1, function (x)  max(table(x[x != 'na'])))
full.150$sMaxCon <- apply(full.150[, col.nam$s150], 1, function (x) max(table(x[x != 'na'])))
full.125$dMaxCon <- apply(full.125[, col.nam$d125], 1, function (x) max(table(x[x != 'na'])))
full.150$dMaxCon <- apply(full.150[, col.nam$d150], 1, function (x) max(table(x[x != 'na'])))

# look at the summaries of these
table(full.125$cMaxCon) # so the minimum consensus is 19% (5/26)
table(full.150$cMaxCon)
table(full.125$sMaxCon) # for the slides, it is 18% (3/17)
table(full.150$sMaxCon)
table(full.125$dMaxCon) # for the digital analysis it would be 22% (2/9)
table(full.150$dMaxCon)

# calculate consensus ID
# I'm calculating this as the minimum consensus. So take most frequent name - if there are multiple take the first alphabetically (this is based on abbreviations)
# where the maximum is 'na', ignore it
# calculate for the full dataset, and slide / digital separated
full.125$cCID <- apply(full.125[, col.nam$c125], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][1])
full.150$cCID <- apply(full.150[, col.nam$c150], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][1])
full.125$sCID <- apply(full.125[, col.nam$s125], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][1])
full.150$sCID <- apply(full.150[, col.nam$s150], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][1])
full.125$dCID <- apply(full.125[, col.nam$d125], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][1])
full.150$dCID <- apply(full.150[, col.nam$d150], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][1])

# 2c. Consider the consensus agreement ------------------------------------------
# show how the data is split into one or multiple equally good sets
size125[, c("cCon1", "cCon2", "cAgreement")] <- NA
size150[, c("cCon1", "cCon2", "cAgreement")] <- NA

for (i in 1:2) {
  # add consensus values for 125 / 150
  size125[, grep("^cCon", names(size125))[i]] <- apply(full.125[, col.nam$c125], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][i])
  size150[, grep("^cCon", names(size150))[i]] <- apply(full.150[, col.nam$c150], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][i])
}
rm(i)

# add in the specimen agreement values
size125$cAgreement <- full.125$cMaxCon/length(col.nam$c125)*100
size150$cAgreement <- full.150$cMaxCon/length(col.nam$c150)*100

# what about for the slides / digital
size125[, c("sCon1", "sCon2", "sCon3", "sCon4", "sAgreement")] <- NA
size150[, c("sCon1", "sCon2", "sAgreement")] <- NA
size125[, c("dCon1", "dCon2", "dCon3", "dCon4", "dAgreement")] <- NA
size150[, c("dCon1", "dCon2", "dAgreement")] <- NA

for (i in 1:4) {
  # add consensus values for 125 / 150
  size125[, grep("^sCon", names(size125))[i]] <- apply(full.125[, col.nam$s125], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][i])
  size125[, grep("^dCon", names(size125))[i]] <- apply(full.125[, col.nam$d125], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][i])
  if (i < 3) {
    size150[, grep("^sCon", names(size150))[i]] <- apply(full.150[, col.nam$s150], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][i])
    size150[, grep("^dCon", names(size150))[i]] <- apply(full.150[, col.nam$d150], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][i])
  }
}
rm(i)

# add in the specimen agreement values
size125$sAgreement <- full.125$sMaxCon/length(col.nam$s125)*100
size150$sAgreement <- full.150$sMaxCon/length(col.nam$s150)*100
size125$dAgreement <- full.125$dMaxCon/length(col.nam$d125)*100
size150$dAgreement <- full.150$dMaxCon/length(col.nam$d150)*100


# there are more cases where the slide wins out than where the digital wins out. Sometimes including both helps break a tie, other times it picks a new consensus

# what fraction of specimens don't have a single consensus
sum(!is.na(size125$cCon2)); sum(!is.na(size125$sCon2)); sum(!is.na(size125$dCon2)) # 6 for combined vs 16/34 for slide / digital separately
sum(!is.na(size150$cCon2)); sum(!is.na(size150$sCon2)); sum(!is.na(size150$dCon2)) # 5 for combined vs 11/18 for slide / digital separately 

# have three possibilities
sum(!is.na(size125$sCon3)); sum(!is.na(size125$dCon3)) # 0 for combined vs 1/4 for slide / digital separately
# 150 have 0 for any of them

# have four possibilities
sum(!is.na(size125$sCon4)); sum(!is.na(size125$dCon4)) # 0 for combined vs 1/1 for slide / digital separately
# 150 have 0 for any of them

# 2d. Generate a species total dataframe ----------------------------------
# create a blank data frame
full.125sp <- data.frame(species = c(sp.abb$Abbreviation[1:2], sort(sp.abb$Abbreviation[3:nrow(sp.abb)])))
# ignore the MaxCon columns and the specimen one
full.125sp[, names(full.125)[!grepl("Specimen|MaxCon", names(full.125))]] <- NA
head(full.125sp)
# fill that dataframe
for (i in names(full.125sp)[2:ncol(full.125sp)]) {
  # table each column
  tmp <- table(full.125[,i])
  # add them in in the right order
  full.125sp[, i] <- tmp[match(full.125sp$species, names(tmp))]  
}

full.125sp[is.na(full.125sp)] <- 0
rm(i, tmp)

# repeat for 150 size fraction
# create a blank data frame
full.150sp <- data.frame(species = c(sp.abb$Abbreviation[1:2], sort(sp.abb$Abbreviation[3:nrow(sp.abb)])))
# ignore the MaxCon columns and the specimen one
full.150sp[, names(full.150)[!grepl("Specimen|MaxCon", names(full.150))]] <- NA
head(full.150sp)
# fill that dataframe
for (i in names(full.150sp)[2:ncol(full.150sp)]) {
  # table each column
  tmp <- table(full.150[,i])
  # add them in in the right order
  full.150sp[, i] <- tmp[match(full.150sp$species, names(tmp))]  
}

full.150sp[is.na(full.150sp)] <- 0
rm(i, tmp)

# 3. Agreement between workers (pairwise comparisons) ---------------------

# 3a. percentage agreement for each participant ----------------------------------
accuracy <- data.frame(PersonID = names(full.125)[col.nam$c125], Analysis = c(rep("Slide", 17), rep("Digital", 9)), stringsAsFactors = FALSE)

# based on the consensus for the combined
accuracy$cPtA125 <- apply(full.125[,col.nam$c125], 2, function(x) sum(x == full.125$cCID) / 300 * 100)
accuracy$cPtA150 <- apply(full.150[,col.nam$c150], 2, function(x) sum(x == full.150$cCID) / 300 * 100)

# then split by slide / digital
accuracy$dPtA150 <- accuracy$dPtA125 <- accuracy$sPtA150 <- accuracy$sPtA125 <- NA 
accuracy$sPtA125[accuracy$Analysis == "Slide"] <- apply(full.125[,col.nam$s125], 2, function(x) sum(x == full.125$sCID) / 300 * 100)
accuracy$sPtA150[accuracy$Analysis == "Slide"] <- apply(full.150[,col.nam$s150], 2, function(x) sum(x == full.150$sCID) / 300 * 100)
accuracy$dPtA125[accuracy$Analysis == "Digital"] <- apply(full.125[,col.nam$d125], 2, function(x) sum(x == full.125$dCID) / 300 * 100)
accuracy$dPtA150[accuracy$Analysis == "Digital"] <- apply(full.150[,col.nam$d150], 2, function(x) sum(x == full.150$dCID) / 300 * 100)

# looking at the summary of this data
tapply(accuracy$cPtA125, accuracy$Analysis, summary)
tapply(accuracy$cPtA150, accuracy$Analysis, summary)
summary(accuracy$sPtA125)
summary(accuracy$sPtA150)
summary(accuracy$dPtA125)
summary(accuracy$dPtA150)


# 3b. CID mean / sd of agreement -----------------------------------------------------
# calculate the mean percentage agreement for each of the four analyses, using initially the combined consensus
CID_mn <- list()
CID_mn$csPtA125 <- mean(accuracy$cPtA125[accuracy$Analysis == "Slide"])
CID_mn$csPtA150 <- mean(accuracy$cPtA150[accuracy$Analysis == "Slide"])
CID_mn$cdPtA125 <- mean(accuracy$cPtA125[accuracy$Analysis == "Digital"])
CID_mn$cdPtA150 <- mean(accuracy$cPtA150[accuracy$Analysis == "Digital"])
# and with the separate consensus'
CID_mn$ssPtA125 <- mean(accuracy$sPtA125, na.rm = TRUE)
CID_mn$ssPtA150 <- mean(accuracy$sPtA150, na.rm = TRUE)
CID_mn$ddPtA125 <- mean(accuracy$dPtA125, na.rm = TRUE)
CID_mn$ddPtA150 <- mean(accuracy$dPtA150, na.rm = TRUE)

# calculate the sd percentage agreement for each of the four analyses, using initially the combined consensus
CID_sd <- list()
CID_sd$csPtA125 <- sd(accuracy$cPtA125[accuracy$Analysis == "Slide"])
CID_sd$csPtA150 <- sd(accuracy$cPtA150[accuracy$Analysis == "Slide"])
CID_sd$cdPtA125 <- sd(accuracy$cPtA125[accuracy$Analysis == "Digital"])
CID_sd$cdPtA150 <- sd(accuracy$cPtA150[accuracy$Analysis == "Digital"])
# and with the separate consensus'
CID_sd$ssPtA125 <- sd(accuracy$sPtA125, na.rm = TRUE)
CID_sd$ssPtA150 <- sd(accuracy$sPtA150, na.rm = TRUE)
CID_sd$ddPtA125 <- sd(accuracy$dPtA125, na.rm = TRUE)
CID_sd$ddPtA150 <- sd(accuracy$dPtA150, na.rm = TRUE)

# mean consensus value based on doing a full pairwise comparison - these are the same as calculated above, so no need to do them.
sum(full.125$cCID == full.125[,col.nam$s125]) / (300*(length(col.nam$s125))) * 100
CID_mn$csPtA125

# 3c. Look at this agreement as plots ----------------------------
# generate a vector to order the data
ord.div <- c("1a", "1b", "2a", "2b", "A", 3:6, "F", 7:9, "G", 10:15, LETTERS[2:5], LETTERS[8:9])

# Compare digital / slide accuracy to the combined consensus (n.b. order is approx. experience)
png("ASFigures/CombCon_agreement_fullID.png", 800, 1000)
par(mfrow = c(2, 1))
with(accuracy, plot(cPtA125[match(ord.div, PersonID)], ylim = c(40, 90), type = "n", xaxt = "n", ylab = "Percentage agreement", xlab = "Participant", cex.lab = 1.5, las = 2, cex.axis = 1.1))
axis(1, at = 1:26, labels = accuracy$PersonID[match(ord.div, accuracy$PersonID)], cex.axis = 1.1)
# adding mean values
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA125, 2), lty = 1)
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA125 - CID_sd$csPtA125, 2), lty = 4)
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA125 + CID_sd$csPtA125, 2), lty = 4)
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA125, 2), lty = 1, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA125 - CID_sd$cdPtA125, 2), lty = 4, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA125 + CID_sd$cdPtA125, 2), lty = 4, col = "blue")
abline(v = 20.5, col = "grey 50")
# points
with(accuracy, points(1:26, cPtA125[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
with(accuracy, arrows(1, cPtA125[PersonID == "1a"], 2, cPtA125[PersonID == "1b"], length = 0.14))
with(accuracy, arrows(3, cPtA125[PersonID == "2a"], 4, cPtA125[PersonID == "2b"], length = 0.14))
text(26, 90, "125", cex = 1.5)

with(accuracy, plot(cPtA150[match(ord.div, PersonID)], ylim = c(40, 90), type = "n", xaxt = "n", ylab = "Percentage agreement", xlab = "Participant", cex.lab = 1.5, las = 2, cex.axis = 1.1))
axis(1, at = 1:26, labels = accuracy$PersonID[match(ord.div, accuracy$PersonID)], cex.axis = 1.1)
# adding mean values
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA150, 2), lty = 1)
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA150 - CID_sd$csPtA150, 2), lty = 4)
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA150 + CID_sd$csPtA150, 2), lty = 4)
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA150, 2), lty = 1, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA150 - CID_sd$cdPtA150, 2), lty = 4, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA150 + CID_sd$cdPtA150, 2), lty = 4, col = "blue")
abline(v = 20.5, col = "grey 50")
# points
with(accuracy, points(1:26, cPtA150[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
with(accuracy, arrows(1, cPtA150[PersonID == "1a"], 2, cPtA150[PersonID == "1b"], length = 0.14))
with(accuracy, arrows(3, cPtA150[PersonID == "2a"], 4, cPtA150[PersonID == "2b"], length = 0.14))
text(26, 90, "150", cex = 1.5)

par(mfrow = c(1,1))
dev.off()

# Compare digital / slide accuracy to the combined consensus (n.b. order is approx. experience)
png("ASFigures/SepCon_agreement_fullID.png", 800, 1000)
par(mfrow = c(2, 1))
with(accuracy, plot(sPtA125[match(ord.div, PersonID)], ylim = c(40, 90), type = "n", xaxt = "n", ylab = "Percentage agreement", xlab = "Participant", cex.lab = 1.5, las = 2, cex.axis = 1.1))
axis(1, at = 1:26, labels = accuracy$PersonID[match(ord.div, accuracy$PersonID)], cex.axis = 1.1)
# adding mean values
lines(x = c(0, 20.5), y = rep(CID_mn$ssPtA125, 2), lty = 1)
lines(x = c(0, 20.5), y = rep(CID_mn$ssPtA125 - CID_sd$ssPtA125, 2), lty = 4)
lines(x = c(0, 20.5), y = rep(CID_mn$ssPtA125 + CID_sd$ssPtA125, 2), lty = 4)
lines(x = c(20.5, 27), y = rep(CID_mn$ddPtA125, 2), lty = 1, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$ddPtA125 - CID_sd$ddPtA125, 2), lty = 4, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$ddPtA125 + CID_sd$ddPtA125, 2), lty = 4, col = "blue")
abline(v = 20.5, col = "grey 50")
# points
with(accuracy, points(1:26, sPtA125[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
with(accuracy, points(1:26, dPtA125[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
with(accuracy, arrows(1, sPtA125[PersonID == "1a"], 2, sPtA125[PersonID == "1b"], length = 0.14))
with(accuracy, arrows(3, sPtA125[PersonID == "2a"], 4, sPtA125[PersonID == "2b"], length = 0.14))
text(26, 90, "125", cex = 1.5)

with(accuracy, plot(sPtA150[match(ord.div, PersonID)], ylim = c(40, 90), type = "n", xaxt = "n", ylab = "Percentage agreement", xlab = "Participant", cex.lab = 1.5, las = 2, cex.axis = 1.1))
axis(1, at = 1:26, labels = accuracy$PersonID[match(ord.div, accuracy$PersonID)], cex.axis = 1.1)
# adding mean values
lines(x = c(0, 20.5), y = rep(CID_mn$ssPtA150, 2), lty = 1)
lines(x = c(0, 20.5), y = rep(CID_mn$ssPtA150 - CID_sd$ssPtA150, 2), lty = 4)
lines(x = c(0, 20.5), y = rep(CID_mn$ssPtA150 + CID_sd$ssPtA150, 2), lty = 4)
lines(x = c(20.5, 27), y = rep(CID_mn$ddPtA150, 2), lty = 1, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$ddPtA150 - CID_sd$ddPtA150, 2), lty = 4, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$ddPtA150 + CID_sd$ddPtA150, 2), lty = 4, col = "blue")
abline(v = 20.5, col = "grey 50")
# points
with(accuracy, points(1:26, sPtA150[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
with(accuracy, points(1:26, dPtA150[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
with(accuracy, arrows(1, sPtA150[PersonID == "1a"], 2, sPtA150[PersonID == "1b"], length = 0.14))
with(accuracy, arrows(3, sPtA150[PersonID == "2a"], 4, sPtA150[PersonID == "2b"], length = 0.14))
text(26, 90, "150", cex = 1.5)

par(mfrow = c(1,1))
dev.off()


# 4. How does agreement depend on whether workers routinely count specimens? --------
# add in routine to the accuracy info
accuracy$Routine <- people$Routine[match(accuracy$PersonID, people$SlideID)]
accuracy$Routine[accuracy$PersonID %in% c("1a", "2a")] <- people$Routine[1:2]
accuracy$Routine[accuracy$PersonID %in% c("1b", "2b")] <- people$Routine[1:2]
accuracy$Routine[is.na(accuracy$Routine)] <- people$Routine[match(accuracy$PersonID, people$DigitalID)][!is.na(match(accuracy$PersonID, people$DigitalID))]

# check the summaries
tapply(accuracy$cPtA125, paste(accuracy$Routine, accuracy$Analysis), summary)
tapply(accuracy$cPtA150, paste(accuracy$Routine, accuracy$Analysis), summary)
# and the sample sizes
table(accuracy$Routine, accuracy$Analysis) # only one person in the digital analysis did not do routine counting

png("ASFigures/Routine_agreement.png")
par(mfrow = c(2,2))
with(accuracy[accuracy$Analysis == "Slide",], boxplot(cPtA125 ~ Routine, main = "Slide 125"))
with(accuracy[accuracy$Analysis == "Slide",], boxplot(cPtA150 ~ Routine, main = "Slide 150"))
with(accuracy[accuracy$Analysis == "Digital",], boxplot(cPtA125 ~ Routine, main = "Digital 125"))
with(accuracy[accuracy$Analysis == "Digital",], boxplot(cPtA150 ~ Routine, main = "Digital 150"))
par(mfrow = c(1,1))
dev.off()

# 5. Sensitivity to alphabetical order -----------------------------------
# Does the way in which the consensus split is chosen affect the results?
# try taking the last choice rather than the first

# 5a. Recalculate the Consensus ID ----------------------------------------
# for the combined consensus
full.125$cCIDr <- apply(full.125[, col.nam$c125], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1])
full.150$cCIDr <- apply(full.150[, col.nam$c150], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1])
# and split by slide / digital
full.125$sCIDr <- apply(full.125[, col.nam$s125], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1])
full.150$sCIDr <- apply(full.150[, col.nam$s150], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1])
full.125$dCIDr <- apply(full.125[, col.nam$d125], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1])
full.150$dCIDr <- apply(full.150[, col.nam$d150], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1])

# 5b. Consider the effect on agreement -------------------------------------
# based on the consensus for the combined
accuracy$cPtA125r <- apply(full.125[,col.nam$c125], 2, function(x) sum(x == full.125$cCIDr) / 300 * 100)
accuracy$cPtA150r <- apply(full.150[,col.nam$c150], 2, function(x) sum(x == full.150$cCIDr) / 300 * 100)

# then split by slide / digital
accuracy$dPtA150r <- accuracy$dPtA125r <- accuracy$sPtA150r <- accuracy$sPtA125r <- NA 
accuracy$sPtA125r[accuracy$Analysis == "Slide"] <- apply(full.125[,col.nam$s125], 2, function(x) sum(x == full.125$sCIDr) / 300 * 100)
accuracy$sPtA150r[accuracy$Analysis == "Slide"] <- apply(full.150[,col.nam$s150], 2, function(x) sum(x == full.150$sCIDr) / 300 * 100)
accuracy$dPtA125r[accuracy$Analysis == "Digital"] <- apply(full.125[,col.nam$d125], 2, function(x) sum(x == full.125$dCIDr) / 300 * 100)
accuracy$dPtA150r[accuracy$Analysis == "Digital"] <- apply(full.150[,col.nam$d150], 2, function(x) sum(x == full.150$dCIDr) / 300 * 100)

# calculate the mean using the combined consensus
CID_mn$csPtA125r <- mean(accuracy$cPtA125r[accuracy$Analysis == "Slide"])
CID_mn$csPtA150r <- mean(accuracy$cPtA150r[accuracy$Analysis == "Slide"])
CID_mn$cdPtA125r <- mean(accuracy$cPtA125r[accuracy$Analysis == "Digital"])
CID_mn$cdPtA150r <- mean(accuracy$cPtA150r[accuracy$Analysis == "Digital"])
# and with the separate consensus'
CID_mn$ssPtA125r <- mean(accuracy$sPtA125r, na.rm = TRUE)
CID_mn$ssPtA150r <- mean(accuracy$sPtA150r, na.rm = TRUE)
CID_mn$ddPtA125r <- mean(accuracy$dPtA125r, na.rm = TRUE)
CID_mn$ddPtA150r <- mean(accuracy$dPtA150r, na.rm = TRUE)

# calculate the sd percentage agreement for each of the four analyses, using initially the combined consensus
CID_sd$csPtA125r <- sd(accuracy$cPtA125r[accuracy$Analysis == "Slide"])
CID_sd$csPtA150r <- sd(accuracy$cPtA150r[accuracy$Analysis == "Slide"])
CID_sd$cdPtA125r <- sd(accuracy$cPtA125r[accuracy$Analysis == "Digital"])
CID_sd$cdPtA150r <- sd(accuracy$cPtA150r[accuracy$Analysis == "Digital"])
# and with the separate consensus'
CID_sd$ssPtA125r <- sd(accuracy$sPtA125r, na.rm = TRUE)
CID_sd$ssPtA150r <- sd(accuracy$sPtA150r, na.rm = TRUE)
CID_sd$ddPtA125r <- sd(accuracy$dPtA125r, na.rm = TRUE)
CID_sd$ddPtA150r <- sd(accuracy$dPtA150r, na.rm = TRUE)

# 5c. Influence on the agreement plots --------------------------------------------------
png("ASFigures/CombCon_agreement_fullID_rev.png", 800, 1000)
par(mfrow = c(2, 1))
with(accuracy, plot(cPtA125r[match(ord.div, PersonID)], ylim = c(40, 90), type = "n", xaxt = "n", ylab = "Percentage agreement", xlab = "Participant", cex.lab = 1.5, las = 2, cex.axis = 1.1))
axis(1, at = 1:26, labels = accuracy$PersonID[match(ord.div, accuracy$PersonID)], cex.axis = 1.1)
# adding mean values
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA125r, 2), lty = 1)
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA125r - CID_sd$csPtA125r, 2), lty = 4)
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA125r + CID_sd$csPtA125r, 2), lty = 4)
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA125r, 2), lty = 1, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA125r - CID_sd$cdPtA125r, 2), lty = 4, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA125r + CID_sd$cdPtA125r, 2), lty = 4, col = "blue")
abline(v = 20.5, col = "grey 50")
# points
with(accuracy, points(1:26, cPtA125r[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
with(accuracy, arrows(1, cPtA125r[PersonID == "1a"], 2, cPtA125r[PersonID == "1b"], length = 0.14))
with(accuracy, arrows(3, cPtA125r[PersonID == "2a"], 4, cPtA125r[PersonID == "2b"], length = 0.14))
text(26, 90, "125", cex = 1.5)

with(accuracy, plot(cPtA150r[match(ord.div, PersonID)], ylim = c(40, 90), type = "n", xaxt = "n", ylab = "Percentage agreement", xlab = "Participant", cex.lab = 1.5, las = 2, cex.axis = 1.1))
axis(1, at = 1:26, labels = accuracy$PersonID[match(ord.div, accuracy$PersonID)], cex.axis = 1.1)
# adding mean values
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA150r, 2), lty = 1)
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA150r - CID_sd$csPtA150r, 2), lty = 4)
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA150r + CID_sd$csPtA150r, 2), lty = 4)
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA150r, 2), lty = 1, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA150r - CID_sd$cdPtA150r, 2), lty = 4, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA150r + CID_sd$cdPtA150r, 2), lty = 4, col = "blue")
abline(v = 20.5, col = "grey 50")
# points
with(accuracy, points(1:26, cPtA150r[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
with(accuracy, arrows(1, cPtA150r[PersonID == "1a"], 2, cPtA150r[PersonID == "1b"], length = 0.14))
with(accuracy, arrows(3, cPtA150r[PersonID == "2a"], 4, cPtA150r[PersonID == "2b"], length = 0.14))
text(26, 90, "150", cex = 1.5)

par(mfrow = c(1,1))
dev.off()


# 6. Confusion matrix ----------------------------------------------------

# 6a. Confusion matrices for the main results -----------------------------
# initially with slide 125
# this requires the data to be in long format
long <- list()
long$s125 <- reshape(full.125[, -col.nam$d125], varying = list(names(full.125)[col.nam$s125]), direction = "long", times = names(full.125)[col.nam$s125], timevar = "Person")
rownames(long$s125) <- 1:nrow(long$s125)
long$s125 <- long$s125[, (names(long$s125) != "id")]
names(long$s125)[names(long$s125) == "1a"] <- "origID"
head(long$s125)
tail(long$s125)

png("ASFigures/CombCon_conf_slide125.png", 1000, 700)
conf_mat(long$s125, "origID", "cCID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Slide 125", 1, cex = 2, adj = -0.8,  line = -5)
dev.off()

# same for the other datasets
# so slide 150
long$s150 <- reshape(full.150[, -col.nam$d150], varying = list(names(full.150)[col.nam$s150]), direction = "long", times = names(full.150)[col.nam$s150], timevar = "Person")
rownames(long$s150) <- 1:nrow(long$s150)
long$s150 <- long$s150[, (names(long$s150) != "id")]
names(long$s150)[names(long$s150) == "1a"] <- "origID"
head(long$s150)
tail(long$s150)

png("ASFigures/CombCon_conf_slide150.png", 1000, 610)
conf_mat(long$s150, "origID", "cCID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Slide 150", 1, cex = 2, adj = -.8,  line = -5)
dev.off()

# digital 125
long$d125 <- reshape(full.125[, -col.nam$s125], varying = list(names(full.125)[col.nam$d125]), direction = "long", times = names(full.125)[col.nam$d125], timevar = "Person")
rownames(long$d125) <- 1:nrow(long$d125)
long$d125 <- long$d125[, (names(long$d125) != "id")]
names(long$d125)[names(long$d125) == "A"] <- "origID"
head(long$d125)
tail(long$d125)

png("ASFigures/CombCon_conf_digital125.png", 1000, 700)
conf_mat(long$d125, "origID", "cCID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Digital 125", 1, cex = 2, adj = -.8,  line = -5)
dev.off()

# digital 150
long$d150 <- reshape(full.150, varying = list(names(full.150)[col.nam$d150]), direction = "long", times = names(full.150)[col.nam$d150], timevar = "Person")
rownames(long$d150) <- 1:nrow(long$d150)
long$d150 <- long$d150[, (names(long$d150) != "id")]
names(long$d150)[names(long$d150) == "A"] <- "origID"
head(long$d150)
tail(long$d150)

png("ASFigures/CombCon_conf_digital150.png", 1000, 610)
conf_mat(long$d150, "origID", "cCID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Digital 150", 1, cex = 2, adj = -.8,  line = -5)
dev.off()

# What about for the separate consensus values
png("ASFigures/SepCon_conf_slide125.png", 1000, 700)
conf_mat(long$s125, "origID", "sCID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Sep Slide 125", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

png("ASFigures/SepCon_conf_slide150.png", 1000, 610)
conf_mat(long$s150, "origID", "sCID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Sep Slide 150", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

png("ASFigures/SepCon_conf_digital125.png", 1000, 700)
conf_mat(long$d125, "origID", "dCID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Sep Digital 125", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

png("ASFigures/SepCon_conf_digital150.png", 1000, 610)
conf_mat(long$d150, "origID", "dCID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Sep Digital 150", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

# 6b. Influence of reversed order on the confusion matrices ---------------------------------
# influence of the alphabetical ordering
png("ASFigures/CombCon_conf_slide125r.png", 1000, 700)
conf_mat(long$s125, "origID", "cCIDr", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Slide 125", 1, cex = 2, adj = -0.8,  line = -5)
dev.off()

png("ASFigures/CombCon_conf_slide150r.png", 1000, 610)
conf_mat(long$s150, "origID", "cCIDr", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Slide 150", 1, cex = 2, adj = -.8,  line = -5)
dev.off()

png("ASFigures/CombCon_conf_digital125r.png", 1000, 700)
conf_mat(long$d125, "origID", "cCIDr", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Digital 125", 1, cex = 2, adj = -.8,  line = -5)
dev.off()

png("ASFigures/CombCon_conf_digital150r.png", 1000, 610)
conf_mat(long$d150, "origID", "cCIDr", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Digital 150", 1, cex = 2, adj = -.8,  line = -5)
dev.off()

# and for the separate consensus values
png("ASFigures/SepCon_conf_slide125r.png", 1000, 700)
conf_mat(long$s125, "origID", "sCIDr", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Sep Slide 125", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

png("ASFigures/SepCon_conf_slide150r.png", 1000, 610)
conf_mat(long$s150, "origID", "sCIDr", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Sep Slide 150", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

png("ASFigures/SepCon_conf_digital125r.png", 1000, 700)
conf_mat(long$d125, "origID", "dCIDr", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Sep Digital 125", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

png("ASFigures/SepCon_conf_digital150r.png", 1000, 610)
conf_mat(long$d150, "origID", "dCIDr", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Sep Digital 150", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()
# # influence of the alphabetical ordering
# png("ASFigures/confusion_slide125r.png", 1000, 700)
# conf_mat(long$s125, "origID", "CIDR", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
# mtext("Slide 125", 1, cex = 2, adj = -0.8,  line = -5)
# dev.off()
# 
# png("ASFigures/confusion_slide150r.png", 1000, 600)
# conf_mat(long$s150, "origID", "CIDR", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
# mtext("Slide 150", 1, cex = 2, adj = -.8,  line = -5)
# dev.off()
# 
# png("ASFigures/confusion_digital125r.png", 1000, 700)
# conf_mat(long$d125, "origID", "CIDR", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
# mtext("Digital 125", 1, cex = 2, adj = -.8,  line = -5)
# dev.off()
# 
# png("ASFigures/confusion_digital150r.png", 1000, 600)
# conf_mat(long$d150, "origID", "CIDR", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
# mtext("Digital 150", 1, cex = 2, adj = -.8,  line = -5)
# dev.off()

# slide125$CIDR <- apply(slide125[, col.nam$s125], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1]) 
# slide150$CIDR <- apply(slide150[, col.nam$s150], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1]) 
# digital125$CIDR <- apply(digital125[, col.nam$d125], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1]) 
# digital150$CIDR <- apply(digital150[, col.nam$d150], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1]) 
# 
# accuracySlide$IF_revPtAc125 <- apply(slide125[,col.nam$s125], 2, function(x) sum(x == slide125$CIDR) / 300 * 100)[accuracySlide$PersonID]
# accuracySlide$IF_revPtAc150 <- apply(slide150[,col.nam$s150], 2, function(x) sum(x == slide150$CIDR) / 300 * 100)[accuracySlide$PersonID]
# accuracyDigital$IF_revPtAc125 <- apply(digital125[,col.nam$d125], 2, function(x) sum(x == digital125$CIDR) / 300 * 100)
# accuracyDigital$IF_revPtAc150 <- apply(digital150[,col.nam$d150], 2, function(x) sum(x == digital150$CIDR) / 300 * 100)
# 
# CID_mn$IF_revPtAc125s <- mean(accuracySlide$IF_revPtAc125)
# CID_mn$IF_revPtAc150s <- mean(accuracySlide$IF_revPtAc150)
# CID_mn$IF_revPtAc125d <- mean(accuracyDigital$IF_revPtAc125)
# CID_mn$IF_revPtAc150d <- mean(accuracyDigital$IF_revPtAc150)
# 
# # My value
# CID_sd$IF_revPtAc125s <- sd(accuracySlide$IF_revPtAc125)
# CID_sd$IF_revPtAc150s <- sd(accuracySlide$IF_revPtAc150)
# CID_sd$IF_revPtAc125d <- sd(accuracyDigital$IF_revPtAc125)
# CID_sd$IF_revPtAc150d <- sd(accuracyDigital$IF_revPtAc150)
# 
# # plot it up
# ord.div <- c("1a", "1b", "2a", "2b", "A", 3:6, "F", 7:9, "G", 10:15, LETTERS[2:5], LETTERS[8:9])
# accuracy <- rbind(accuracySlide, accuracyDigital)
# accuracy$Analysis <- c(rep("Slide", 17), rep("Digital", 9))
# 
# png("ASFigures/Fig3_Consensus_agreement_rev.png", 800, 1000)
# par(mfrow = c(2, 1))
# with(accuracy, plot(IF_revPtAc125[match(ord.div, PersonID)], ylim = c(43, 90), type = "n", xaxt = "n", ylab = "Percentage accuracy", xlab = "Participant", cex.lab = 1.5, las = 2, cex.axis = 1.1))
# axis(1, at = 1:26, labels = accuracy$PersonID[match(ord.div, accuracy$PersonID)], cex.axis = 1.1)
# abline(h = CID_mn$IF_revPtAc125s, lty = 1)
# abline(h = c(CID_mn$IF_revPtAc125s - CID_sd$IF_revPtAc125s, CID_mn$IF_revPtAc125s + CID_sd$IF_revPtAc125s), lty = 4)
# abline(h = CID_mn$IF_revPtAc125d, lty = 1, col = "blue")
# abline(h = c(CID_mn$IF_revPtAc125d - CID_sd$IF_revPtAc125d, CID_mn$IF_revPtAc125d + CID_sd$IF_revPtAc125d), lty = 4, col = "blue")
# with(accuracy, points(1:26, IF_revPtAc125[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
# with(accuracy, arrows(1, IF_revPtAc125[PersonID == "1a"], 2, IF_revPtAc125[PersonID == "1b"], length = 0.14))
# with(accuracy, arrows(3, IF_revPtAc125[PersonID == "2a"], 4, IF_revPtAc125[PersonID == "2b"], length = 0.14))
# text(26, 90, "125", cex = 1.5)
# 
# with(accuracy, plot(IF_revPtAc150[match(ord.div, PersonID)], ylim = c(43, 90), type = "n", xaxt = "n", ylab = "Percentage accuracy", xlab = "Participant", cex.lab = 1.5, las = 2, cex.axis = 1.1))
# axis(1, at = 1:26, labels = accuracy$PersonID[match(ord.div, accuracy$PersonID)], cex.axis = 1.1)
# abline(h = CID_mn$IF_revPtAc150s, lty = 1)
# abline(h = c(CID_mn$IF_revPtAc150s - CID_sd$IF_revPtAc150s, CID_mn$IF_revPtAc150s + CID_sd$IF_revPtAc150s), lty = 4)
# abline(h = CID_mn$IF_revPtAc150d, lty = 1, col = "blue")
# abline(h = c(CID_mn$IF_revPtAc150d - CID_sd$IF_revPtAc150d, CID_mn$IF_revPtAc150d + CID_sd$IF_revPtAc150d), lty = 4, col = "blue")
# with(accuracy, points(1:26, IF_revPtAc150[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
# with(accuracy, arrows(1, IF_revPtAc150[PersonID == "1a"], 2, IF_revPtAc150[PersonID == "1b"], length = 0.14))
# with(accuracy, arrows(3, IF_revPtAc150[PersonID == "2a"], 4, IF_revPtAc150[PersonID == "2b"], length = 0.14))
# text(26, 90, "150", cex = 1.5)
# 
# par(mfrow = c(1,1))
# dev.off()
# 
# 7. Species level accuracy ------------------------------------------------
# What if instead of looking at the overall accuracy, we look at a species level version.
# only do this for the combined consensus
# for 125
acc.sp.125 <- data.frame(Species = c(sp.abb$Abbreviation[1:2], sort(sp.abb$Abbreviation[3:nrow(sp.abb)])))
acc.sp.125[, names(full.150)[col.nam$c125]] <- NA

# how often did they get that species right?
for (i in 1:nrow(acc.sp.125)) {
  acc.sp.125[i, col.nam$c125] <- apply(full.125[,col.nam$c125], 2, function (x) sum(x[full.125$cCID == acc.sp.125$Species[i]] == full.125$cCID[full.125$cCID == acc.sp.125$Species[i]])/sum(full.125$cCID == acc.sp.125$Species[i]))
}

# look at some summary stats - split by person
apply(acc.sp.125[,col.nam$c125], 2, summary)
# what fraction of species did they always get right
apply(acc.sp.125[,col.nam$c125], 2, function(x) sum(x == 1, na.rm = TRUE) / sum(!is.na(x)))

# how does the mean vary between slide and digital
mean(apply(acc.sp.125[,col.nam$s125], 2, mean, na.rm = TRUE)) * 100
mean(apply(acc.sp.125[,col.nam$d125], 2, mean, na.rm = TRUE)) * 100
# how does this compare to the accuracy at the specimen level?
mean(accuracy$cPtA125)

# for 150
acc.sp.150 <- data.frame(Species = c(sp.abb$Abbreviation[1:2], sort(sp.abb$Abbreviation[3:nrow(sp.abb)])))
acc.sp.150[, names(full.150)[col.nam$c150]] <- NA

# how often did they get that species right?
for (i in 1:nrow(acc.sp.150)) {
  acc.sp.150[i, col.nam$c150] <- apply(full.150[,col.nam$c150], 2, function (x) sum(x[full.150$cCID == acc.sp.150$Species[i]] == full.150$cCID[full.150$cCID == acc.sp.150$Species[i]])/sum(full.150$cCID == acc.sp.150$Species[i]))
}

# look at some summary stats - split by person
apply(acc.sp.150[,col.nam$c150], 2, summary)
# what fraction of species did they always get right
apply(acc.sp.150[,col.nam$c150], 2, function(x) sum(x == 1, na.rm = TRUE) / sum(!is.na(x)))

# how does the mean vary between slide and digital
mean(apply(acc.sp.150[,col.nam$s150], 2, mean, na.rm = TRUE)) * 100
mean(apply(acc.sp.150[,col.nam$d150], 2, mean, na.rm = TRUE)) * 100
# how does this compare to the accuracy at the specimen level?
mean(accuracy$cPtA150)

# answer appears to be that it doesn't make things any better

# Are participants consistent on species?
# par(ask = TRUE)
for (i in 2:ncol(acc.sp.125)) {
  plot(acc.sp.125[, i], acc.sp.150[, i], pch = 16, main = names(acc.sp.125)[i])
  abline(lm(acc.sp.150[, i] ~ acc.sp.125[, i]))
}
# par(ask = FALSE)

# # 3l. Removing <150 from <125 for comparison ------------------------------
# size125[which(size125$Length < 150),]
# summary(size125$slideAgreement[which(size125$Length < 150)])
# summary(size125$slideAgreement[which(size125$Length > 150)])
# summary(size150$slideAgreement)
# 
# summary(size125$digitalAgreement[which(size125$Length < 150)])
# summary(size125$digitalAgreement[which(size125$Length > 150)])
# summary(size150$digitalAgreement)
# 
# tmp <- data.frame(size = seq(130, 710, by = 10))
# tmp$c125[tmp$size %in% names(table(round(size125$Length, -1)))] <- table(round(size125$Length, -1))
# tmp$c150[tmp$size %in% names(table(round(size150$Length, -1)))] <- table(round(size150$Length, -1))
# 
# png("ASFigures/SizeComparison.png", 800, 500)
# barplot(t(tmp[, 2:3]), beside = TRUE, legend.text = c("125", "150"), names.arg = tmp$size, las = 2)
# dev.off()
# 
# # 4. NMDS -----------------------------------------------------------------
# 
# # 4a. Recreating Nadia's NMDS ---------------------------------------------
# # I think this was originally run on the species counts (i.e. working at the community level)
# # I attempted to do this and failed to produce the same results
# 
# # Figure 2
# # for size fraction 125 (slide and digital are combined, with the consensus values)
# # the data needs to be the other way round, so transpose it
# sp.125 <- merge(slide125sp, digital125sp, by = "species")
# 
# trsp <- list()
# trsp$NA125sp <- data.frame(t(sp.125[, !(names(sp.125) %in% c("species", "MaxCon.x", "MaxCon.y", "SC50.x", "CID.x", "SC50.y", "CID.y"))]))
# names(trsp$NA125sp) <- sp.125$species
# rownames(trsp$NA125sp)[nchar(rownames(trsp$NA125sp)) >= 5] <- c("Sc50", "Sc20", "Dc50", "Dc20")
# # re-running this creates a certain amount of movement
# nmds <- list()
# nmds$NA125sp <- metaMDS(trsp$NA125sp)
# 
# # create a dataframe for the colours
# mds.col <- data.frame(person = rownames(trsp$NA125sp))
# # add in the school
# mds.col$school <- people$School[match(mds.col$person, people$SlideID)]
# mds.col$school[mds.col$person %in% people$DigitalID] <- people$School[order(people$DigitalID)[1:sum(!is.na(people$DigitalID))]]
# mds.col$school[grep("1[a-z]", mds.col$person)] <- people$School[people$SlideID == "1"][1]
# mds.col$school[grep("2[a-z]", mds.col$person)] <- people$School[people$SlideID == "2"][1]
# mds.col$school[is.na(mds.col$school)] <- as.character(mds.col$person[is.na(mds.col$school)])
# # create a colour for the school
# mds.col$sch.col <- NA
# mds.col$sch.col[grep("^[1-9]", mds.col$school)] <- mds.col$school[grep("^[1-9]", mds.col$school)]
# mds.col$sch.col <- gsub("_.*", "", mds.col$sch.col)
# mds.col$sch.col <- as.numeric(mds.col$sch.col)
# 
# # add in the paired analyses
# mds.col$pair <- NA
# for (i in which(!is.na(people$SlideID) & !is.na(people$DigitalID))) {
#   mds.col$pair[mds.col$person %in% people$DigitalID[i]] <- i
#   mds.col$pair[grep(paste("^",people$SlideID[i], sep = ""), mds.col$person)] <- i
# }
# mds.col$pair[mds.col$person == "2a"] <- NA
# rm(i)
# 
# # plot the NMDS
# plot(nmds$NA125sp, type = "n", display = "sites")
# points(nmds$NA125sp, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
# text(nmds$NA125sp, labels = rownames(trsp$NA125sp))
# legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16)
# 
# # for 150
# sp.150 <- merge(slide150sp, digital150sp, by = "species")
# trsp$NA150sp <- data.frame(t(sp.150[, !(names(sp.150) %in% c("species", "MaxCon.x", "MaxCon.y", "SC50.x", "CID.x", "SC50.y", "CID.y"))]))
# names(trsp$NA150sp) <- sp.150$species
# rownames(trsp$NA150sp)[nchar(rownames(trsp$NA150sp)) >= 5] <- c("Sc50", "Sc20", "Dc50", "Dc20")
# nmds$NA150sp <- metaMDS(trsp$NA150sp)
# plot(nmds$NA150sp, type = "n", display = "sites")
# points(nmds$NA150sp, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
# text(nmds$NA150sp, labels = rownames(trsp$NA150sp))
# legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16)
# # these are currently very different from Nadia's data
# 
# # 4b. Using raw data rather than species counts ---------------------------
# # again this is run initially for Nadia's consensus values
# # I think instead it is more sensible to use the raw data
# # try using the original data instead (n.b. I can't seem to do this in PAST)
# 
# # for 125
# full.125 <- merge(slide125, digital125, by = "Specimen")
# trsp$NA125f <- data.frame(t(full.125[, !(names(full.125) %in% c("Specimen", "MaxCon.x", "MaxCon.y", "SC50.x", "CID.x", "CIDR.x", "SC50.y", "CID.y", "CIDR.y"))]))
# rownames(trsp$NA125f)[nchar(rownames(trsp$NA125f)) >= 5] <- c("Sc50", "Sc20", "Dc50", "Dc20")
# nmds$NA125f <- metaMDS(daisy(trsp$NA125f))
# 
# # consider the stress of the NMDS
# stress <- list()
# stress$NA125f <- rep(NA, 10)
# for (i in 1:10) {
#   stress$NA125f[i] <- metaMDS(daisy(trsp$NA125f), k = i)$stress
# }
# plot(stress$NA125f, type = "b")  
# rm(i)
# stressplot(nmds$NA125f)
# # looks like between 2 and 3 dimensions would be reasonable
# 
# plot(nmds$NA125f, type = "n", display = "sites")
# points(nmds$NA125f, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
# text(nmds$NA125f, labels = rownames(trsp$NA125f))
# legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16)
# 
# 
# # and for 150
# full.150 <- merge(slide150, digital150, by = "Specimen")
# trsp$NA150f <- data.frame(t(full.150[, !(names(full.150) %in% c("Specimen", "MaxCon.x", "MaxCon.y", "SC50.x", "CID.x", "CIDR.x", "SC50.y", "CID.y", "CIDR.y"))]))
# rownames(trsp$NA150f)[nchar(rownames(trsp$NA150f)) >= 5 & !grepl("MinR", rownames(trsp$NA125f))] <- c("Sc50", "Sc20", "Dc50", "Dc20")
# nmds$NA150f <- metaMDS(daisy(trsp$NA150f))
# 
# # consider the stress of the NMDS
# stress$NA150f <- rep(NA, 10)
# for (i in 1:10) {
#   stress$NA150f[i] <- metaMDS(daisy(trsp$NA150f), k = i)$stress
# }
# plot(stress$NA150f, type = "b")  
# rm(i)
# stressplot(nmds$NA150f)
# # again between 2 and 3 dimensions is probably reasonable
# 
# # the full plot
# plot(nmds$NA150f, type = "n", display = "sites", cex = 1)
# points(nmds$NA150f, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
# text(nmds$NA150f, labels = rownames(trsp$NA150f))
# legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16)
# 
# # Given the influence of the outliers, a more informative relationship between these points can be obtained by running the analysis excluding the outliers. 
# nmds$NA150z <- metaMDS(daisy(trsp$NA150f[!(rownames(trsp$NA150f) %in% c("3", "C", "E", "G")),]))
# 
# # consider the stress of the NMDS
# stress$NA150z <- rep(NA, 10)
# for (i in 1:10) {
#   stress$NA150z[i] <- metaMDS(daisy(trsp$NA150f[!(rownames(trsp$NA150f) %in% c("3", "C", "E", "G")),]), k = i)$stress
# }
# plot(stress$NA150z, type = "b")
# rm(i)
# stressplot(nmds$NA150z)
# # again between 2 and 3 dimensions is probably reasonable
# 
# # the full plot
# plot(nmds$NA150z, type = "n", display = "sites", cex = 1)
# points(nmds$NA150z, pch = 21, cex = 4, col = mds.col$pair[!(mds.col$person %in% c("3", "C", "E", "G"))], bg = brewer.pal(5, "Set2")[mds.col$sch.col[!(mds.col$person %in% c("3", "C", "E", "G"))]], lwd = 2)
# text(nmds$NA150z, labels = rownames(trsp$NA150f[!(rownames(trsp$NA150f) %in% c("3", "C", "E", "G")),]))
# legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16)
# 
# # 4c. plotting using my consensus estimates -------------------------------
# # for 125
# trsp$IF125f <- data.frame(t(full.125[, !(names(full.125) %in% c("Specimen", "consensus50.x", "consensus20.x", "consensus50.y", "consensus20.y", "MaxCon.x", "CIDR.x", "MaxCon.y", "CIDR.y"))]))
# rownames(trsp$IF125f)[nchar(rownames(trsp$IF125f)) >= 5] <- c("SsC", "SCID", "DsC", "DCID")
# 
# # consider the stress of the NMDS
# stress$IF125f <- rep(NA, 10)
# for (i in 1:10) {
#   stress$IF125f[i] <- metaMDS(daisy(trsp$IF125f), k = i)$stress
# }
# plot(stress$IF125f, type = "b")  
# rm(i)
# # looks like between 2 and 3 dimensions would be reasonable
# 
# # check for variation
# tmp <- metaMDS(daisy(trsp$IF125f))
# tmp$stress
# plot(tmp, display = "sites", type = "t", cex = 1.5)
# # some variation, so optimise
# 
# stress$IFop125f <- 1
# # find the nmds plot with the lowest stress out of 20000 runs
# # for (i in 1:1000) {
# #   tmp <- metaMDS(daisy(trsp$IF125f))
# #   if (stress$IFop125f > tmp$stress) {
# #     nmds$IF125f <- tmp
# #     stress$IFop125f <- tmp$stress
# #   }
# # }
# # rm(i, tmp)
# # save(stress, nmds, file = "Outputs/NMDS_125.RData")
# load("Outputs/NMDS_125.RData")
# 
# stressplot(nmds$IF125f)
# 
# png("ASFigures/IF_NMDS_125.png", 600, 600)
# plot(nmds$IF125f, type = "n", display = "sites", cex = 1, xlab = "Axis 1", ylab = "Axis 2", las = 1, cex.lab = 1.5, cex.axis = 1.2, main = "NMDS 125")
# points(nmds$IF125f, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
# text(nmds$IF125f$points[nchar(rownames(nmds$IF125f$points)) <3, ], labels = rownames(nmds$IF125f$points)[nchar(rownames(nmds$IF125f$points)) <3])
# points(nmds$IF125f$points[nchar(rownames(nmds$IF125f$points)) >2, ], pch = "+")
# text(sweep(data.matrix(nmds$IF125f$points[nchar(rownames(nmds$IF125f$points)) >2 & rownames(nmds$IF125f$points) != "SCID", ]), 2, c(-0.03, 0)), labels = rownames(nmds$IF125f$points)[nchar(rownames(nmds$IF125f$points)) >2 & rownames(nmds$IF125f$points) != "SCID"])
# text(nmds$IF125f$points[rownames(nmds$IF125f$points) == "SCID", 1]+0.04, nmds$IF125f$points[rownames(nmds$IF125f$points) == "SCID", 2]+0.01, labels = "SCID")
# legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16, cex = 1.3)
# dev.off()
# 
# # for 150
# full.150 <- merge(slide150, digital150, by = "Specimen")
# trsp$IF150f <- data.frame(t(full.150[, !(names(full.150) %in% c("Specimen", "consensus50.x", "consensus20.x", "consensus50.y", "consensus20.y", "MaxCon.x", "CIDR.x", "MaxCon.y", "CIDR.y"))]))
# rownames(trsp$IF150f)[nchar(rownames(trsp$IF150f)) >= 5] <- c("SsC", "SCID", "DsC", "DCID")
# 
# # consider the stress of the NMDS
# stress$IF150f <- rep(NA, 10)
# for (i in 1:10) {
#   stress$IF150f[i] <- metaMDS(daisy(trsp$IF150f), k = i)$stress
# }
# plot(stress$IF150f, type = "b")
# rm(i)
# # looks like between 2 and 3 dimensions would be reasonable, although as with Nadia's data, the breakpoint is less obvious for 150 than it is for 125. 
# 
# # check for variation
# tmp <- metaMDS(daisy(trsp$IF150f))
# tmp$stress
# plot(tmp, display = "sites", type = "t", cex = 1.5)
# # some variation, so optimise
# 
# stress$IFop150f <- 1
# # find the nmds plot with the lowest stress out of 1000 runs
# # for (i in 1:1000) {
# #   tmp <- metaMDS(daisy(trsp$IF150f))
# #   if (stress$IFop150f > tmp$stress) {
# #     nmds$IF150f <- tmp
# #     stress$IFop150f <- tmp$stress
# #   }
# # }
# # rm(i, tmp)
# # save(stress, nmds, file = "Outputs/NMDS_150f.RData")
# load("Outputs/NMDS_150f.RData")
# 
# 
# 
# stressplot(nmds$IF150f)
# 
# # the full plot
# png("ASFigures/IF_NMDS_150.png", 600, 600)
# plot(nmds$IF150f, type = "n", display = "sites", cex = 1, xlab = "Axis 1", ylab = "Axis 2", las = 1, cex.lab = 1.5, cex.axis = 1.2, main = "NMDS 150")
# points(nmds$IF150f, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
# text(nmds$IF150f$points[nchar(rownames(nmds$IF150f$points)) <3, ], labels = rownames(nmds$IF150f$points)[nchar(rownames(nmds$IF150f$points)) <3])
# points(nmds$IF150f$points[nchar(rownames(nmds$IF150f$points)) >2, ], pch = "+")
# text(sweep(data.matrix(nmds$IF150f$points[nchar(rownames(nmds$IF150f$points)) >2 & rownames(nmds$IF150f$points) != "DsC", ]), 2, c(-0.025, 0)), labels = rownames(nmds$IF150f$points)[nchar(rownames(nmds$IF150f$points)) >2 & rownames(nmds$IF150f$points) != "DsC"])
# text(nmds$IF150f$points[rownames(nmds$IF150f$points) == "DsC", 1], nmds$IF150f$points[rownames(nmds$IF150f$points) == "DsC", 2]-0.02, labels = "DsC")
# legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16, cex = 1.3)
# dev.off()
# 
# # focussing on the main section. As noted above, it is better to run this as a new analysis rather than just zoom in, as the influence of the outliers means that the stability of the central points hasn't been tested
# nmds$IF150z <- metaMDS(daisy(trsp$IF150f[!(rownames(trsp$IF150f) %in% c("3", "C", "E", "G")),]))
# 
# # consider the stress of the NMDS
# stress$IF150z <- rep(NA, 10)
# for (i in 1:10) {
#   stress$IF150z[i] <- metaMDS(daisy(trsp$IF150f[!(rownames(trsp$IF150f) %in% c("3", "C", "E", "G")),]), k = i)$stress
# }
# plot(stress$IF150z, type = "b")
# rm(i)
# stressplot(nmds$IF150z)
# # again between 2 and 3 dimensions is probably reasonable
# 
# # check for variation
# tmp <- metaMDS(daisy(trsp$IF150f[!(rownames(trsp$IF150f) %in% c("3", "C", "E", "G")),]))
# tmp$stress
# plot(tmp, display = "sites", type = "t", cex = 1.5)
# # some variation, so optimise
# 
# stress$IFop150z <- 1
# # find the nmds plot with the lowest stress out of 1000 runs
# # for (i in 1:1000) {
# #   tmp <- metaMDS(daisy(trsp$IF150f[!(rownames(trsp$IF150f) %in% c("3", "C", "E", "G")),]))
# #   if (stress$IFop150z > tmp$stress) {
# #     nmds$IF150z <- tmp
# #     stress$IFop150z <- tmp$stress
# #   }
# # }
# # rm(i, tmp)
# # save(stress, nmds, file = "Outputs/NMDS_150z.RData")
# load("Outputs/NMDS_150z.RData")
# 
# stress$IFop150z
# stressplot(nmds$IF150z)
# 
# # the zoomed plot
# png("ASFigures/IF_NMDS_150_zoom.png", 600, 600)
# plot(nmds$IF150z, type = "n", display = "sites", cex = 1, xlab = "Axis 1", ylab = "Axis 2", las = 1, cex.lab = 1.5, cex.axis = 1.2, main = "NMDS 150 zoomed")
# points(nmds$IF150z, pch = 21, cex = 4, col = mds.col$pair[!(mds.col$person %in% c("3", "C", "E", "G"))], bg = brewer.pal(5, "Set2")[mds.col$sch.col[!(mds.col$person %in% c("3", "C", "E", "G"))]], lwd = 2)
# text(nmds$IF150z$points[nchar(rownames(nmds$IF150z$points)) <3, ], labels = rownames(nmds$IF150z$points)[nchar(rownames(nmds$IF150z$points)) <3])
# points(nmds$IF150z$points[nchar(rownames(nmds$IF150z$points)) >2, ], pch = "+")
# text(sweep(data.matrix(nmds$IF150z$points[nchar(rownames(nmds$IF150z$points)) >2 & rownames(nmds$IF150z$points) != "SsC", ]), 2, c(-0.012, -0.008)), labels = rownames(nmds$IF150z$points)[nchar(rownames(nmds$IF150z$points)) >2 & rownames(nmds$IF150z$points) != "SsC"])
# text(nmds$IF150z$points[rownames(nmds$IF150z$points) == "SsC", 1]-0.012, nmds$IF150z$points[rownames(nmds$IF150z$points) == "SsC", 2]-0.008, labels = "SsC")
# legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16, cex = 1.3)
# dev.off()
# 
# # output the scree plots
# png("ASFigures/Scree plots.png")
# par(mfrow = c(2, 2))
# plot(stress$IF125f, type = "b", bty = "l", las = 1, xlab = "Dimensions", ylab = "Stress", main = "NMDS 125")
# plot(stress$IF150f, type = "b", bty = "l", las = 1, xlab = "Dimensions", ylab = "Stress", main = "NMDS 150")
# plot(stress$IF150z, type = "b", bty = "l", las = 1, xlab = "Dimensions", ylab = "Stress", main = "NMDS 150 zoomed")
# par(mfrow = c(1,1))
# dev.off()
# 
# # 4b. Dendrogram ----------------------------------------------------------
# # an alternative way of plotting this is as a dendrogram
# # based on the community
# plot(hclust(daisy(data.frame(trsp$IF150f))))
# 
# # based on the original IDs
# plot(hclust(daisy(data.frame(t(merge(slide125, digital125, by = "Specimen")[, !(names(sp.125) %in% c("species", "MaxCon.x", "MaxCon.y", "SC50.x", "CID.x", "SC50.y", "CID.y"))])))))
# 
# # but I don't think this is as helpful
# 
# # 5. Repeated analysis by workers -----------------------------------------
# # ex figure 6
# # again, do this as a confusion matrix
# 
# # for 1a / 1b slide 125
# head(long$s125)
# 
# # I considered different ways of plotting this:
# conf_mat(long$s125, "origID", axis.col = "Person", axis1 = "1b", axis2 = "1a", spec.abb = sp.abb, abb.end = c("na", "nc"), key = FALSE)
# # or
# conf_mat(long$s125, "origID", axis.col = "Person", axis1 = "1a", axis2 = "1b", spec.abb = sp.abb, abb.end = c("na", "nc"), key = FALSE)
# # but these only highlight changes, not increasing accuracy, so instead I'm plotting them both against the consensus
# 
# sum(slide125$`1a` == slide125$'1b') # 183 or 61% similarity
# sum(slide125$`1a` == slide125$CID) # 198 or 66% accuracy
# sum(slide125$`1b` == slide125$CID) # 228 or 76% accuracy
# 
# png("ASFigures/Time/confusion_125_1aCon.png", 1000, 700)
# conf_mat(long$s125[long$s125$Person == "1a", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "1a", ylab = "Consensus")
# dev.off() 
# png("ASFigures/Time/confusion_125_1bCon.png", 1000, 700)
# conf_mat(long$s125[long$s125$Person == "1b", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "1b", ylab = "Consensus")
# dev.off() 
# 
# # and for 2a / 2b
# sum(slide125$`2a` == slide125$'2b') # 241 or 80% similarity
# sum(slide125$`2a` == slide125$CID) # 208 or 69% accuracy
# sum(slide125$`2b` == slide125$CID) # 237 or 79% accuracy
# 
# # again, do this as a confusion matrix
# png("ASFigures/Time/confusion_125_2aCon.png", 1000, 700)
# conf_mat(long$s125[long$s125$Person == "2a", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2a", ylab = "Consensus")
# dev.off() 
# png("ASFigures/Time/confusion_125_2bCon.png", 1000, 700)
# conf_mat(long$s125[long$s125$Person == "2b", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2b", ylab = "Consensus")
# dev.off() 
# 
# # And for 150
# sum(slide150$`1a` == slide150$'1b') # 204 or 68% similarity
# sum(slide150$`1a` == slide150$CID) # 221 or 74% accuracy
# sum(slide150$`1b` == slide150$CID) # 223 or 74% accuracy
# 
# png("ASFigures/Time/confusion_150_1aCon.png", 1000, 700)
# conf_mat(long$s150[long$s150$Person == "1a", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "1a", ylab = "Consensus")
# dev.off() 
# png("ASFigures/Time/confusion_150_1bCon.png", 1000, 700)
# conf_mat(long$s150[long$s150$Person == "1b", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "1b", ylab = "Consensus")
# dev.off() 
# 
# # 2a/2b
# sum(slide150$`2a` == slide150$'2b') # 288 or 96% similarity
# sum(slide150$`2a` == slide150$CID) # 255 or 85% accuracy
# sum(slide150$`2b` == slide150$CID) # 255 or 85% accuracy
# 
# png("ASFigures/Time/confusion_150_2aCon.png", 1000, 700)
# conf_mat(long$s150[long$s150$Person == "2a", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2a", ylab = "Consensus")
# dev.off()
# png("ASFigures/Time/confusion_150_2bCon.png", 1000, 700)
# conf_mat(long$s150[long$s150$Person == "2b", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2b", ylab = "Consensus")
# dev.off()
# 
# # 6. Digital vs. slides ---------------------------------------------------
# 
# # 6a. Individual comparisons ----------------------------------------------
# # Table 6
# long$f125 <- rbind(long$s125, long$d125)
# 
# # 125 2b vs. A
# sum(full.125$`2b` == full.125$'A') # 171 or 57% similarity
# sum(full.125$`2b` == full.125$CID.x) # 237 or 79% accuracy
# sum(full.125$`A` == full.125$CID.y) # 181 or 60% accuracy
# 
# png("ASFigures/DigitalSlide/confusion_125_2bA.png", 1000, 700)
# conf_mat(long$f125, "origID", axis.col = "Person", axis1 = "2b", axis2 = "A", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
# dev.off() 
# 
# png("ASFigures/DigitalSlide/confusion_125_2bCon.png", 1000, 700)
# conf_mat(long$f125[long$f125$Person == "2b", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2b", ylab = "Consensus")
# dev.off() 
# png("ASFigures/DigitalSlide/confusion_125_ACon.png", 1000, 700)
# conf_mat(long$f125[long$f125$Person == "A", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "A", ylab = "Consensus")
# dev.off() 
# 
# # 125 6 vs. F
# sum(full.125$`6` == full.125$'F') # 187 or 62% similarity
# sum(full.125$`6` == full.125$CID.x) # 211 or 70% accuracy
# sum(full.125$`F` == full.125$CID.y) # 226 or 75% accuracy
# 
# png("ASFigures/DigitalSlide/confusion_125_6F.png", 1000, 700)
# conf_mat(long$f125, "origID", axis.col = "Person", axis1 = "6", axis2 = "F", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
# dev.off() 
# png("ASFigures/DigitalSlide/confusion_125_FCon.png", 1000, 700)
# conf_mat(long$f125[long$f125$Person == "F", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "F", ylab = "Consensus")
# dev.off() 
# png("ASFigures/DigitalSlide/confusion_125_6Con.png", 1000, 700)
# conf_mat(long$f125[long$f125$Person == "6", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "6", ylab = "Consensus")
# dev.off() 
# 
# # 125 9 vs. G
# sum(full.125$`9` == full.125$'G') # 232 or 77% similarity
# sum(full.125$`9` == full.125$CID.x) # 206 or 69% accuracy
# sum(full.125$`G` == full.125$CID.y) # 157 or 52% accuracy
# 
# png("ASFigures/DigitalSlide/confusion_125_9G.png", 1000, 700)
# conf_mat(long$f125, "origID", axis.col = "Person", axis1 = "9", axis2 = "G", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
# dev.off() 
# 
# png("ASFigures/DigitalSlide/confusion_125_GCon.png", 1000, 700)
# conf_mat(long$f125[long$f125$Person == "G", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "G", ylab = "Consensus")
# dev.off() 
# png("ASFigures/DigitalSlide/confusion_125_9Con.png", 1000, 700)
# conf_mat(long$f125[long$f125$Person == "9", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "9", ylab = "Consensus")
# dev.off() 
# 
# # for 150
# long$f150 <- rbind(long$s150, long$d150)
# 
# # 150 2b vs. A
# sum(full.150$`2b` == full.150$'A') # 232 or 77% similarity
# sum(full.150$`2b` == full.150$CID.x) # 255 or 85% accuracy
# sum(full.150$`A` == full.150$CID.y) # 246 or 82% accuracy
# 
# png("ASFigures/DigitalSlide/confusion_150_2bA.png", 1000, 700)
# conf_mat(long$f150, "origID", axis.col = "Person", axis1 = "2b", axis2 = "A", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
# dev.off() 
# 
# png("ASFigures/DigitalSlide/confusion_150_2bCon.png", 1000, 700)
# conf_mat(long$f150[long$f150$Person == "2b", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2b", ylab = "Consensus")
# dev.off() 
# png("ASFigures/DigitalSlide/confusion_150_ACon.png", 1000, 700)
# conf_mat(long$f150[long$f150$Person == "A", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "A", ylab = "Consensus")
# dev.off() 
# 
# # 150 6 vs. F
# sum(full.150$`6` == full.150$'F') # 216 or 72% similarity
# sum(full.150$`6` == full.150$CID.x) # 244 or 81% accuracy
# sum(full.150$`F` == full.150$CID.y) # 242 or 81% accuracy
# 
# png("ASFigures/DigitalSlide/confusion_150_6F.png", 1000, 700)
# conf_mat(long$f150, "origID", axis.col = "Person", axis1 = "6", axis2 = "F", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
# dev.off() 
# 
# png("ASFigures/DigitalSlide/confusion_150_FCon.png", 1000, 700)
# conf_mat(long$f150[long$f150$Person == "F", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "F", ylab = "Consensus")
# dev.off() 
# png("ASFigures/DigitalSlide/confusion_150_6Con.png", 1000, 700)
# conf_mat(long$f150[long$f150$Person == "6", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "6", ylab = "Consensus")
# dev.off() 
# 
# # 150 9 vs. G
# sum(full.150$`9` == full.150$'G') # 212 or 71% similarity
# sum(full.150$`9` == full.150$CID.x) # 220 or 73% accuracy
# sum(full.150$`G` == full.150$CID.y) # 149 or 50% accuracy
# 
# png("ASFigures/DigitalSlide/confusion_150_9G.png", 1000, 700)
# conf_mat(long$f150, "origID", axis.col = "Person", axis1 = "9", axis2 = "G", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
# dev.off() 
# 
# png("ASFigures/DigitalSlide/confusion_150_GCon.png", 1000, 700)
# conf_mat(long$f150[long$f150$Person == "G", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "G", ylab = "Consensus")
# dev.off() 
# png("ASFigures/DigitalSlide/confusion_150_9Con.png", 1000, 700)
# conf_mat(long$f150[long$f150$Person == "9", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "9", ylab = "Consensus")
# dev.off() 
# 
# # 6b. Consensus comparisons -----------------------------------------------
# 
# sum(full.125$CID.x == full.125$CID.y) # 234 or 78% accuracy
# sum(full.150$CID.x == full.150$CID.y) # 248 or 83% accuracy
# 
# # plotting the consensus' against each other
# png("ASFigures/DigitalSlide/confusion_125_Con.png", 1000, 700)
# conf_mat(long$f125, "CID", axis.col = "Person", axis1 = "A", axis2 = "1a", spec.abb = sp.abb, abb.end = c("na", "nc"), xlab = "Digital", ylab = "Slide")
# text(-0.4, 1, "Consensus 125", cex = 1.5)
# dev.off()
# 
# png("ASFigures/DigitalSlide/confusion_150_Con.png", 1000, 700)
# conf_mat(long$f150, "CID", axis.col = "Person", axis1 = "A", axis2 = "1a", spec.abb = sp.abb, abb.end = c("na", "nc"), xlab = "Digital", ylab = "Slide")
# text(-0.4, 1, "Consensus 150", cex = 1.5)
# dev.off()
# 
# # 7. SST ------------------------------------------------------------------
# # this was only done 150 size fraction.
# # it is also not perfect as it currently uses Nadia's consensus values not mine. But that's because I can't currently rerun the ANN analysis. 
# # Figure 6
# head(divTemp)
# row.nam <- list()
# row.nam$div <- which(nchar(divTemp$Person) < 5)
# row.nam$s125 <- which(nchar(divTemp$Person) < 5 & divTemp$Analysis == "Slide" & divTemp$Size == 125)
# row.nam$s150 <- which(nchar(divTemp$Person) < 5 & divTemp$Analysis == "Slide" & divTemp$Size == 150)
# row.nam$d125 <- which(nchar(divTemp$Person) < 5 & divTemp$Analysis == "Digital" & divTemp$Size == 125)
# row.nam$d150 <- which(nchar(divTemp$Person) < 5 & divTemp$Analysis == "Digital" & divTemp$Size == 150)
# row.nam$s125c <- which(divTemp$Person == "consensus" & divTemp$Analysis == "Slide" & divTemp$Size == 125)
# row.nam$s150c <- which(divTemp$Person == "consensus" & divTemp$Analysis == "Slide" & divTemp$Size == 150)
# row.nam$d125c <- which(divTemp$Person == "consensus" & divTemp$Analysis == "Digital" & divTemp$Size == 125)
# row.nam$d150c <- which(divTemp$Person == "consensus" & divTemp$Analysis == "Digital" & divTemp$Size == 150)
# row.nam$sSST <- which(divTemp$Analysis == "Slide" & divTemp$Size == 150)
# row.nam$dSST <- which(divTemp$Analysis == "Digital" & divTemp$Size == 150)
# 
# png("ASFigures/Fig6_SST.png", 900, 800)
# par(mfrow = c(2, 1), mar = c(2.5, 4.1, .5, 1))
# tmp <- divTemp[row.nam$sSST,]
# plot(1:17, tmp$SST10m[tmp$Person != "consensus"], pch = 16, type = "n", xaxt = "n", xlab = "Person", ylab = expression(paste("SST / ", degree, "C")), ylim = c(20, 24))
# axis(1, at = 1:17, labels = tmp$Person[tmp$Person != "consensus"])
# with(tmp[tmp$Person == "consensus", ], rect(0, SST10m - SD, 18, SST10m + SD, col = rgb(0, .8, .2, alpha = .5)))
# abline(h = tmp$SST10m[tmp$Person == "consensus"], lty = 4)
# points(1:17, tmp$SST10m[match(accuracySlide$PersonID, tmp$Person)], pch = 16)
# with(tmp[match(accuracySlide$PersonID, tmp$Person), ], err_bar(SST10m, SD, 1:17))
# with(tmp, lines(c(1:2), c(SST10m[Person == "1a"], SST10m[Person == "1b"])))
# with(tmp, lines(c(3:4), c(SST10m[Person == "2a"], SST10m[Person == "2b"])))
# text(16.5, 20, "Slide 150", cex = 1.5)
# abline(h = 21.76, col = 4)
# text(16.5, 21.65, "WOA 1998", cex = 1.3, col = 4)
# 
# tmp <- divTemp[row.nam$dSST,]
# plot(1:9, tmp$SST10m[tmp$Person != "consensus"], pch = 16, type = "n", xaxt = "n", xlab = "Person", ylab = expression(paste("SST / ", degree, "C")), ylim = c(20, 24))
# axis(1, at = 1:9, labels = tmp$Person[tmp$Person != "consensus"])
# with(tmp[tmp$Person == "consensus", ], rect(0, SST10m - SD, 18, SST10m + SD, col = rgb(0, .8, .2, alpha = .5)))
# abline(h = tmp$SST10m[tmp$Person == "consensus"], lty = 4)
# points(1:9, tmp$SST10m[tmp$Person != "consensus"], pch = 16)
# with(tmp[tmp$Person != "consensus", ], err_bar(SST10m, SD, 1:9))
# text(8.75, 20, "Digital 150", cex = 1.5)
# abline(h = 21.76, col = 4)
# text(8.75, 21.65, "WOA 1998", cex = 1.3, col = 4)
# 
# par(mfrow = c(1,1))
# dev.off()
# rm(tmp)
# 
# png("ASFigures/Fig6_SST_comb.png", 800, 500)
# par(mar = c(5.1, 5.1, 4.1, 2.1))
# with(divTemp[divTemp$Size == 150,], plot(1:26, SST10m[match(ord.div, Person)], pch = 16, xaxt = "n", xlab = "Participant", ylab = expression(paste("SST / ", degree, "C")), col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(20, 24), cex.lab = 1.5, las = 1, cex.axis = 1.1))
# axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
# with(divTemp[c(row.nam$s150c, row.nam$d150c),], abline(h = c(SST10m - SD, SST10m + SD), col = ((Analysis != "Slide")*3 + 1), lty = 4))
# with(divTemp[c(row.nam$s150c, row.nam$d150c),], abline(h = SST10m, col = ((Analysis != "Slide")*3 + 1)))
# with(divTemp[divTemp$Size == 150,], err_bar(SST10m[match(ord.div, Person)], SD[match(ord.div, Person)], 1:26, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
# abline(h = 21.76, col = "green4")
# text(25, 21.65, "WOA 1998", cex = 1.3, col = "green4")
# legend("topleft", legend = c("Slide 150", "Digital 150"), pch = 16, col = c(1, 4))
# par(mar = c(5.1, 4.1, 4.1, 2.1))
# dev.off()
# 
# # 8. Diversity ------------------------------------------------------------
# 
# # 8a. Calculating diversity -----------------------------------------------
# # Table 5
# # compare the results I calculate with those that Nadia has
# head(slide125sp)
# 
# # these all included 'na' (and 'nc') as a species (to get the match. )
# # what do I get if I calculate these without those 
# # richness Slide 125
# divTemp$IF_Richness <- NA
# divTemp$IF_Richness[row.nam$s125] <- specnumber(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))[1:17] # richness
# divTemp$IF_Richness[row.nam$s125c] <- specnumber(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))["CID"]
# 
# # ShannonWiener
# divTemp$IF_ShannonWiener <- NA
# divTemp$IF_ShannonWiener[row.nam$s125] <- diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))[1:17] 
# divTemp$IF_ShannonWiener[row.nam$s125c] <- diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))["CID"]
# 
# # Dominance
# divTemp$IF_Dominance <- NA
# divTemp$IF_Dominance[row.nam$s125] <- (1 - diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]), index = "simpson"))[1:17]  
# divTemp$IF_Dominance[row.nam$s125c] <- (1 - diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]), index = "simpson"))["CID"]  
# 
# # Evenness
# divTemp$IF_Evenness <- NA
# divTemp$IF_Evenness[row.nam$s125] <- (exp(diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))) / specnumber(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)])))[1:17]
# divTemp$IF_Evenness[row.nam$s125c] <- (exp(diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))) / specnumber(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)])))["CID"]
# 
# 
# # Slide 150
# # richness 
# divTemp$IF_Richness[row.nam$s150] <- specnumber(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]))[1:17] # richness
# divTemp$IF_Richness[row.nam$s150c] <- specnumber(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]))["CID"]
# 
# # ShannonWiener
# divTemp$IF_ShannonWiener[row.nam$s150] <- diversity(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]))[1:17] 
# divTemp$IF_ShannonWiener[row.nam$s150c] <- diversity(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]))["CID"]
# 
# # Dominance
# divTemp$IF_Dominance[row.nam$s150] <- (1 - diversity(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]), index = "simpson"))[1:17]  
# divTemp$IF_Dominance[row.nam$s150c] <- (1 - diversity(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]), index = "simpson"))["CID"]  
# 
# # Evenness
# divTemp$IF_Evenness[row.nam$s150] <- (exp(diversity(t(slide150sp[3:nrow(slide125sp),2:ncol(slide150sp)]))) / specnumber(t(slide150sp[3:nrow(slide125sp),2:ncol(slide150sp)])))[1:17]
# divTemp$IF_Evenness[row.nam$s150c] <- (exp(diversity(t(slide150sp[3:nrow(slide125sp),2:ncol(slide150sp)]))) / specnumber(t(slide150sp[3:nrow(slide125sp),2:ncol(slide150sp)])))["CID"]
# 
# # Digital 125
# # richness 
# divTemp$IF_Richness[row.nam$d125] <- specnumber(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]))[1:9] # richness
# divTemp$IF_Richness[row.nam$d125c] <- specnumber(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]))["CID"]
# 
# # ShannonWiener
# divTemp$IF_ShannonWiener[row.nam$d125] <- diversity(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]))[1:9] 
# divTemp$IF_ShannonWiener[row.nam$d125c] <- diversity(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]))["CID"]
# 
# # Dominance
# divTemp$IF_Dominance[row.nam$d125] <- (1 - diversity(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]), index = "simpson"))[1:9]  
# divTemp$IF_Dominance[row.nam$d125c] <- (1 - diversity(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]), index = "simpson"))["CID"]  
# 
# # Evenness
# divTemp$IF_Evenness[row.nam$d125] <- (exp(diversity(t(digital125sp[3:nrow(slide125sp),2:ncol(digital125sp)]))) / specnumber(t(digital125sp[3:nrow(slide125sp),2:ncol(digital125sp)])))[1:9]
# divTemp$IF_Evenness[row.nam$d125c] <- (exp(diversity(t(digital125sp[3:nrow(slide125sp),2:ncol(digital125sp)]))) / specnumber(t(digital125sp[3:nrow(slide125sp),2:ncol(digital125sp)])))["CID"]
# 
# # Digital 150
# # richness 
# divTemp$IF_Richness[row.nam$d150] <- specnumber(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]))[1:9] # richness
# divTemp$IF_Richness[row.nam$d150c] <- specnumber(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]))["CID"]
# 
# # ShannonWiener
# divTemp$IF_ShannonWiener[row.nam$d150] <- diversity(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]))[1:9] 
# divTemp$IF_ShannonWiener[row.nam$d150c] <- diversity(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]))["CID"]
# 
# # Dominance
# divTemp$IF_Dominance[row.nam$d150] <- (1 - diversity(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]), index = "simpson"))[1:9]  
# divTemp$IF_Dominance[row.nam$d150c] <- (1 - diversity(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]), index = "simpson"))["CID"]  
# 
# # Evenness
# divTemp$IF_Evenness[row.nam$d150] <- (exp(diversity(t(digital150sp[3:nrow(slide125sp),2:ncol(digital150sp)]))) / specnumber(t(digital150sp[3:nrow(slide125sp),2:ncol(digital150sp)])))[1:9]
# divTemp$IF_Evenness[row.nam$d150c] <- (exp(diversity(t(digital150sp[3:nrow(slide125sp),2:ncol(digital150sp)]))) / specnumber(t(digital150sp[3:nrow(slide125sp),2:ncol(digital150sp)])))["CID"]
# 
# # 8b. Plotting diversity --------------------------------------------------
# ord.div <- c("1a", "1b", "2a", "2b", "A", 3:6, "F", 7:9, "G", 10:15, LETTERS[2:5], LETTERS[8:9])
# 
# # Figure 7
# png("ASFigures/Fig7_richness.png", 800, 500)
# # 125
# par(mar = c(5.1, 5.1, 4.1, 2.1))
# with(divTemp[divTemp$Size == 125,], plot(1:26, IF_Richness[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Richness", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(14, 30), cex.lab = 1.5, las = 1, cex.axis = 1.1))
# axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
# with(divTemp[c(row.nam$s125c, row.nam$d125c),], abline(h = IF_Richness, lty = 2, col = ((Analysis != "Slide")*3 + 1)))
# # 150
# with(divTemp[divTemp$Size == 150,], points(1:26, IF_Richness[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
# with(divTemp[c(row.nam$s150c, row.nam$d150c),], abline(h = IF_Richness, col = ((Analysis != "Slide")*3 + 1)))
# legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
# par(mar = c(5.1, 4.1, 4.1, 2.1))
# dev.off()
# 
# # comparing the values
# plot(seq(2, 18, by = 2), sort(divTemp$IF_Richness[divTemp$Analysis == "Digital" & divTemp$Size == 125 & divTemp$Person != "consensus"]), type = "b", col = "blue", lty = 2, ylim = c(14, 30))
# points(seq(2, 18, by = 2), sort(divTemp$IF_Richness[divTemp$Analysis == "Digital" & divTemp$Size == 150 & divTemp$Person != "consensus"]), type = "b", col = "blue")
# points(1:17, sort(divTemp$IF_Richness[divTemp$Analysis == "Slide" & divTemp$Size == 150 & divTemp$Person != "consensus"]), type = "b")
# points(1:17, sort(divTemp$IF_Richness[divTemp$Analysis == "Slide" & divTemp$Size == 125 & divTemp$Person != "consensus"]), type = "b", lty = 2)
# abline(h = divTemp$IF_Richness[row.nam$d125c], lty = 2, col = "blue")
# abline(h = divTemp$IF_Richness[row.nam$d150c], col = "blue")
# abline(h = divTemp$IF_Richness[row.nam$s125c], lty = 2)
# abline(h = divTemp$IF_Richness[row.nam$s150c])
# 
# png("ASFigures/Fig7_Dominance.png", 800, 500)
# # 125
# par(mar = c(5.1, 5.1, 4.1, 2.1))
# with(divTemp[divTemp$Size == 125,], plot(1:26, IF_Dominance[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Dominance", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(0.1, 0.22), cex.lab = 1.5, las = 1, cex.axis = 1.1))
# axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
# with(divTemp[c(row.nam$s125c, row.nam$d125c),], abline(h = IF_Dominance, lty = 2, col = ((Analysis != "Slide")*3 + 1)))
# # 150
# with(divTemp[divTemp$Size == 150,], points(1:26, IF_Dominance[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
# with(divTemp[c(row.nam$s150c, row.nam$d150c),], abline(h = IF_Dominance, col = ((Analysis != "Slide")*3 + 1)))
# legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
# par(mar = c(5.1, 4.1, 4.1, 2.1))
# dev.off()
# 
# png("ASFigures/Fig7_ShannonWiener.png", 800, 500)
# # 125
# par(mar = c(5.1, 5.1, 4.1, 2.1))
# with(divTemp[divTemp$Size == 125,], plot(1:26, IF_ShannonWiener[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Shannon-Wiener Diversity", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(1.95, 2.6), cex.lab = 1.5, las = 1, cex.axis = 1.1))
# axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
# with(divTemp[c(row.nam$s125c, row.nam$d125c),], abline(h = IF_ShannonWiener, lty = 2, col = ((Analysis != "Slide")*3 + 1)))
# # 150
# with(divTemp[divTemp$Size == 150,], points(1:26, IF_ShannonWiener[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
# with(divTemp[c(row.nam$s150c, row.nam$d150c),], abline(h = IF_ShannonWiener, col = ((Analysis != "Slide")*3 + 1)))
# legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
# par(mar = c(5.1, 4.1, 4.1, 2.1))
# dev.off()
# 
# png("ASFigures/Fig7_Evenness.png", 800, 500)
# # 125
# par(mar = c(5.1, 5.1, 4.1, 2.1))
# with(divTemp[divTemp$Size == 125,], plot(1:26, IF_Evenness[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Evenness", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(0.35, 0.6), cex.lab = 1.5, las = 1, cex.axis = 1.1))
# axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
# with(divTemp[c(row.nam$s125c, row.nam$d125c),], abline(h = IF_Evenness, lty = 2, col = ((Analysis != "Slide")*3 + 1)))
# # 150
# with(divTemp[divTemp$Size == 150,], points(1:26, IF_Evenness[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
# with(divTemp[c(row.nam$s150c, row.nam$d150c),], abline(h = IF_Evenness, col = ((Analysis != "Slide")*3 + 1)))
# legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
# par(mar = c(5.1, 4.1, 4.1, 2.1))
# dev.off()
# 
# # 8c. Shannon Wiener comparison -----------------------------------------
# png("ASFigures/Fig8_SWcomp.png")
# # if I run this using Nadia's data, I get similar results, however, this should be plotted with my data.
# plot(divTemp$IF_ShannonWiener[c(row.nam$s125, row.nam$d125)], divTemp$IF_ShannonWiener[c(row.nam$s150, row.nam$d150)], pch = 16, xlab = "Shannon-Wiener 125", ylab = "Shannon-Wiener 150")
# summary(lm(divTemp$IF_ShannonWiener[c(row.nam$s150, row.nam$d150)] ~ divTemp$IF_ShannonWiener[c(row.nam$s125, row.nam$d125)])) # r2 = 0.1558
# 
# abline(lm(divTemp$IF_ShannonWiener[c(row.nam$s150, row.nam$d150)] ~ divTemp$IF_ShannonWiener[c(row.nam$s125, row.nam$d125)]))
# 
# # Figure 8
# points(divTemp$IF_ShannonWiener[divTemp$Size == 125 & divTemp$Person %in% c("10", "13", "15", "B")], divTemp$IF_ShannonWiener[divTemp$Size == 150 & divTemp$Person %in% c("10", "13", "15", "B")], pch = 16, xlab = "Shannon-Wiener 125", ylab = "Shannon-Wiener 150", col = "red")
# 
# 
# summary(lm(divTemp$IF_ShannonWiener[divTemp$Size == 150 & !(divTemp$Person %in% c("10", "13", "15", "B"))] ~ divTemp$IF_ShannonWiener[divTemp$Size == 125 & !(divTemp$Person %in% c("10", "13", "15", "B"))]))
# 
# abline(lm(divTemp$IF_ShannonWiener[divTemp$Size == 150 & !(divTemp$Person %in% c("10", "13", "15", "B"))] ~ divTemp$IF_ShannonWiener[divTemp$Size == 125 & !(divTemp$Person %in% c("10", "13", "15", "B"))]), lty = 2)
# dev.off()
# 
# # the other variables show mixed results
# # richness r2 = 0.1849, p = 0.0283*
# summary(lm(divTemp$IF_Richness[c(row.nam$s150, row.nam$d150)] ~ divTemp$IF_Richness[c(row.nam$s125, row.nam$d125)]))
# # evenness r2 = 0.1891, p = 0.0264*
# summary(lm(divTemp$IF_Evenness[c(row.nam$s150, row.nam$d150)] ~ divTemp$IF_Evenness[c(row.nam$s125, row.nam$d125)]))
# # dominance r2 = 0.0551, p = 0.248
# summary(lm(divTemp$IF_Dominance[c(row.nam$s150, row.nam$d150)] ~ divTemp$IF_Dominance[c(row.nam$s125, row.nam$d125)]))
# 
# # 8d. Sensitivity to alphabetical order ---------------------------------
# # checking diversity
# specnumber(table(slide125$CID)) # 22
# specnumber(table(slide125$CIDR)) # 22
# 
# # ShannonWiener
# diversity(table(slide125$CID)) # 2.28
# diversity(table(slide125$CIDR)) # 2.26
# 
# # Dominance
# 1 - diversity(table(slide125$CID), index = "simpson") # 0.133
# 1 - diversity(table(slide125$CIDR), index = "simpson") # 0.138
# 
# # Evenness
# exp(diversity(table(slide125$CID))) / specnumber(table(slide125$CID)) # 0.445
# exp(diversity(table(slide125$CIDR))) / specnumber(table(slide125$CIDR)) # 0.440
# 
# # slide 150
# specnumber(table(slide150$CID)) # 18
# specnumber(table(slide150$CIDR)) # 19
# 
# # ShannonWiener
# diversity(table(slide150$CID)) # 2.189
# diversity(table(slide150$CIDR)) # 2.206
# 
# # Dominance
# 1 - diversity(table(slide150$CID), index = "simpson") # 0.166
# 1 - diversity(table(slide150$CIDR), index = "simpson") # 0.165
# 
# # Evenness
# exp(diversity(table(slide150$CID))) / specnumber(table(slide150$CID)) # 0.496
# exp(diversity(table(slide150$CIDR))) / specnumber(table(slide150$CIDR)) # 0.478
# 
# # digital 125
# specnumber(table(digital125$CID)) # 23 
# specnumber(table(digital125$CIDR)) # 21
# 
# # ShannonWiener
# diversity(table(digital125$CID)) # 2.257
# diversity(table(digital125$CIDR)) # 2.233
# 
# # Dominance
# 1 - diversity(table(digital125$CID), index = "simpson") # 0.141
# 1 - diversity(table(digital125$CIDR), index = "simpson") # 0.145
# 
# # Evenness
# exp(diversity(table(digital125$CID))) / specnumber(table(digital125$CID)) # 0.415
# exp(diversity(table(digital125$CIDR))) / specnumber(table(digital125$CIDR)) # 0.444
# 
# # digital 150
# specnumber(table(digital150$CID)) # 21
# specnumber(table(digital150$CIDR)) # 21
# 
# # ShannonWiener
# diversity(table(digital150$CID)) # 2.303
# diversity(table(digital150$CIDR)) # 2.289
# 
# # Dominance
# 1 - diversity(table(digital150$CID), index = "simpson") # 0.151
# 1 - diversity(table(digital150$CIDR), index = "simpson") # 0.156
# 
# # Evenness
# exp(diversity(table(digital150$CID))) / specnumber(table(digital150$CID)) # 0.476
# exp(diversity(table(digital150$CIDR))) / specnumber(table(digital150$CIDR)) # 0.470
# 
# 
# # 9. Outliers -------------------------------------------------------------
# # Generate a dataframe of outliers. 
# # Table 7
# 
# # 9a. Outliers for each analysis ------------------------------------------
# # I decided not to try to do this as it was done here. Instead, I will rank species based on their distance from the mean
# outliers <- data.frame(PersonID = accuracy$PersonID, Analysis = accuracy$Analysis)
# 
# # rank for MDS 125
# outliers$MDS125 <- NA
# tmp.pt <- nmds$IF125f$points["SCID", ]
# tmp.tab <- sqrt((nmds$IF125f$points[,1] - tmp.pt[1])^2 + (nmds$IF125f$points[,2] - tmp.pt[2])^2)
# outliers$MDS125[outliers$Analysis == "Slide"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))]))
# tmp.pt <- nmds$IF125f$points["DCID", ]
# tmp.tab <- sqrt((nmds$IF125f$points[,1] - tmp.pt[1])^2 + (nmds$IF125f$points[,2] - tmp.pt[2])^2)
# outliers$MDS125[outliers$Analysis == "Digital"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))]))
# rm(tmp.pt, tmp.tab)
# 
# # and mds 150
# outliers$MDS150 <- NA
# tmp.pt <- nmds$IF150f$points["SCID", ]
# tmp.tab <- sqrt((nmds$IF150f$points[,1] - tmp.pt[1])^2 + (nmds$IF150f$points[,2] - tmp.pt[2])^2)
# outliers$MDS150[outliers$Analysis == "Slide"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))]))
# tmp.pt <- nmds$IF150f$points["DCID", ]
# tmp.tab <- sqrt((nmds$IF150f$points[,1] - tmp.pt[1])^2 + (nmds$IF150f$points[,2] - tmp.pt[2])^2)
# outliers$MDS150[outliers$Analysis == "Digital"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))]))
# rm(tmp.pt, tmp.tab)
# 
# # for percentage accuracy
# outliers$IF_PtAc125 <- NA
# outliers$IF_PtAc125[outliers$Analysis == "Slide"][order(100-accuracySlide$IF_PtAc125[match(outliers$PersonID[outliers$Analysis == "Slide"], accuracySlide$PersonID)])] <- sort(rank(100-accuracySlide$IF_PtAc125))
# outliers$IF_PtAc125[outliers$Analysis == "Digital"][order(100-accuracyDigital$IF_PtAc125[match(outliers$PersonID[outliers$Analysis == "Digital"], accuracyDigital$PersonID)])] <- sort(rank(100-accuracyDigital$IF_PtAc125))
# 
# outliers$IF_PtAc150 <- NA
# outliers$IF_PtAc150[outliers$Analysis == "Slide"][order(100-accuracySlide$IF_PtAc150[match(outliers$PersonID[outliers$Analysis == "Slide"], accuracySlide$PersonID)])] <- sort(rank(100-accuracySlide$IF_PtAc150))
# outliers$IF_PtAc150[outliers$Analysis == "Digital"][order(100-accuracyDigital$IF_PtAc150[match(outliers$PersonID[outliers$Analysis == "Digital"], accuracyDigital$PersonID)])] <- sort(rank(100-accuracyDigital$IF_PtAc150))
# 
# # for mean pairwise agreement
# # for percentage accuracy
# outliers$mnPA125 <- NA
# outliers$mnPA125[outliers$Analysis == "Slide"][order(100-accuracySlide$mnPA125[match(outliers$PersonID[outliers$Analysis == "Slide"], accuracySlide$PersonID)])] <- sort(rank(100-accuracySlide$mnPA125))
# outliers$mnPA125[outliers$Analysis == "Digital"][order(100-accuracyDigital$mnPA125[match(outliers$PersonID[outliers$Analysis == "Digital"], accuracyDigital$PersonID)])] <- sort(rank(100-accuracyDigital$mnPA125))
# 
# outliers$mnPA150 <- NA
# outliers$mnPA150[outliers$Analysis == "Slide"][order(100-accuracySlide$mnPA150[match(outliers$PersonID[outliers$Analysis == "Slide"], accuracySlide$PersonID)])] <- sort(rank(100-accuracySlide$mnPA150))
# outliers$mnPA150[outliers$Analysis == "Digital"][order(100-accuracyDigital$mnPA150[match(outliers$PersonID[outliers$Analysis == "Digital"], accuracyDigital$PersonID)])] <- sort(rank(100-accuracyDigital$mnPA150))
# 
# # SST
# outliers$SST <- NA
# tmp.pt <- abs(divTemp$SST10m[row.nam$s150] - 21.76)
# names(tmp.pt) <- divTemp$Person[row.nam$s150]
# outliers$SST[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# tmp.pt <- abs(divTemp$SST10m[row.nam$d150] - 21.76)
# names(tmp.pt) <- divTemp$Person[row.nam$d150]
# outliers$SST[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# rm(tmp.pt)
# 
# # richness
# outliers$IF_Richness125 <- NA
# tmp.pt <- abs(divTemp$IF_Richness[row.nam$s125] - divTemp$IF_Richness[row.nam$s125c])
# names(tmp.pt) <- divTemp$Person[row.nam$s125]
# outliers$IF_Richness125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# tmp.pt <- abs(divTemp$IF_Richness[row.nam$d125] - divTemp$IF_Richness[row.nam$d125c])
# names(tmp.pt) <- divTemp$Person[row.nam$d125]
# outliers$IF_Richness125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# outliers$IF_Richness150 <- NA
# tmp.pt <- abs(divTemp$IF_Richness[row.nam$s150] - divTemp$IF_Richness[row.nam$s150c])
# names(tmp.pt) <- divTemp$Person[row.nam$s150]
# outliers$IF_Richness150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# tmp.pt <- abs(divTemp$IF_Richness[row.nam$d150] - divTemp$IF_Richness[row.nam$d150c])
# names(tmp.pt) <- divTemp$Person[row.nam$d150]
# outliers$IF_Richness150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# rm(tmp.pt)
# 
# # Dominance
# outliers$IF_Dominance125 <- NA
# tmp.pt <- abs(divTemp$IF_Dominance[row.nam$s125] - divTemp$IF_Dominance[row.nam$s125c])
# names(tmp.pt) <- divTemp$Person[row.nam$s125]
# outliers$IF_Dominance125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# tmp.pt <- abs(divTemp$IF_Dominance[row.nam$d125] - divTemp$IF_Dominance[row.nam$d125c])
# names(tmp.pt) <- divTemp$Person[row.nam$d125]
# outliers$IF_Dominance125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# outliers$IF_Dominance150 <- NA
# tmp.pt <- abs(divTemp$IF_Dominance[row.nam$s150] - divTemp$IF_Dominance[row.nam$s150c])
# names(tmp.pt) <- divTemp$Person[row.nam$s150]
# outliers$IF_Dominance150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# tmp.pt <- abs(divTemp$IF_Dominance[row.nam$d150] - divTemp$IF_Dominance[row.nam$d150c])
# names(tmp.pt) <- divTemp$Person[row.nam$d150]
# outliers$IF_Dominance150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# rm(tmp.pt)
# 
# # Evenness
# outliers$IF_Evenness125 <- NA
# tmp.pt <- abs(divTemp$IF_Evenness[row.nam$s125] - divTemp$IF_Evenness[row.nam$s125c])
# names(tmp.pt) <- divTemp$Person[row.nam$s125]
# outliers$IF_Evenness125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# tmp.pt <- abs(divTemp$IF_Evenness[row.nam$d125] - divTemp$IF_Evenness[row.nam$d125c])
# names(tmp.pt) <- divTemp$Person[row.nam$d125]
# outliers$IF_Evenness125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# outliers$IF_Evenness150 <- NA
# tmp.pt <- abs(divTemp$IF_Evenness[row.nam$s150] - divTemp$IF_Evenness[row.nam$s150c])
# names(tmp.pt) <- divTemp$Person[row.nam$s150]
# outliers$IF_Evenness150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# tmp.pt <- abs(divTemp$IF_Evenness[row.nam$d150] - divTemp$IF_Evenness[row.nam$d150c])
# names(tmp.pt) <- divTemp$Person[row.nam$d150]
# outliers$IF_Evenness150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# rm(tmp.pt)
# 
# # ShannonWiener
# outliers$IF_ShannonWiener125 <- NA
# tmp.pt <- abs(divTemp$IF_ShannonWiener[row.nam$s125] - divTemp$IF_ShannonWiener[row.nam$s125c])
# names(tmp.pt) <- divTemp$Person[row.nam$s125]
# outliers$IF_ShannonWiener125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# tmp.pt <- abs(divTemp$IF_ShannonWiener[row.nam$d125] - divTemp$IF_ShannonWiener[row.nam$d125c])
# names(tmp.pt) <- divTemp$Person[row.nam$d125]
# outliers$IF_ShannonWiener125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# outliers$IF_ShannonWiener150 <- NA
# tmp.pt <- abs(divTemp$IF_ShannonWiener[row.nam$s150] - divTemp$IF_ShannonWiener[row.nam$s150c])
# names(tmp.pt) <- divTemp$Person[row.nam$s150]
# outliers$IF_ShannonWiener150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# tmp.pt <- abs(divTemp$IF_ShannonWiener[row.nam$d150] - divTemp$IF_ShannonWiener[row.nam$d150c])
# names(tmp.pt) <- divTemp$Person[row.nam$d150]
# outliers$IF_ShannonWiener150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
# rm(tmp.pt)
# 
# # weighted sums
# outliers$fullSum <- rowSums(outliers[, 3:17])
# 
# outliers$wtSum <- rowSums(outliers[, grep("MDS", names(outliers))])/2 + rowSums(outliers[, grep("Pt", names(outliers))])/2 + outliers$SST + rowSums(outliers[, grep("Richness", names(outliers))])/8 + rowSums(outliers[, grep("Dominance", names(outliers))])/8 + rowSums(outliers[, grep("Evenness", names(outliers))])/8 + rowSums(outliers[, grep("ShannonWiener", names(outliers))])/8
# 
# outliers$wtSumPt <- outliers$wtSum / c(rep(4*17, 17), rep(4*9, 9)) * 100
# 
# # Table 8 
# # add percentage of specimens identified to accuracy
# accuracySlide$ptID125 <- apply(slide125[,col.nam$s125], 2, function(x) sum(x != "na") / 300 * 100)[accuracySlide$PersonID]
# accuracySlide$ptID150 <- apply(slide150[,col.nam$s150], 2, function(x) sum(x != "na") / 300 * 100)[accuracySlide$PersonID]
# accuracyDigital$ptID125 <- apply(digital125[,col.nam$d125], 2, function(x) sum(x != "na") / 300 * 100)[accuracyDigital$PersonID]
# accuracyDigital$ptID150 <- apply(digital150[,col.nam$d150], 2, function(x) sum(x != "na") / 300 * 100)[accuracyDigital$PersonID]
# accuracy <- rbind(accuracySlide, accuracyDigital)
# accuracy$Analysis <- c(rep("Slide", 17), rep("Digital", 9))
# 
# # looking for correlations
# pairs(outliers[, 3:ncol(outliers)])
# pairs(outliers[, grep("125", names(outliers))])
# pairs(outliers[, grep("150", names(outliers))])
# 
# # Create a .csv file for the accuracy data
# # create a subset of the data for this
# tmp.sub <- accuracy[, c("PersonID", "Analysis", "Experience", "Routine", "IF_PtAc125", "IF_PtAc150", "ptID125", "ptID150")]
# tmp.sub[, grep("125|150", names(tmp.sub))] <- round(tmp.sub[, grep("125|150", names(tmp.sub))], 2)
# write.csv(tmp.sub, "Outputs/Accuracy.csv", row.names = FALSE)
# rm(tmp.sub)
# 
# # add the person level info to this
# tmp.s <- people[!is.na(people$SlideID), !grepl("Digital", names(people))]
# tmp.d <- people[!is.na(people$DigitalID), !grepl("Slide", names(people))]
# names(tmp.s)[1] <- names(tmp.d)[1] <- "PersonID"
# tmp.s <- rbind(tmp.s, tmp.s[!is.na(tmp.s$ExperienceSlideB), ])
# tmp.s$PersonID[tmp.s$PersonID %in% tmp.s$PersonID[duplicated(tmp.s$PersonID)]] <- paste(tmp.s$PersonID[tmp.s$PersonID %in% tmp.s$PersonID[duplicated(tmp.s$PersonID)]], rep(c("a", "b"), each = 2), sep = "")
# tmp.s$ExperienceSlideA[grep("b", tmp.s$PersonID)] <- tmp.s$ExperienceSlideB[grep("b", tmp.s$PersonID)]
# tmp.s <- tmp.s[, names(tmp.s) != "ExperienceSlideB"]
# names(tmp.d) <- names(tmp.s) <- gsub("SlideA", "", names(tmp.s))
# tmp <- rbind(tmp.s, tmp.d)
# tmp.out <- merge(outliers, tmp)
# tmp.out <- merge(tmp.out, accuracy[, c("PersonID", "ptID125", "ptID150")])
# 
# #par(ask = TRUE)
# par(mfrow = c(6, 3))
# par(mar = c(2,2,1,1))
# for (i in names(tmp.out)[c(21:25, 27:28)]) {
#   for (j in names(tmp.out)[3:20]) 
#     plot(factor(tmp.out[tmp.out$Analysis == "Slide",i]), tmp.out[tmp.out$Analysis == "Slide",j], main = j, col = "red")
#   for (j in names(tmp.out)[3:20]) 
#     plot(factor(tmp.out[tmp.out$Analysis == "Digital",i]), tmp.out[tmp.out$Analysis == "Digital",j], main = j, col = "blue")
# }
# rm(i, j)
# #par(ask = FALSE)
# par(mfrow = c(1,1))
# par(mar = c(5.1, 4.1, 4.1, 2.1))
# rm(tmp.s, tmp.out, tmp.d, tmp)
# 
# 
# # 9b. Outliers for digital and slide combined -----------------------------
# # I decided not to try to do this as it was done here. Instead, I will rank species based on their distance from the mean
# # rank for MDS 125
# outliers$cMDS125 <- NA
# tmp.pt <- nmds$IF125f$points["SCID", ]
# tmp.tab <- sqrt((nmds$IF125f$points[,1] - tmp.pt[1])^2 + (nmds$IF125f$points[,2] - tmp.pt[2])^2)
# tmp.pt <- nmds$IF125f$points["DCID", ]
# tmp.tab[grep("^[A-M]", names(tmp.tab))] <- sqrt((nmds$IF125f$points[,1] - tmp.pt[1])^2 + (nmds$IF125f$points[,2] - tmp.pt[2])^2)[grep("^[A-M]", names(tmp.tab))]
# outliers$cMDS125[order(tmp.tab[match(outliers$PersonID, names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID, names(tmp.tab))]))
# rm(tmp.pt, tmp.tab)
# 
# # and mds 150
# outliers$cMDS150 <- NA
# tmp.pt <- nmds$IF150f$points["SCID", ]
# tmp.tab <- sqrt((nmds$IF150f$points[,1] - tmp.pt[1])^2 + (nmds$IF150f$points[,2] - tmp.pt[2])^2)
# tmp.pt <- nmds$IF150f$points["DCID", ]
# tmp.tab[grep("^[A-M]", names(tmp.tab))] <- sqrt((nmds$IF150f$points[,1] - tmp.pt[1])^2 + (nmds$IF150f$points[,2] - tmp.pt[2])^2)[grep("^[A-M]", names(tmp.tab))]
# outliers$cMDS150[order(tmp.tab[match(outliers$PersonID, names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID, names(tmp.tab))]))
# rm(tmp.pt, tmp.tab)
# 
# # for percentage accuracy
# outliers$cIF_PtAc125 <- NA
# outliers$cIF_PtAc125[order(100-accuracy$IF_PtAc125[match(outliers$PersonID, accuracy$PersonID)])] <- sort(rank(100-accuracy$IF_PtAc125))
# 
# outliers$cIF_PtAc150 <- NA
# outliers$cIF_PtAc150[order(100-accuracy$IF_PtAc150[match(outliers$PersonID, accuracy$PersonID)])] <- sort(rank(100-accuracy$IF_PtAc150))
# 
# # for mean pairwise agreement
# # for percentage accuracy
# outliers$cmnPA125 <- NA
# outliers$cmnPA125[order(100-accuracy$mnPA125[match(outliers$PersonID, accuracy$PersonID)])] <- sort(rank(100-accuracy$mnPA125))
# 
# outliers$cmnPA150 <- NA
# outliers$cmnPA150[order(100-accuracy$mnPA150[match(outliers$PersonID, accuracy$PersonID)])] <- sort(rank(100-accuracy$mnPA150))
# 
# # SST
# outliers$cSST <- NA
# tmp.pt <- abs(divTemp$SST10m[c(row.nam$s150, row.nam$d150)] - 21.76)
# names(tmp.pt) <- divTemp$Person[c(row.nam$s150, row.nam$d150)]
# outliers$cSST[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
# rm(tmp.pt)
# 
# # richness
# outliers$cIF_Richness125 <- NA
# tmp.pt <- abs(divTemp$IF_Richness[c(row.nam$s125, row.nam$d125)] - divTemp$IF_Richness[row.nam$s125c])
# names(tmp.pt) <- divTemp$Person[c(row.nam$s125, row.nam$d125)]
# tmp.pt[grep("^[A-M]", names(tmp.pt))] <- abs(divTemp$IF_Richness[row.nam$d125] - divTemp$IF_Richness[row.nam$d125c])
# outliers$cIF_Richness125[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# outliers$cIF_Richness150 <- NA
# tmp.pt <- abs(divTemp$IF_Richness[c(row.nam$s150, row.nam$d150)] - divTemp$IF_Richness[row.nam$s150c])
# names(tmp.pt) <- divTemp$Person[c(row.nam$s150, row.nam$d150)]
# tmp.pt[grep("^[A-M]", names(tmp.pt))] <- abs(divTemp$IF_Richness[row.nam$d150] - divTemp$IF_Richness[row.nam$d150c])
# outliers$cIF_Richness150[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
# rm(tmp.pt)
# 
# # Dominance
# outliers$cIF_Dominance125 <- NA
# tmp.pt <- abs(divTemp$IF_Dominance[c(row.nam$s125, row.nam$d125)] - divTemp$IF_Dominance[row.nam$s125c])
# names(tmp.pt) <- divTemp$Person[c(row.nam$s125, row.nam$d125)]
# tmp.pt[grep("^[A-M]", names(tmp.pt))] <- abs(divTemp$IF_Dominance[row.nam$d125] - divTemp$IF_Dominance[row.nam$d125c])
# outliers$cIF_Dominance125[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# outliers$cIF_Dominance150 <- NA
# tmp.pt <- abs(divTemp$IF_Dominance[c(row.nam$s150, row.nam$d150)] - divTemp$IF_Dominance[row.nam$s150c])
# names(tmp.pt) <- divTemp$Person[c(row.nam$s150, row.nam$d150)]
# tmp.pt[grep("^[A-M]", names(tmp.pt))] <- abs(divTemp$IF_Dominance[row.nam$d150] - divTemp$IF_Dominance[row.nam$d150c])
# outliers$cIF_Dominance150[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
# rm(tmp.pt)
# 
# # Evenness
# outliers$cIF_Evenness125 <- NA
# tmp.pt <- abs(divTemp$IF_Evenness[c(row.nam$s125, row.nam$d125)] - divTemp$IF_Evenness[row.nam$s125c])
# names(tmp.pt) <- divTemp$Person[c(row.nam$s125, row.nam$d125)]
# tmp.pt[grep("^[A-M]", names(tmp.pt))] <- abs(divTemp$IF_Evenness[row.nam$d125] - divTemp$IF_Evenness[row.nam$d125c])
# outliers$cIF_Evenness125[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# outliers$cIF_Evenness150 <- NA
# tmp.pt <- abs(divTemp$IF_Evenness[c(row.nam$s150, row.nam$d150)] - divTemp$IF_Evenness[row.nam$s150c])
# names(tmp.pt) <- divTemp$Person[c(row.nam$s150, row.nam$d150)]
# tmp.pt[grep("^[A-M]", names(tmp.pt))] <- abs(divTemp$IF_Evenness[row.nam$d150] - divTemp$IF_Evenness[row.nam$d150c])
# outliers$cIF_Evenness150[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
# rm(tmp.pt)
# 
# # ShannonWiener
# outliers$cIF_ShannonWiener125 <- NA
# tmp.pt <- abs(divTemp$IF_ShannonWiener[c(row.nam$s125, row.nam$d125)] - divTemp$IF_ShannonWiener[row.nam$s125c])
# names(tmp.pt) <- divTemp$Person[c(row.nam$s125, row.nam$d125)]
# tmp.pt[grep("^[A-M]", names(tmp.pt))] <- abs(divTemp$IF_ShannonWiener[row.nam$d125] - divTemp$IF_ShannonWiener[row.nam$d125c])
# outliers$cIF_ShannonWiener125[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
# 
# outliers$cIF_ShannonWiener150 <- NA
# tmp.pt <- abs(divTemp$IF_ShannonWiener[c(row.nam$s150, row.nam$d150)] - divTemp$IF_ShannonWiener[row.nam$s150c])
# names(tmp.pt) <- divTemp$Person[c(row.nam$s150, row.nam$d150)]
# tmp.pt[grep("^[A-M]", names(tmp.pt))] <- abs(divTemp$IF_ShannonWiener[row.nam$d150] - divTemp$IF_ShannonWiener[row.nam$d150c])
# outliers$cIF_ShannonWiener150[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
# rm(tmp.pt)
# 
# # weighted sums
# outliers$cfullSum <- rowSums(outliers[, grep("^c", names(outliers))])
# 
# outliers$cwtSum <- rowSums(outliers[, grep("cMDS", names(outliers))])/2 + rowSums(outliers[, grep("^c.*Pt", names(outliers))])/2 + outliers$cSST + rowSums(outliers[, grep("cIF_Richness", names(outliers))])/8 + rowSums(outliers[, grep("cIF_Dominance", names(outliers))])/8 + rowSums(outliers[, grep("cIF_Evenness", names(outliers))])/8 + rowSums(outliers[, grep("cIF_ShannonWiener", names(outliers))])/8
# 
# outliers$cwtSumPt <- outliers$cwtSum / 26 * 100
# 
# # looking for correlations
# pairs(outliers[, grep("^c", names(outliers))])
# pairs(outliers[, grep("^c.*125", names(outliers))])
# pairs(outliers[, grep("^c.*150", names(outliers))])
# 
# # add the person level info to this
# tmp.s <- people[!is.na(people$SlideID), !grepl("Digital", names(people))]
# tmp.d <- people[!is.na(people$DigitalID), !grepl("Slide", names(people))]
# names(tmp.s)[1] <- names(tmp.d)[1] <- "PersonID"
# tmp.s <- rbind(tmp.s, tmp.s[!is.na(tmp.s$ExperienceSlideB), ])
# tmp.s$PersonID[tmp.s$PersonID %in% tmp.s$PersonID[duplicated(tmp.s$PersonID)]] <- paste(tmp.s$PersonID[tmp.s$PersonID %in% tmp.s$PersonID[duplicated(tmp.s$PersonID)]], rep(c("a", "b"), each = 2), sep = "")
# tmp.s$ExperienceSlideA[grep("b", tmp.s$PersonID)] <- tmp.s$ExperienceSlideB[grep("b", tmp.s$PersonID)]
# tmp.s <- tmp.s[, names(tmp.s) != "ExperienceSlideB"]
# names(tmp.d) <- names(tmp.s) <- gsub("SlideA", "", names(tmp.s))
# tmp <- rbind(tmp.s, tmp.d)
# tmp.out <- merge(outliers, tmp)
# tmp.out <- merge(tmp.out, accuracy[, c("PersonID", "ptID125", "ptID150")])
# 
# #par(ask = TRUE)
# par(mfrow = c(6, 3))
# par(mar = c(2,2,1,1))
# for (i in names(tmp.out)[c(39:43, 45:45)]) {
#   for (j in grep("^c", names(tmp.out), value = TRUE)) 
#     plot(factor(tmp.out[,i]), tmp.out[,j], main = j, col = "red")
# }
# rm(i, j)
# #par(ask = FALSE)
# par(mfrow = c(1,1))
# par(mar = c(5.1, 4.1, 4.1, 2.1))
# 
# plot(tmp.out$Analysis, tmp.out$cwtSumPt)
# 
# rm(tmp.s, tmp.out, tmp.d, tmp)
# 
# # 10. Size vs. maximum agreement ------------------------------------------
# # Figure 5
# png("ASFigures/Fig5_size_agreement_125.png")
# with(size125, plot(slideAgreement, Length, pch = 16, main = "> 125", las = 1, xlab = "Agreement", ylab = expression(paste("Maximum diamter / ", mu, "m"))))
# lines(names(tapply(size125$Length, size125$slideAgreement, max)), tapply(size125$Length, size125$slideAgreement, max), pch = 16)
# with(size125, points(digitalAgreement, Length, pch = 16, col = "blue"))
# lines(names(tapply(size125$Length, size125$digitalAgreement, max)), tapply(size125$Length, size125$digitalAgreement, max), pch = 16, col = 4)
# legend("topleft", col = c(1, 4), pch = 16, legend = c("Slide", "Digital"))
# dev.off()
# 
# png("ASFigures/Fig5_size_agreement_150.png")
# with(size150, plot(slideAgreement, Length, pch = 16, main = "> 150", las = 1, xlab = "Agreement", ylab = expression(paste("Maximum diamter / ", mu, "m"))))
# lines(names(tapply(size150$Length, size150$slideAgreement, max)), tapply(size150$Length, size150$slideAgreement, max), pch = 16)
# with(size150, points(digitalAgreement, Length, pch = 16, col = "blue"))
# lines(names(tapply(size150$Length, size150$digitalAgreement, max)), tapply(size150$Length, size150$digitalAgreement, max), pch = 16, col = 4)
# legend("topleft", col = c(1, 4), pch = 16, legend = c("Slide", "Digital"))
# dev.off()
# 
# # ex. figure 7
# # run a  linear model
# lm_size_125 <- lm(size125$slideAgreement ~ size125$Length)
# par(mfrow = c(2,2))
# plot(lm_size_125)
# par(mfrow = c(1,1))
# 
# # try logging to see if that helps
# lm_size_125 <- lm(size125$slideAgreement ~ log(size125$Length))
# par(mfrow = c(2,2))
# plot(lm_size_125)
# par(mfrow = c(1,1))
# 
# # slightly better
# with(size125, plot(log(Length), slideAgreement, pch = 16))
# abline(lm_size_125)
#  # however this doesn't really make sense, as it predicts slide agreement above 100%. Given the y-value is bounded between 0 and 100 (or 0 and 1), I think this should be a binomial model. 
# 
# glm_size_125 <- glm(slideAgreement/100 ~ log(Length), data = size125, family = "binomial")
# par(mfrow = c(2,2))
# plot(glm_size_125) # still not a good QQ plot
# par(mfrow = c(1,1))
# 
# png("ASFigures/exFig7_sizeAgg125_glm.png")
# with(size125, plot(Length, slideAgreement/100, pch = 16))
# pred <- predict(glm_size_125, newdata = data.frame(Length = 100:600), type = "response")
# points(100:600, pred, type = "l")
# dev.off()
# rm(pred)
# 
# glm_size_150 <- glm(slideAgreement/100 ~ log(Length), data = size150, family = "binomial")
# par(mfrow = c(2,2))
# plot(glm_size_150) # still not a good QQ plot
# par(mfrow = c(1,1))
# 
# png("ASFigures/exFig7_sizeAgg150_glm.png")
# with(size150, plot(Length, slideAgreement/100, pch = 16))
# pred <- predict(glm_size_150, newdata = data.frame(Length = 100:700), type = "response")
# points(100:700, pred, type = "l")
# dev.off()
# rm(pred)
# 
# # is size correlated with abundance?
# # what is the mean species size
# sp.size <- list()
# sp.size$ssC <- tapply(c(size125$Length, size150$Length), c(size125$Con1, size150$Con1), mean)
# sp.size$ss125 <-tapply(size125$Length, size125$Con1, mean)
# sp.size$ss150 <-tapply(size150$Length, size150$Con1, mean)
# 
# sp.abun <- list()
# sp.abun$ssC <- table(c(size125$Con1, size150$Con1))
# sp.abun$ss125 <- table(size125$Con1)
# sp.abun$ss150 <- table(size150$Con1)
# 
# # plot species mean against abundance
# plot(as.numeric(sp.size$ssC), log(as.numeric(sp.abun$ssC)), pch = 16)
# points(as.numeric(sp.size$ss125), log(as.numeric(sp.abun$ss125)), pch = 16, col = "blue")
# points(as.numeric(sp.size$ss150), log(as.numeric(sp.abun$ss150)), pch = 16, col = "red")
# abline(lm(log(as.numeric(sp.abun$ssC)) ~ as.numeric(sp.size$ssC)))
# abline(lm(log(as.numeric(sp.abun$ss125)) ~ as.numeric(sp.size$ss125)), col = "blue")
# abline(lm(log(as.numeric(sp.abun$ss150)) ~ as.numeric(sp.size$ss150)), col = "red")
# 
# sort(sp.size)
# 
# # 11. Comparison of different tests ---------------------------------------
# # ex. Figure 12
# png("ASFigures/exFig12_Consensus frequency.png", 500, 700)
# par(mfrow = c(2, 1))
# # for slide
# barplot(t(cbind(table(slide125$MaxCon), table(slide150$MaxCon))), beside = TRUE, xlab = "Maximum consensus", legend.text = c("125", "150"), args.legend = c(x = 7, y = 70), main = "Slide")
# abline(v = 18.5, lty = 2)
# 
# # for digital
# barplot(t(cbind(table(digital125$MaxCon), table(digital150$MaxCon))), beside = TRUE, xlab = "Maximum consensus", legend.text = c("125", "150"), args.legend = c(x = 4.25, y = 70), main = "Digital")
# abline(v = 9.5, lty = 2)
# par(mfrow = c(1, 1))
# dev.off()
# 
# # accuracy as percentage
# table(slide125$MaxCon)/300 * 100
# table(slide150$MaxCon)/300 * 100
# table(digital125$MaxCon)/300 * 100
# table(digital150$MaxCon)/300 * 100
# 
# # fraction with > 50% aggreement
# sum(slide125$MaxCon > c50_cutoff$slide)/300 * 100 # 72.3%
# sum(slide150$MaxCon > c50_cutoff$slide)/300 * 100 # 85.3%
# sum(digital125$MaxCon > c50_cutoff$digital)/300 * 100 # 66.3%
# sum(digital150$MaxCon > c50_cutoff$digital)/300 * 100 # 76.0%
# 
# # sum of unidentified
# sum(slide125$MaxCon <= c50_cutoff$slide) # 83
# sum(slide150$MaxCon <= c50_cutoff$slide) # 44
# sum(digital125$MaxCon <= c50_cutoff$digital) # 101
# sum(digital150$MaxCon <= c50_cutoff$digital) # 72
# 
# # 12. Interpreting the diversity metrics ----------------------------------
# 
# # 12a. Comparison with ForCenS -----------------------------------------------
# # see how the datasets compare with the ForCenS values for that latitude
# ForCenSred <- as.data.frame(read_excel("Data/Siccha_ForCenS.xlsx", sheet = "ForCenSred", na = "N/A")) 
# ForCenSred[, 22:62][is.na(ForCenSred[, 22:62])] <- 0
# str(ForCenSred)
# 
# png("ASFigures/Div_cf_ForCenSred.png", 550, 600)
# par(mfrow = c(2,2))
# # species richness
# tmp.rich <- specnumber(ForCenSred[, 22:62]) # richness(ForCenSred[,22:62])
# plot(ForCenSred$Latitude, tmp.rich, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Richness")
# points(rep(30.2, nrow(divTemp)), divTemp$IF_Richness, col = "blue", pch = 16)
# points(rep(30.2, 4), divTemp$IF_Richness[divTemp$Person == "consensus"], col = "red", pch = 16)
# 
# # ShannonWiener
# tmp.sw <- diversity(ForCenSred[,22:62])
# plot(ForCenSred$Latitude, tmp.sw, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Shannon Wiener")
# points(rep(30.2, nrow(divTemp)), divTemp$IF_ShannonWiener, col = "blue", pch = 16)
# points(rep(30.2, 4), divTemp$IF_ShannonWiener[divTemp$Person == "consensus"], col = "red", pch = 16)
# 
# # Dominance
# tmp.dom <- (1 - diversity(ForCenSred[,22:62], index = "simpson"))
# plot(ForCenSred$Latitude, tmp.dom, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Dominance")
# points(rep(30.2, nrow(divTemp)), divTemp$IF_Dominance, col = "blue", pch = 16)
# points(rep(30.2, 4), divTemp$IF_Dominance[divTemp$Person == "consensus"], col = "red", pch = 16)
# 
# # Evenness
# tmp.eve <- (exp(diversity(ForCenSred[,22:62])) / specnumber(ForCenSred[,22:62]))
# plot(ForCenSred$Latitude, tmp.eve, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Evenness")
# points(rep(30.2, nrow(divTemp)), divTemp$IF_Evenness, col = "blue", pch = 16)
# points(rep(30.2, 4), divTemp$IF_Evenness[divTemp$Person == "consensus"], col = "red", pch = 16)
# par(mfrow = c(1,1))
# dev.off()
# 
# rm(tmp, i, tmp.rich, tmp.sw, tmp.dom, tmp.eve, ForCenSred)
# 
# 
# # 12b. Simulations for studying changes -----------------------------------
# tmp <- data.frame("a" = c(10, 10, 0, 20, 12, 10, 10, 5, 9, 10, 8, 8, 10, 11), "b" = c(10, 10, 10, 0, 10, 0, 10, 10, 10, 10, 12, 10, 10, 10), "c" = c(2, 2, 2, 2, 0, 12, 4, 2, 2, 1, 2, 4, 3, 1), "d" = c(2, 0, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 1, 2), "e" = c(0, 2, 10, 0, 0, 0, 0, 5, 1, 1, 0, 0, 0, 0))
# rownames(tmp) <- c("Orig", "ChangeSpListR", "ChangeSpListC", "LumpCc", "LumpCr", "LumpRc", "LumpRr", "SplitCe", "SplitCu", "SplitR", "BoundShiftCc", "BoundShift2Rc", "BoundShift3Rr", "BoundShift4Cr")
# tmp$Rich <- specnumber(tmp[, 1:5]) # richness
# tmp$SW <- diversity(tmp[, 1:5]) # Shannon-Wiener
# tmp$Dom <- (1 - diversity(tmp[, 1:5], index = "simpson")) # dominance
# tmp$Eve <- exp(diversity(tmp[, 1:5])) / specnumber(tmp[, 1:5]) # evenness
# tmp
# rm(tmp)
# 
# 
# # 12c. Direction of change ------------------------------------------------
# # add some columns to the diversity dataframe to investigate the magnitude of the change. 
# divTemp$Cen_IFE <- divTemp$Cen_IFD <- divTemp$Cen_IFSW <- divTemp$Cen_IFR <- divTemp$Cen_SST <- NA
# divTemp[row.nam$s125, grep("Cen", names(divTemp))] <- sweep(data.matrix(divTemp[row.nam$s125, c(8, 10:13)]), 2, as.numeric(divTemp[row.nam$s125c, c(8, 10:13)]))
# 
# divTemp[row.nam$s150, grep("Cen", names(divTemp))] <- sweep(data.matrix(divTemp[row.nam$s150, c(8, 10:13)]), 2, as.numeric(divTemp[row.nam$s150c, c(8, 10:13)]))
# 
# divTemp[row.nam$d125, grep("Cen", names(divTemp))] <- sweep(data.matrix(divTemp[row.nam$d125, c(8, 10:13)]), 2, as.numeric(divTemp[row.nam$d125c, c(8, 10:13)]))
# 
# divTemp[row.nam$d150, grep("Cen", names(divTemp))] <- sweep(data.matrix(divTemp[row.nam$d150, c(8, 10:13)]), 2, as.numeric(divTemp[row.nam$d150c, c(8, 10:13)]))
# 
# pairs(divTemp[, grep("Cen_", names(divTemp))], col = factor(paste(divTemp$Analysis, divTemp$Size, sep = "_")), pch = 16)
# 
# # add some columns to the diversity dataframe to investigate the direction of the change. 
# divTemp$Dir_IFE <- divTemp$Dir_IFD <- divTemp$Dir_IFSW <- divTemp$Dir_IFR <- divTemp$Dir_SST <- NA
# divTemp[row.nam$s125, grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[row.nam$s125, c(8, 10:13)]), 2, as.numeric(divTemp[row.nam$s125c, c(8, 10:13)])) > 0, 1, -1)
# 
# divTemp[row.nam$s150, grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[row.nam$s150, c(8, 10:13)]), 2, as.numeric(divTemp[row.nam$s150c, c(8, 10:13)])) > 0, 1, -1)
# 
# divTemp[row.nam$d125, grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[row.nam$d125, c(8, 10:13)]), 2, as.numeric(divTemp[row.nam$d125c, c(8, 10:13)])) > 0, 1, -1)
# 
# divTemp[row.nam$d150, grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[row.nam$d150, c(8, 10:13)]), 2, as.numeric(divTemp[row.nam$d150c, c(8, 10:13)])) > 0, 1, -1)
# 
# table(divTemp$Dir_SST, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
# table(divTemp$Dir_IFR, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
# table(divTemp$Dir_IFSW, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
# table(divTemp$Dir_IFD, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
# table(divTemp$Dir_IFE, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
# 
# tmp.div <- divTemp[, grep("Analysis|Size|Person|SST10m|SD|IF_", names(divTemp))]
# tmp.div$ID <- paste(tmp.div$Analysis, tmp.div$Person)
# head(tmp.div)
# tmp.div <- reshape(tmp.div, direction = "wide", v.names = grep("SST10m|SD|IF_", names(tmp.div), value = TRUE), timevar = "Size", idvar = "ID")
# tmp.div <- tmp.div[, names(tmp.div) != "ID"]
# tmp.div[, grepl("1", names(tmp.div)) & !grepl("Rich", names(tmp.div))] <- round(tmp.div[, grepl("1", names(tmp.div)) & !grepl("Rich", names(tmp.div))], 2)
# write.csv(tmp.div, file = "Outputs/Diversity_temperature.csv", row.names = FALSE)
# rm(tmp.div)
# 
# # 12d. Without incomplete data --------------------------------------------
# summary(divTemp$IF_Richness[divTemp$Person %in% accuracySlide$PersonID[accuracySlide$ptID125 == 100] & divTemp$Analysis == "Slide"& divTemp$Size == 125])
# summary(divTemp$IF_Richness[divTemp$Person %in% accuracySlide$PersonID[accuracySlide$ptID150 == 100] & divTemp$Analysis == "Slide"& divTemp$Size == 150])
# summary(divTemp$IF_Richness[divTemp$Person %in% accuracyDigital$PersonID[accuracyDigital$ptID125 == 100] & divTemp$Analysis == "Digital"& divTemp$Size == 125])
# summary(divTemp$IF_Richness[divTemp$Person %in% accuracyDigital$PersonID[accuracyDigital$ptID150 == 100] & divTemp$Analysis == "Digital"& divTemp$Size == 150])
# 
# 
# 
# 
# 13. What if we use a combined consensus for the analysis? ---------------
# 13c. NMDS analysis ------------------------------------------------------
# for 125
full.125 <- merge(slide125, digital125, by = "Specimen")

trsp$CCIF125f <- data.frame(t(full.125[, !(names(full.125) %in% c("Specimen", "consensus50.x", "consensus20.x", "consensus50.y", "consensus20.y", "MaxCon.x", "CIDR.x", "MaxCon.y", "CIDR.y"))]))
rownames(trsp$CCIF125f)[nchar(rownames(trsp$CCIF125f)) >= 5] <- c("SsC", "SCID", "DsC", "DCID")

# consider the stress of the NMDS
stress$CCIF125f <- rep(NA, 10)
for (i in 1:10) {
  stress$CCIF125f[i] <- metaMDS(daisy(trsp$CCIF125f), k = i)$stress
}
plot(stress$CCIF125f, type = "b")
rm(i)
# looks like between 2 and 3 dimensions would be reasonable

# check for variation
tmp <- metaMDS(daisy(trsp$CCIF125f))
tmp$stress
plot(tmp, display = "sites", type = "t", cex = 1.5)
# some variation, so optimise

stress$CCIFop125f <- 1
# find the nmds plot with the lowest stress out of 20000 runs
# for (i in 1:1000) {
#   tmp <- metaMDS(daisy(trsp$CCIF125f))
#   if (stress$CCIFop125f > tmp$stress) {
#     nmds$CCIF125f <- tmp
#     stress$CCIFop125f <- tmp$stress
#   }
# }
# rm(i, tmp)
# save(stress, nmds, file = "Outputs/CCNMDS_125.RData")
load("Outputs/CCNMDS_125.RData")

stressplot(nmds$CCIF125f)

png("ASFigures/CombCon/IF_NMDS_125.png", 600, 600)
plot(nmds$CCIF125f, type = "n", display = "sites", cex = 1, xlab = "Axis 1", ylab = "Axis 2", las = 1, cex.lab = 1.5, cex.axis = 1.2, main = "NMDS 125")
points(nmds$CCIF125f, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
text(nmds$CCIF125f$points[nchar(rownames(nmds$CCIF125f$points)) <3, ], labels = rownames(nmds$CCIF125f$points)[nchar(rownames(nmds$CCIF125f$points)) <3])
points(nmds$CCIF125f$points[nchar(rownames(nmds$CCIF125f$points)) >2, ], pch = "+")
text(nmds$CCIF125f$points[rownames(nmds$CCIF125f$points) == "SCID", 1]+0.03, nmds$CCIF125f$points[rownames(nmds$CCIF125f$points) == "SCID", 2]+0.00, labels = "CID")
text(nmds$CCIF125f$points[rownames(nmds$CCIF125f$points) == "SsC", 1]+0.02, nmds$CCIF125f$points[rownames(nmds$CCIF125f$points) == "SsC", 2]+0.00, labels = "sC")

legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16, cex = 1.3)
dev.off()

# for 150
full.150 <- merge(slide150, digital150, by = "Specimen")
trsp$CCIF150f <- data.frame(t(full.150[, !(names(full.150) %in% c("Specimen", "consensus50.x", "consensus20.x", "consensus50.y", "consensus20.y", "MaxCon.x", "CIDR.x", "MaxCon.y", "CIDR.y"))]))
rownames(trsp$CCIF150f)[nchar(rownames(trsp$CCIF150f)) >= 5] <- c("SsC", "SCID", "DsC", "DCID")

# consider the stress of the NMDS
stress$CCIF150f <- rep(NA, 10)
for (i in 1:10) {
  stress$CCIF150f[i] <- metaMDS(daisy(trsp$CCIF150f), k = i)$stress
}
plot(stress$CCIF150f, type = "b")
rm(i)
# looks like between 2 and 3 dimensions would be reasonable, although as with Nadia's data, the breakpoint is less obvious for 150 than it is for 125.

# check for variation
tmp <- metaMDS(daisy(trsp$CCIF150f))
tmp$stress
plot(tmp, display = "sites", type = "t", cex = 1.5)
# some variation, so optimise

stress$CCIFop150f <- 1
# find the nmds plot with the lowest stress out of 1000 runs
# for (i in 1:1000) {
#   tmp <- metaMDS(daisy(trsp$CCIF150f))
#   if (stress$CCIFop150f > tmp$stress) {
#     nmds$CCIF150f <- tmp
#     stress$CCIFop150f <- tmp$stress
#   }
# }
# rm(i, tmp)
# save(stress, nmds, file = "Outputs/CCNMDS_150f.RData")
load("Outputs/CCNMDS_150f.RData")

stressplot(nmds$CCIF150f)

# the full plot
png("ASFigures/CombCon/IF_NMDS_150.png", 600, 600)
plot(nmds$CCIF150f, type = "n", display = "sites", cex = 1, xlab = "Axis 1", ylab = "Axis 2", las = 1, cex.lab = 1.5, cex.axis = 1.2, main = "NMDS 150")
points(nmds$CCIF150f, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
text(nmds$CCIF150f$points[nchar(rownames(nmds$CCIF150f$points)) <3, ], labels = rownames(nmds$CCIF150f$points)[nchar(rownames(nmds$CCIF150f$points)) <3])
points(nmds$CCIF150f$points[nchar(rownames(nmds$CCIF150f$points)) >2, ], pch = "+")
text(nmds$CCIF150f$points[rownames(nmds$CCIF150f$points) == "SCID", 1]+0.02, nmds$CCIF150f$points[rownames(nmds$CCIF150f$points) == "SCID", 2], labels = "CID")
text(nmds$CCIF150f$points[rownames(nmds$CCIF150f$points) == "SsC", 1]+0.025, nmds$CCIF150f$points[rownames(nmds$CCIF150f$points) == "SsC", 2], labels = "sC")
legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16, cex = 1.3)
dev.off()

# focussing on the main section. As noted above, it is better to run this as a new analysis rather than just zoom in, as the influence of the outliers means that the stability of the central points hasn't been tested
nmds$CCIF150z <- metaMDS(daisy(trsp$CCIF150f[!(rownames(trsp$CCIF150f) %in% c("3", "C", "E", "G")),]))

# consider the stress of the NMDS
stress$CCIF150z <- rep(NA, 10)
for (i in 1:10) {
  stress$CCIF150z[i] <- metaMDS(daisy(trsp$CCIF150f[!(rownames(trsp$CCIF150f) %in% c("3", "C", "E", "G")),]), k = i)$stress
}
plot(stress$CCIF150z, type = "b")
rm(i)
stressplot(nmds$CCIF150z)
# again between 2 and 3 dimensions is probably reasonable

# check for variation
tmp <- metaMDS(daisy(trsp$CCIF150f[!(rownames(trsp$CCIF150f) %in% c("3", "C", "E", "G")),]))
tmp$stress
plot(tmp, display = "sites", type = "t", cex = 1.5)
# some variation, so optimise

stress$CCIFop150z <- 1
# find the nmds plot with the lowest stress out of 1000 runs
# for (i in 1:1000) {
#   tmp <- metaMDS(daisy(trsp$CCIF150f[!(rownames(trsp$CCIF150f) %in% c("3", "C", "E", "G")),]))
#   if (stress$CCIFop150z > tmp$stress) {
#     nmds$CCIF150z <- tmp
#     stress$CCIFop150z <- tmp$stress
#   }
# }
# rm(i, tmp)
# save(stress, nmds, file = "Outputs/CCNMDS_150z.RData")
load("Outputs/CCNMDS_150z.RData")

stress$CCIFop150z
stressplot(nmds$CCIF150z)

# the zoomed plot
png("ASFigures/CombCon/IF_NMDS_150_zoom.png", 600, 600)
plot(nmds$CCIF150z, type = "n", display = "sites", cex = 1, xlab = "Axis 1", ylab = "Axis 2", las = 1, cex.lab = 1.5, cex.axis = 1.2, main = "NMDS 150 zoomed")
points(nmds$CCIF150z, pch = 21, cex = 4, col = mds.col$pair[!(mds.col$person %in% c("3", "C", "E", "G"))], bg = brewer.pal(5, "Set2")[mds.col$sch.col[!(mds.col$person %in% c("3", "C", "E", "G"))]], lwd = 2)
text(nmds$CCIF150z$points[nchar(rownames(nmds$CCIF150z$points)) <3, ], labels = rownames(nmds$CCIF150z$points)[nchar(rownames(nmds$CCIF150z$points)) <3])
points(nmds$CCIF150z$points[nchar(rownames(nmds$CCIF150z$points)) >2, ], pch = "+")
text(nmds$CCIF150z$points[rownames(nmds$CCIF150z$points) == "SCID", 1]+0.012, nmds$CCIF150z$points[rownames(nmds$CCIF150z$points) == "SCID", 2]+0.008, labels = "CID")
text(nmds$CCIF150z$points[rownames(nmds$CCIF150z$points) == "SsC", 1]+0.012, nmds$CCIF150z$points[rownames(nmds$CCIF150z$points) == "SsC", 2]+0.008, labels = "sC")
legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16, cex = 1.3)
dev.off()

# output the scree plots
png("ASFigures/CombCon/Scree plots.png")
par(mfrow = c(2, 2))
plot(stress$CCIF125f, type = "b", bty = "l", las = 1, xlab = "Dimensions", ylab = "Stress", main = "NMDS 125")
plot(stress$CCIF150f, type = "b", bty = "l", las = 1, xlab = "Dimensions", ylab = "Stress", main = "NMDS 150")
plot(stress$CCIF150z, type = "b", bty = "l", las = 1, xlab = "Dimensions", ylab = "Stress", main = "NMDS 150 zoomed")
par(mfrow = c(1,1))
dev.off()

# 13d. Repeated analysis by workers -----------------------------------------
# for 1a / 1b slide 125
head(combconlong$s125)

sum(slide125$`1a` == slide125$'1b') # 183 or 61% similarity
sum(slide125$`1a` == slide125$CID) # 192 or 64% accuracy (cf. 198)
sum(slide125$`1b` == slide125$CID) # 227 or 76% accuracy (cf. 228)

png("ASFigures/Combcon/Time/confusion_125_1aCon.png", 1000, 700)
conf_mat(combconlong$s125[combconlong$s125$Person == "1a", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "1a", ylab = "Consensus")
dev.off()
png("ASFigures/Combcon/Time/confusion_125_1bCon.png", 1000, 700)
conf_mat(combconlong$s125[combconlong$s125$Person == "1b", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "1b", ylab = "Consensus")
dev.off()

# and for 2a / 2b
sum(slide125$`2a` == slide125$'2b') # 241 or 80% similarity
sum(slide125$`2a` == slide125$CID) # 206 or 69% accuracy (cf. 208)
sum(slide125$`2b` == slide125$CID) # 240 or 80% accuracy (cf. 237)

# again, do this as a confusion matrix
png("ASFigures/Combcon/Time/confusion_125_2aCon.png", 1000, 700)
conf_mat(combconlong$s125[combconlong$s125$Person == "2a", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2a", ylab = "Consensus")
dev.off()
png("ASFigures/Combcon/Time/confusion_125_2bCon.png", 1000, 700)
conf_mat(combconlong$s125[combconlong$s125$Person == "2b", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2b", ylab = "Consensus")
dev.off()

# And for 150
sum(slide150$`1a` == slide150$'1b') # 204 or 68% similarity
sum(slide150$`1a` == slide150$CID) # 221 or 74% accuracy (cf. 221)
sum(slide150$`1b` == slide150$CID) # 219 or 73% accuracy (cf. 223)

png("ASFigures/Combcon/Time/confusion_150_1aCon.png", 1000, 700)
conf_mat(combconlong$s150[combconlong$s150$Person == "1a", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "1a", ylab = "Consensus")
dev.off()
png("ASFigures/Combcon/Time/confusion_150_1bCon.png", 1000, 700)
conf_mat(combconlong$s150[combconlong$s150$Person == "1b", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "1b", ylab = "Consensus")
dev.off()

# 2a/2b
sum(slide150$`2a` == slide150$'2b') # 288 or 96% similarity
sum(slide150$`2a` == slide150$CID) # 257 or 86% accuracy (cf. 255)
sum(slide150$`2b` == slide150$CID) # 257 or 86% accuracy (cf. 255)

png("ASFigures/Combcon/Time/confusion_150_2aCon.png", 1000, 700)
conf_mat(combconlong$s150[combconlong$s150$Person == "2a", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2a", ylab = "Consensus")
dev.off()
png("ASFigures/Combcon/Time/confusion_150_2bCon.png", 1000, 700)
conf_mat(combconlong$s150[combconlong$s150$Person == "2b", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2b", ylab = "Consensus")
dev.off()

# 13e. Digital vs. slides ---------------------------------------------------
combconlong$f125 <- rbind(combconlong$s125, combconlong$d125)

# 125 2b vs. A
sum(full.125$`2b` == full.125$'A') # 171 or 57% similarity
sum(full.125$`2b` == full.125$CID.x) # 240 or 80% (cf. 237)
sum(full.125$`A` == full.125$CID.y) # 175 or 58% accuracy (cf. 181)

png("ASFigures/CombCon/DigitalSlide/confusion_125_2bA.png", 1000, 700)
conf_mat(combconlong$f125, "origID", axis.col = "Person", axis1 = "2b", axis2 = "A", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
dev.off()

png("ASFigures/CombCon/DigitalSlide/confusion_125_2bCon.png", 1000, 700)
conf_mat(combconlong$f125[combconlong$f125$Person == "2b", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2b", ylab = "Consensus")
dev.off()
png("ASFigures/CombCon/DigitalSlide/confusion_125_ACon.png", 1000, 700)
conf_mat(combconlong$f125[combconlong$f125$Person == "A", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "A", ylab = "Consensus")
dev.off()

# 125 6 vs. F
sum(full.125$`6` == full.125$'F') # 187 or 62% similarity
sum(full.125$`6` == full.125$CID.x) # 212 or 71% (cf. 211)
sum(full.125$`F` == full.125$CID.y) # 229 or 76% accuracy (cf. 226)

png("ASFigures/CombCon/DigitalSlide/confusion_125_6F.png", 1000, 700)
conf_mat(combconlong$f125, "origID", axis.col = "Person", axis1 = "6", axis2 = "F", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
dev.off()
png("ASFigures/CombCon/DigitalSlide/confusion_125_FCon.png", 1000, 700)
conf_mat(combconlong$f125[combconlong$f125$Person == "F", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "F", ylab = "Consensus")
dev.off()
png("ASFigures/CombCon/DigitalSlide/confusion_125_6Con.png", 1000, 700)
conf_mat(combconlong$f125[combconlong$f125$Person == "6", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "6", ylab = "Consensus")
dev.off()

# 125 9 vs. G
sum(full.125$`9` == full.125$'G') # 232 or 77% similarity
sum(full.125$`9` == full.125$CID.x) # 212 or 71% accuracy (cf. 206)
sum(full.125$`G` == full.125$CID.y) # 160 or 53% accuracy (cf. 157)

png("ASFigures/CombCon/DigitalSlide/confusion_125_9G.png", 1000, 700)
conf_mat(combconlong$f125, "origID", axis.col = "Person", axis1 = "9", axis2 = "G", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
dev.off()

png("ASFigures/CombCon/DigitalSlide/confusion_125_GCon.png", 1000, 700)
conf_mat(combconlong$f125[combconlong$f125$Person == "G", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "G", ylab = "Consensus")
dev.off()
png("ASFigures/CombCon/DigitalSlide/confusion_125_9Con.png", 1000, 700)
conf_mat(combconlong$f125[combconlong$f125$Person == "9", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "9", ylab = "Consensus")
dev.off()

# for 150
combconlong$f150 <- rbind(combconlong$s150, combconlong$d150)

# 150 2b vs. A
sum(full.150$`2b` == full.150$'A') # 232 or 77% similarity
sum(full.150$`2b` == full.150$CID.x) # 257 or 86% accuracy (cf. 255)
sum(full.150$`A` == full.150$CID.y) # 241 or 80% accuracy  (cf. 246)

png("ASFigures/CombCon/DigitalSlide/confusion_150_2bA.png", 1000, 700)
conf_mat(combconlong$f150, "origID", axis.col = "Person", axis1 = "2b", axis2 = "A", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
dev.off()

png("ASFigures/CombCon/DigitalSlide/confusion_150_2bCon.png", 1000, 700)
conf_mat(combconlong$f150[combconlong$f150$Person == "2b", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2b", ylab = "Consensus")
dev.off()
png("ASFigures/CombCon/DigitalSlide/confusion_150_ACon.png", 1000, 700)
conf_mat(combconlong$f150[combconlong$f150$Person == "A", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "A", ylab = "Consensus")
dev.off()

# 150 6 vs. F
sum(full.150$`6` == full.150$'F') # 216 or 72% similarity
sum(full.150$`6` == full.150$CID.x) # 246 or 82% accuracy (cf. 244)
sum(full.150$`F` == full.150$CID.y) # 244 or 81% accuracy (cf. 242)

png("ASFigures/CombCon/DigitalSlide/confusion_150_6F.png", 1000, 700)
conf_mat(combconlong$f150, "origID", axis.col = "Person", axis1 = "6", axis2 = "F", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
dev.off()

png("ASFigures/CombCon/DigitalSlide/confusion_150_FCon.png", 1000, 700)
conf_mat(combconlong$f150[combconlong$f150$Person == "F", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "F", ylab = "Consensus")
dev.off()
png("ASFigures/CombCon/DigitalSlide/confusion_150_6Con.png", 1000, 700)
conf_mat(combconlong$f150[combconlong$f150$Person == "6", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "6", ylab = "Consensus")
dev.off()

# 150 9 vs. G
sum(full.150$`9` == full.150$'G') # 212 or 71% similarity
sum(full.150$`9` == full.150$CID.x) # 216 or 72% accuracy (cf. 220)
sum(full.150$`G` == full.150$CID.y) # 148 or 49% accuracy (cf. 149)

png("ASFigures/CombCon/DigitalSlide/confusion_150_9G.png", 1000, 700)
conf_mat(combconlong$f150, "origID", axis.col = "Person", axis1 = "9", axis2 = "G", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
dev.off()

png("ASFigures/CombCon/DigitalSlide/confusion_150_GCon.png", 1000, 700)
conf_mat(combconlong$f150[combconlong$f150$Person == "G", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "G", ylab = "Consensus")
dev.off()
png("ASFigures/CombCon/DigitalSlide/confusion_150_9Con.png", 1000, 700)
conf_mat(combconlong$f150[combconlong$f150$Person == "9", ], "origID", "CID", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "9", ylab = "Consensus")
dev.off()

# 13f. SST ----------------------------------------------------------------
# this was only done 150 size fraction.
# it is also not perfect as it currently uses Nadia's consensus values not mine. But that's because I can't currently rerun the ANN analysis.
# Figure 6
divTemp <- divTemp[, !grepl("Cen_|Dir_", names(divTemp))]

png("ASFigures/CombCon/Fig6_SST_comb.png", 800, 500)
par(mar = c(5.1, 5.1, 4.1, 2.1))
# points
with(divTemp[divTemp$Size == 150,], plot(1:26, SST10m[match(ord.div, Person)], pch = 16, xaxt = "n", xlab = "Participant", ylab = expression(paste("SST / ", degree, "C")), col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(20, 24), cex.lab = 1.5, las = 1, cex.axis = 1.1))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
# consensus values - this should be one line, but haven't been able to re-run yet.
with(divTemp[row.nam$s150c,], lines(x = c(0, 20.5), y = rep(SST10m - SD, each = 2), col = "green4", lty = 4))
with(divTemp[row.nam$s150c,], lines(x = c(0, 20.5), y = rep(SST10m + SD, each = 2), col = "green4", lty = 4))
with(divTemp[row.nam$d150c,], lines(x = c(20.5, 27), y = rep(SST10m - SD, each = 2), col = "green4", lty = 4))
with(divTemp[row.nam$d150c,], lines(x = c(20.5, 27), y = rep(SST10m + SD, each = 2), col = "green4", lty = 4))
with(divTemp[row.nam$s150c,], lines(x = c(0, 20.5), y = rep(SST10m, each = 2), col = "green4"))
with(divTemp[row.nam$d150c,], lines(x = c(20.5, 27), y = rep(SST10m, each = 2), col = "green4"))
# mean values
with(divTemp[row.nam$s150c,], lines(x = c(0, 20.5), y = rep(mean(divTemp$SST10m[row.nam$s150]), each = 2)))
with(divTemp[row.nam$d150c,], lines(x = c(20.5, 27), y = rep(mean(divTemp$SST10m[row.nam$d150]), each = 2), col = 4))
abline(v = 20.5, col = "grey 50")

# error bars / actual  / legend
with(divTemp[divTemp$Size == 150,], err_bar(SST10m[match(ord.div, Person)], SD[match(ord.div, Person)], 1:26, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
abline(h = 21.76, col = "green4", lwd = 2)
text(25, 21.65, "WOA 1998", cex = 1.3, col = "green4")
legend("topleft", legend = c("Slide 150", "Digital 150"), pch = 16, col = c(1, 4))
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# 13g. Diversity ------------------------------------------------------------
# Calculating diversity
# richness
divTemp$IF_Richness[row.nam$s125c] <- specnumber(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))["CID"]
divTemp$IF_Richness[row.nam$s150c] <- specnumber(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]))["CID"]
divTemp$IF_Richness[row.nam$d125c] <- specnumber(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]))["CID"]
divTemp$IF_Richness[row.nam$d150c] <- specnumber(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]))["CID"]

# ShannonWiener
divTemp$IF_ShannonWiener[row.nam$s125c] <- diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))["CID"]
divTemp$IF_ShannonWiener[row.nam$s150c] <- diversity(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]))["CID"]
divTemp$IF_ShannonWiener[row.nam$d125c] <- diversity(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]))["CID"]
divTemp$IF_ShannonWiener[row.nam$d150c] <- diversity(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]))["CID"]

# Dominance
divTemp$IF_Dominance[row.nam$s125c] <- (1 - diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]), index = "simpson"))["CID"]
divTemp$IF_Dominance[row.nam$s150c] <- (1 - diversity(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]), index = "simpson"))["CID"]
divTemp$IF_Dominance[row.nam$d125c] <- (1 - diversity(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]), index = "simpson"))["CID"]
divTemp$IF_Dominance[row.nam$d150c] <- (1 - diversity(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]), index = "simpson"))["CID"]

# Evenness
divTemp$IF_Evenness[row.nam$s125c] <- (exp(diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))) / specnumber(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)])))["CID"]
divTemp$IF_Evenness[row.nam$s150c] <- (exp(diversity(t(slide150sp[3:nrow(slide125sp),2:ncol(slide150sp)]))) / specnumber(t(slide150sp[3:nrow(slide125sp),2:ncol(slide150sp)])))["CID"]
divTemp$IF_Evenness[row.nam$d125c] <- (exp(diversity(t(digital125sp[3:nrow(slide125sp),2:ncol(digital125sp)]))) / specnumber(t(digital125sp[3:nrow(slide125sp),2:ncol(digital125sp)])))["CID"]
divTemp$IF_Evenness[row.nam$d150c] <- (exp(diversity(t(digital150sp[3:nrow(slide125sp),2:ncol(digital150sp)]))) / specnumber(t(digital150sp[3:nrow(slide125sp),2:ncol(digital150sp)])))["CID"]

# Plotting diversity
ord.div <- c("1a", "1b", "2a", "2b", "A", 3:6, "F", 7:9, "G", 10:15, LETTERS[2:5], LETTERS[8:9])

# Figure 7
png("ASFigures/CombCon/Fig7_richness.png", 800, 500)
par(mar = c(5.1, 5.1, 4.1, 2.1))
# points
with(divTemp[divTemp$Size == 125,], plot(1:26, IF_Richness[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Richness", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(14, 30), cex.lab = 1.5, las = 1, cex.axis = 1.1))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
with(divTemp[divTemp$Size == 150,], points(1:26, IF_Richness[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
with(divTemp[c(row.nam$s150c, row.nam$d150c),], abline(h = IF_Richness, col = 1))
# consensus lines
abline(h = divTemp$IF_Richness[row.nam$s125c], col = "green4", lty = 2)
abline(h = divTemp$IF_Richness[row.nam$s150c], col = "green4")
# mean values
lines(x = c(20.5, 27), y = rep(mean(divTemp$IF_Richness[row.nam$d125]), 2), col = "blue", lty = 2)
lines(x = c(0, 20.5), y = rep(mean(divTemp$IF_Richness[row.nam$s125]), 2), lty = 2)
lines(x = c(20.5, 27), y = rep(mean(divTemp$IF_Richness[row.nam$d150]), 2), col = "blue")
lines(x = c(0, 20.5), y = rep(mean(divTemp$IF_Richness[row.nam$s150]), 2))
abline(v = 20.5, col = "grey 50")
# legend
legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("ASFigures/CombCon/Fig7_Dominance.png", 800, 500)
# points
par(mar = c(5.1, 5.1, 4.1, 2.1))
with(divTemp[divTemp$Size == 125,], plot(1:26, IF_Dominance[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Dominance", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(0.1, 0.22), cex.lab = 1.5, las = 1, cex.axis = 1.1))
with(divTemp[divTemp$Size == 150,], points(1:26, IF_Dominance[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
# consensus lines
abline(h = divTemp$IF_Dominance[row.nam$s125c], col = "green4", lty = 2)
abline(h = divTemp$IF_Dominance[row.nam$s150c], col = "green4")
# mean values
lines(x = c(20.5, 27), y = rep(mean(divTemp$IF_Dominance[row.nam$d125]), 2), col = "blue", lty = 2)
lines(x = c(0, 20.5), y = rep(mean(divTemp$IF_Dominance[row.nam$s125]), 2), lty = 2)
lines(x = c(20.5, 27), y = rep(mean(divTemp$IF_Dominance[row.nam$d150]), 2), col = "blue")
lines(x = c(0, 20.5), y = rep(mean(divTemp$IF_Dominance[row.nam$s150]), 2))
abline(v = 20.5, col = "grey 50")
# legend
legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("ASFigures/CombCon/Fig7_ShannonWiener.png", 800, 500)
# 125
par(mar = c(5.1, 5.1, 4.1, 2.1))
# points
with(divTemp[divTemp$Size == 125,], plot(1:26, IF_ShannonWiener[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Shannon-Wiener Diversity", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(1.95, 2.6), cex.lab = 1.5, las = 1, cex.axis = 1.1))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
with(divTemp[divTemp$Size == 150,], points(1:26, IF_ShannonWiener[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
# consensus lines
abline(h = divTemp$IF_ShannonWiener[row.nam$s125c], col = "green4", lty = 2)
abline(h = divTemp$IF_ShannonWiener[row.nam$s150c], col = "green4")
# mean values
lines(x = c(20.5, 27), y = rep(mean(divTemp$IF_ShannonWiener[row.nam$d125]), 2), col = "blue", lty = 2)
lines(x = c(0, 20.5), y = rep(mean(divTemp$IF_ShannonWiener[row.nam$s125]), 2), lty = 2)
lines(x = c(20.5, 27), y = rep(mean(divTemp$IF_ShannonWiener[row.nam$d150]), 2), col = "blue")
lines(x = c(0, 20.5), y = rep(mean(divTemp$IF_ShannonWiener[row.nam$s150]), 2))
abline(v = 20.5, col = "grey 50")
# legend
legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("ASFigures/CombCon/Fig7_Evenness.png", 800, 500)
par(mar = c(5.1, 5.1, 4.1, 2.1))
with(divTemp[divTemp$Size == 125,], plot(1:26, IF_Evenness[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Evenness", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(0.35, 0.6), cex.lab = 1.5, las = 1, cex.axis = 1.1))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
with(divTemp[divTemp$Size == 150,], points(1:26, IF_Evenness[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
# consensus lines
abline(h = divTemp$IF_Evenness[row.nam$s125c], col = "green4", lty = 2)
abline(h = divTemp$IF_Evenness[row.nam$s150c], col = "green4")
# mean values
lines(x = c(20.5, 27), y = rep(mean(divTemp$IF_Evenness[row.nam$d125]), 2), col = "blue", lty = 2)
lines(x = c(0, 20.5), y = rep(mean(divTemp$IF_Evenness[row.nam$s125]), 2), lty = 2)
lines(x = c(20.5, 27), y = rep(mean(divTemp$IF_Evenness[row.nam$d150]), 2), col = "blue")
lines(x = c(0, 20.5), y = rep(mean(divTemp$IF_Evenness[row.nam$s150]), 2))
abline(v = 20.5, col = "grey 50")
# legend
legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# Sensitivity to alphabetical order
# checking diversity
# 125 (both slide and digital are the same)
specnumber(table(slide125$CID)) # 22 (s - same, d - 23)
specnumber(table(slide125$CIDR)) # 22 (s - same, d - 21)

# ShannonWiener
diversity(table(slide125$CID)) # 2.28 (s - same, d - 2.257)
diversity(table(slide125$CIDR)) # 2.26 (s - same, d - 2.233)

# Dominance
1 - diversity(table(slide125$CID), index = "simpson") # 0.134 (s - 0.133, d - 0.141)
1 - diversity(table(slide125$CIDR), index = "simpson") # 0.137 (s - 0.138, d - 0.145)

# Evenness
exp(diversity(table(slide125$CID))) / specnumber(table(slide125$CID)) # 0.446 (s - 0.445, d - 0.415)
exp(diversity(table(slide125$CIDR))) / specnumber(table(slide125$CIDR)) # 0.436 (s - 0.440, d - 0.444)

# slide 150
specnumber(table(slide150$CID)) # 20 (s - 18, d - 21)
specnumber(table(slide150$CIDR)) # 20 (s - 19, d - 21)

# ShannonWiener
diversity(table(slide150$CID)) # 2.213 (s - 2.189, d - 2.303)
diversity(table(slide150$CIDR)) # 2.213 (s - 2.206, d - 2.289)

# Dominance
1 - diversity(table(slide150$CID), index = "simpson") # 0.164 (s - 0.166, d - 0.151)
1 - diversity(table(slide150$CIDR), index = "simpson") # 0.166 (s - 0.165, d - 0.156)

# Evenness
exp(diversity(table(slide150$CID))) / specnumber(table(slide150$CID)) # 0.457 (s - 0.496, d - 0.476)
exp(diversity(table(slide150$CIDR))) / specnumber(table(slide150$CIDR)) # 0.457 (s - 0.478, d - 0.470)

tmp.div <- divTemp[, grep("Analysis|Size|Person|SST10m|SD|IF_", names(divTemp))]
tmp.div$ID <- paste(tmp.div$Analysis, tmp.div$Person)
head(tmp.div)
tmp.div <- reshape(tmp.div, direction = "wide", v.names = grep("SST10m|SD|IF_", names(tmp.div), value = TRUE), timevar = "Size", idvar = "ID")
tmp.div <- tmp.div[, names(tmp.div) != "ID"]
tmp.div[, grepl("1", names(tmp.div)) & !grepl("Rich", names(tmp.div))] <- round(tmp.div[, grepl("1", names(tmp.div)) & !grepl("Rich", names(tmp.div))], 2)
write.csv(tmp.div, file = "Outputs/combconDiversity_temperature.csv", row.names = FALSE)
rm(tmp.div)

# add some columns to the diversity dataframe to investigate the direction of the change.
divTemp$Dir_IFE <- divTemp$Dir_IFD <- divTemp$Dir_IFSW <- divTemp$Dir_IFR <- divTemp$Dir_SST <- NA
divTemp[row.nam$s125, grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[row.nam$s125, c(8, 10:13)]), 2, as.numeric(divTemp[row.nam$s125c, c(8, 10:13)])) > 0, 1, -1)

divTemp[row.nam$s150, grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[row.nam$s150, c(8, 10:13)]), 2, as.numeric(divTemp[row.nam$s150c, c(8, 10:13)])) > 0, 1, -1)

divTemp[row.nam$d125, grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[row.nam$d125, c(8, 10:13)]), 2, as.numeric(divTemp[row.nam$d125c, c(8, 10:13)])) > 0, 1, -1)

divTemp[row.nam$d150, grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[row.nam$d150, c(8, 10:13)]), 2, as.numeric(divTemp[row.nam$d150c, c(8, 10:13)])) > 0, 1, -1)


table(divTemp$Dir_SST, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
table(divTemp$Dir_IFR, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
table(divTemp$Dir_IFSW, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
table(divTemp$Dir_IFD, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
table(divTemp$Dir_IFE, paste(divTemp$Analysis, divTemp$Size, sep = "_"))

table(divTemp$Dir_IFR, divTemp$Dir_IFE, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
table(divTemp$Dir_IFR, divTemp$Dir_IFSW, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
table(divTemp$Dir_IFR, divTemp$Dir_IFD, paste(divTemp$Analysis, divTemp$Size, sep = "_"))

# 13h. Outliers -------------------------------------------------------------
# Outliers for each analysis
outliers <- data.frame(PersonID = accuracy$PersonID, Analysis = accuracy$Analysis)

# rank for MDS 125
outliers$MDS125 <- NA
tmp.pt <- nmds$CCIF125f$points["SCID", ]
tmp.tab <- sqrt((nmds$CCIF125f$points[,1] - tmp.pt[1])^2 + (nmds$CCIF125f$points[,2] - tmp.pt[2])^2)
outliers$MDS125[outliers$Analysis == "Slide"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))]))
outliers$MDS125[outliers$Analysis == "Digital"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))]))
rm(tmp.pt, tmp.tab)

# and mds 150
outliers$MDS150 <- NA
tmp.pt <- nmds$CCIF150f$points["SCID", ]
tmp.tab <- sqrt((nmds$CCIF150f$points[,1] - tmp.pt[1])^2 + (nmds$CCIF150f$points[,2] - tmp.pt[2])^2)
outliers$MDS150[outliers$Analysis == "Slide"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))]))
outliers$MDS150[outliers$Analysis == "Digital"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))]))
rm(tmp.pt, tmp.tab)

# for percentage accuracy
outliers$PtA125 <- NA
outliers$PtA125[outliers$Analysis == "Slide"][order(100-accuracySlide$PtA125[match(outliers$PersonID[outliers$Analysis == "Slide"], accuracySlide$PersonID)])] <- sort(rank(100-accuracySlide$PtA125))
outliers$PtA125[outliers$Analysis == "Digital"][order(100-accuracyDigital$PtA125[match(outliers$PersonID[outliers$Analysis == "Digital"], accuracyDigital$PersonID)])] <- sort(rank(100-accuracyDigital$PtA125))

outliers$PtA150 <- NA
outliers$PtA150[outliers$Analysis == "Slide"][order(100-accuracySlide$PtA150[match(outliers$PersonID[outliers$Analysis == "Slide"], accuracySlide$PersonID)])] <- sort(rank(100-accuracySlide$PtA150))
outliers$PtA150[outliers$Analysis == "Digital"][order(100-accuracyDigital$PtA150[match(outliers$PersonID[outliers$Analysis == "Digital"], accuracyDigital$PersonID)])] <- sort(rank(100-accuracyDigital$PtA150))

# SST
outliers$SST <- NA
tmp.pt <- abs(divTemp$SST10m[row.nam$s150] - 21.76)
names(tmp.pt) <- divTemp$Person[row.nam$s150]
outliers$SST[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))
tmp.pt <- abs(divTemp$SST10m[row.nam$d150] - 21.76)
names(tmp.pt) <- divTemp$Person[row.nam$d150]
outliers$SST[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# richness
outliers$IF_Richness125 <- NA
tmp.pt <- abs(divTemp$IF_Richness[row.nam$s125] - divTemp$IF_Richness[row.nam$s125c])
names(tmp.pt) <- divTemp$Person[row.nam$s125]
outliers$IF_Richness125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$IF_Richness[row.nam$d125] - divTemp$IF_Richness[row.nam$d125c])
names(tmp.pt) <- divTemp$Person[row.nam$d125]
outliers$IF_Richness125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$IF_Richness150 <- NA
tmp.pt <- abs(divTemp$IF_Richness[row.nam$s150] - divTemp$IF_Richness[row.nam$s150c])
names(tmp.pt) <- divTemp$Person[row.nam$s150]
outliers$IF_Richness150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$IF_Richness[row.nam$d150] - divTemp$IF_Richness[row.nam$d150c])
names(tmp.pt) <- divTemp$Person[row.nam$d150]
outliers$IF_Richness150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# Dominance
outliers$IF_Dominance125 <- NA
tmp.pt <- abs(divTemp$IF_Dominance[row.nam$s125] - divTemp$IF_Dominance[row.nam$s125c])
names(tmp.pt) <- divTemp$Person[row.nam$s125]
outliers$IF_Dominance125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$IF_Dominance[row.nam$d125] - divTemp$IF_Dominance[row.nam$d125c])
names(tmp.pt) <- divTemp$Person[row.nam$d125]
outliers$IF_Dominance125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$IF_Dominance150 <- NA
tmp.pt <- abs(divTemp$IF_Dominance[row.nam$s150] - divTemp$IF_Dominance[row.nam$s150c])
names(tmp.pt) <- divTemp$Person[row.nam$s150]
outliers$IF_Dominance150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$IF_Dominance[row.nam$d150] - divTemp$IF_Dominance[row.nam$d150c])
names(tmp.pt) <- divTemp$Person[row.nam$d150]
outliers$IF_Dominance150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# Evenness
outliers$IF_Evenness125 <- NA
tmp.pt <- abs(divTemp$IF_Evenness[row.nam$s125] - divTemp$IF_Evenness[row.nam$s125c])
names(tmp.pt) <- divTemp$Person[row.nam$s125]
outliers$IF_Evenness125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$IF_Evenness[row.nam$d125] - divTemp$IF_Evenness[row.nam$d125c])
names(tmp.pt) <- divTemp$Person[row.nam$d125]
outliers$IF_Evenness125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$IF_Evenness150 <- NA
tmp.pt <- abs(divTemp$IF_Evenness[row.nam$s150] - divTemp$IF_Evenness[row.nam$s150c])
names(tmp.pt) <- divTemp$Person[row.nam$s150]
outliers$IF_Evenness150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$IF_Evenness[row.nam$d150] - divTemp$IF_Evenness[row.nam$d150c])
names(tmp.pt) <- divTemp$Person[row.nam$d150]
outliers$IF_Evenness150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# ShannonWiener
outliers$IF_ShannonWiener125 <- NA
tmp.pt <- abs(divTemp$IF_ShannonWiener[row.nam$s125] - divTemp$IF_ShannonWiener[row.nam$s125c])
names(tmp.pt) <- divTemp$Person[row.nam$s125]
outliers$IF_ShannonWiener125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$IF_ShannonWiener[row.nam$d125] - divTemp$IF_ShannonWiener[row.nam$d125c])
names(tmp.pt) <- divTemp$Person[row.nam$d125]
outliers$IF_ShannonWiener125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$IF_ShannonWiener150 <- NA
tmp.pt <- abs(divTemp$IF_ShannonWiener[row.nam$s150] - divTemp$IF_ShannonWiener[row.nam$s150c])
names(tmp.pt) <- divTemp$Person[row.nam$s150]
outliers$IF_ShannonWiener150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$IF_ShannonWiener[row.nam$d150] - divTemp$IF_ShannonWiener[row.nam$d150c])
names(tmp.pt) <- divTemp$Person[row.nam$d150]
outliers$IF_ShannonWiener150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# weighted sums
outliers$fullSum <- rowSums(outliers[, 3:15])

outliers$wtSum <- rowSums(outliers[, grep("MDS", names(outliers))])/2 + rowSums(outliers[, grep("Pt", names(outliers))])/2 + outliers$SST + rowSums(outliers[, grep("Richness", names(outliers))])/8 + rowSums(outliers[, grep("Dominance", names(outliers))])/8 + rowSums(outliers[, grep("Evenness", names(outliers))])/8 + rowSums(outliers[, grep("ShannonWiener", names(outliers))])/8

outliers$wtSumPt <- outliers$wtSum / c(rep(4*17, 17), rep(4*9, 9)) * 100

# Table 8
# add percentage of specimens identified to accuracy
accuracySlide$ptID125 <- apply(slide125[,col.nam$s125], 2, function(x) sum(x != "na") / 300 * 100)[accuracySlide$PersonID]
accuracySlide$ptID150 <- apply(slide150[,col.nam$s150], 2, function(x) sum(x != "na") / 300 * 100)[accuracySlide$PersonID]
accuracyDigital$ptID125 <- apply(digital125[,col.nam$d125], 2, function(x) sum(x != "na") / 300 * 100)[accuracyDigital$PersonID]
accuracyDigital$ptID150 <- apply(digital150[,col.nam$d150], 2, function(x) sum(x != "na") / 300 * 100)[accuracyDigital$PersonID]

accuracySlide$Experience <- people$ExperienceSlideA[match(accuracySlide$PersonID, people$SlideID)]
accuracySlide$Experience[accuracySlide$PersonID %in% c("1a", "2a")] <- people$ExperienceSlideA[1:2]
accuracySlide$Experience[accuracySlide$PersonID %in% c("1b", "2b")] <- people$ExperienceSlideB[1:2]
accuracyDigital$Experience <- people$ExperienceDigital[match(accuracyDigital$PersonID, people$DigitalID)]

accuracy <- rbind(accuracySlide, accuracyDigital)
accuracy$Analysis <- c(rep("Slide", 17), rep("Digital", 9))

# Create a .csv file for the accuracy data

# create a subset of the data for this
tmp.sub <- accuracy[, c("PersonID", "Analysis", "Experience", "Routine", "PtA125", "PtA150", "ptID125", "ptID150")]
tmp.sub[, grep("125|150", names(tmp.sub))] <- round(tmp.sub[, grep("125|150", names(tmp.sub))], 2)
write.csv(tmp.sub, "Outputs/combconAccuracy.csv", row.names = FALSE)
rm(tmp.sub)

# add the person level info to this
tmp.s <- people[!is.na(people$SlideID), !grepl("Digital", names(people))]
tmp.d <- people[!is.na(people$DigitalID), !grepl("Slide", names(people))]
names(tmp.s)[1] <- names(tmp.d)[1] <- "PersonID"
tmp.s <- rbind(tmp.s, tmp.s[!is.na(tmp.s$ExperienceSlideB), ])
tmp.s$PersonID[tmp.s$PersonID %in% tmp.s$PersonID[duplicated(tmp.s$PersonID)]] <- paste(tmp.s$PersonID[tmp.s$PersonID %in% tmp.s$PersonID[duplicated(tmp.s$PersonID)]], rep(c("a", "b"), each = 2), sep = "")
tmp.s$ExperienceSlideA[grep("b", tmp.s$PersonID)] <- tmp.s$ExperienceSlideB[grep("b", tmp.s$PersonID)]
tmp.s <- tmp.s[, names(tmp.s) != "ExperienceSlideB"]
names(tmp.d) <- names(tmp.s) <- gsub("SlideA", "", names(tmp.s))
tmp <- rbind(tmp.s, tmp.d)
tmp.out <- merge(outliers, tmp)
tmp.out <- merge(tmp.out, accuracy[, c("PersonID", "ptID125", "ptID150")])

#par(ask = TRUE)
par(mfrow = c(4, 4))
par(mar = c(2,2,1,1))
for (i in names(tmp.out)[c(19:23, 25:26)]) {
  for (j in names(tmp.out)[3:18])
    plot(factor(tmp.out[tmp.out$Analysis == "Slide",i]), tmp.out[tmp.out$Analysis == "Slide",j], main = i, col = "red")
  for (j in names(tmp.out)[3:18])
    plot(factor(tmp.out[tmp.out$Analysis == "Digital",i]), tmp.out[tmp.out$Analysis == "Digital",j], main = j, col = "blue")
}
rm(i, j)
#par(ask = FALSE)
par(mfrow = c(1,1))
par(mar = c(5.1, 4.1, 4.1, 2.1))
rm(tmp.s, tmp.out, tmp.d, tmp)


# outliers for digital and slide combined
# rank for MDS 125
outliers$cMDS125 <- NA
tmp.pt <- nmds$CCIF125f$points["SCID", ]
tmp.tab <- sqrt((nmds$CCIF125f$points[,1] - tmp.pt[1])^2 + (nmds$CCIF125f$points[,2] - tmp.pt[2])^2)
outliers$cMDS125[order(tmp.tab[match(outliers$PersonID, names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID, names(tmp.tab))]))
rm(tmp.pt, tmp.tab)

# and mds 150
outliers$cMDS150 <- NA
tmp.pt <- nmds$CCIF150f$points["SCID", ]
tmp.tab <- sqrt((nmds$CCIF150f$points[,1] - tmp.pt[1])^2 + (nmds$CCIF150f$points[,2] - tmp.pt[2])^2)
outliers$cMDS150[order(tmp.tab[match(outliers$PersonID, names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID, names(tmp.tab))]))
rm(tmp.pt, tmp.tab)

# for percentage accuracy
outliers$cPtA125 <- NA
outliers$cPtA125[order(100-accuracy$PtA125[match(outliers$PersonID, accuracy$PersonID)])] <- sort(rank(100-accuracy$PtA125))

outliers$cPtA150 <- NA
outliers$cPtA150[order(100-accuracy$PtA150[match(outliers$PersonID, accuracy$PersonID)])] <- sort(rank(100-accuracy$PtA150))

# SST
outliers$cSST <- NA
tmp.pt <- abs(divTemp$SST10m[c(row.nam$s150, row.nam$d150)] - 21.76)
names(tmp.pt) <- divTemp$Person[c(row.nam$s150, row.nam$d150)]
outliers$cSST[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# richness
outliers$cIF_Richness125 <- NA
tmp.pt <- abs(divTemp$IF_Richness[c(row.nam$s125, row.nam$d125)] - divTemp$IF_Richness[row.nam$s125c])
names(tmp.pt) <- divTemp$Person[c(row.nam$s125, row.nam$d125)]
outliers$cIF_Richness125[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$cIF_Richness150 <- NA
tmp.pt <- abs(divTemp$IF_Richness[c(row.nam$s150, row.nam$d150)] - divTemp$IF_Richness[row.nam$s150c])
names(tmp.pt) <- divTemp$Person[c(row.nam$s150, row.nam$d150)]
outliers$cIF_Richness150[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# Dominance
outliers$cIF_Dominance125 <- NA
tmp.pt <- abs(divTemp$IF_Dominance[c(row.nam$s125, row.nam$d125)] - divTemp$IF_Dominance[row.nam$s125c])
names(tmp.pt) <- divTemp$Person[c(row.nam$s125, row.nam$d125)]
outliers$cIF_Dominance125[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$cIF_Dominance150 <- NA
tmp.pt <- abs(divTemp$IF_Dominance[c(row.nam$s150, row.nam$d150)] - divTemp$IF_Dominance[row.nam$s150c])
names(tmp.pt) <- divTemp$Person[c(row.nam$s150, row.nam$d150)]
outliers$cIF_Dominance150[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# Evenness
outliers$cIF_Evenness125 <- NA
tmp.pt <- abs(divTemp$IF_Evenness[c(row.nam$s125, row.nam$d125)] - divTemp$IF_Evenness[row.nam$s125c])
names(tmp.pt) <- divTemp$Person[c(row.nam$s125, row.nam$d125)]
outliers$cIF_Evenness125[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$cIF_Evenness150 <- NA
tmp.pt <- abs(divTemp$IF_Evenness[c(row.nam$s150, row.nam$d150)] - divTemp$IF_Evenness[row.nam$s150c])
names(tmp.pt) <- divTemp$Person[c(row.nam$s150, row.nam$d150)]
outliers$cIF_Evenness150[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# ShannonWiener
outliers$cIF_ShannonWiener125 <- NA
tmp.pt <- abs(divTemp$IF_ShannonWiener[c(row.nam$s125, row.nam$d125)] - divTemp$IF_ShannonWiener[row.nam$s125c])
names(tmp.pt) <- divTemp$Person[c(row.nam$s125, row.nam$d125)]
outliers$cIF_ShannonWiener125[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$cIF_ShannonWiener150 <- NA
tmp.pt <- abs(divTemp$IF_ShannonWiener[c(row.nam$s150, row.nam$d150)] - divTemp$IF_ShannonWiener[row.nam$s150c])
names(tmp.pt) <- divTemp$Person[c(row.nam$s150, row.nam$d150)]
outliers$cIF_ShannonWiener150[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# weighted sums
outliers$cfullSum <- rowSums(outliers[, grep("^c", names(outliers))])

outliers$cwtSum <- rowSums(outliers[, grep("cMDS", names(outliers))])/2 + rowSums(outliers[, grep("^c.*Pt", names(outliers))])/2 + outliers$cSST + rowSums(outliers[, grep("cIF_Richness", names(outliers))])/8 + rowSums(outliers[, grep("cIF_Dominance", names(outliers))])/8 + rowSums(outliers[, grep("cIF_Evenness", names(outliers))])/8 + rowSums(outliers[, grep("cIF_ShannonWiener", names(outliers))])/8

outliers$cwtSumPt <- outliers$cwtSum / 26 * 100

# looking for correlations
pairs(outliers[, grep("^c", names(outliers))])
pairs(outliers[, grep("^c.*125", names(outliers))])
pairs(outliers[, grep("^c.*150", names(outliers))])

# add the person level info to this
tmp.s <- people[!is.na(people$SlideID), !grepl("Digital", names(people))]
tmp.d <- people[!is.na(people$DigitalID), !grepl("Slide", names(people))]
names(tmp.s)[1] <- names(tmp.d)[1] <- "PersonID"
tmp.s <- rbind(tmp.s, tmp.s[!is.na(tmp.s$ExperienceSlideB), ])
tmp.s$PersonID[tmp.s$PersonID %in% tmp.s$PersonID[duplicated(tmp.s$PersonID)]] <- paste(tmp.s$PersonID[tmp.s$PersonID %in% tmp.s$PersonID[duplicated(tmp.s$PersonID)]], rep(c("a", "b"), each = 2), sep = "")
tmp.s$ExperienceSlideA[grep("b", tmp.s$PersonID)] <- tmp.s$ExperienceSlideB[grep("b", tmp.s$PersonID)]
tmp.s <- tmp.s[, names(tmp.s) != "ExperienceSlideB"]
names(tmp.d) <- names(tmp.s) <- gsub("SlideA", "", names(tmp.s))
tmp <- rbind(tmp.s, tmp.d)
tmp.out <- merge(outliers, tmp)
tmp.out <- merge(tmp.out, accuracy[, c("PersonID", "ptID125", "ptID150")])

#par(ask = TRUE)
par(mfrow = c(4, 4))
par(mar = c(2,2,1,1))
for (i in names(tmp.out)[c(35:39, 41:42)]) {
  for (j in grep("^c", names(tmp.out), value = TRUE))
    plot(factor(tmp.out[,i]), tmp.out[,j], main = j, col = "red")
}
rm(i, j)
#par(ask = FALSE)
par(mfrow = c(1,1))
par(mar = c(5.1, 4.1, 4.1, 2.1))

plot(tmp.out$Analysis, tmp.out$cwtSumPt)

rm(tmp.s, tmp.out, tmp.d, tmp)



# 13i. Size vs. maximum agreement ------------------------------------------
# Figure 5
png("ASFigures/CombCon/Fig5_size_agreement_sd125.png")
with(size125, plot(slideAgreement, Length, pch = 16, main = "> 125", las = 1, xlab = "Agreement", ylab = expression(paste("Maximum diamter / ", mu, "m"))))
lines(names(tapply(size125$Length, size125$slideAgreement, max)), tapply(size125$Length, size125$slideAgreement, max), pch = 16)
with(size125, points(digitalAgreement, Length, pch = 16, col = "blue"))
lines(names(tapply(size125$Length, size125$digitalAgreement, max)), tapply(size125$Length, size125$digitalAgreement, max), pch = 16, col = 4)
legend("topleft", col = c(1, 4), pch = 16, legend = c("Slide", "Digital"))
dev.off()

png("ASFigures/CombCon/Fig5_size_agreement_sd150.png")
with(size150, plot(slideAgreement, Length, pch = 16, main = "> 150", las = 1, xlab = "Agreement", ylab = expression(paste("Maximum diamter / ", mu, "m"))))
lines(names(tapply(size150$Length, size150$slideAgreement, max)), tapply(size150$Length, size150$slideAgreement, max), pch = 16)
with(size150, points(digitalAgreement, Length, pch = 16, col = "blue"))
lines(names(tapply(size150$Length, size150$digitalAgreement, max)), tapply(size150$Length, size150$digitalAgreement, max), pch = 16, col = 4)
legend("topleft", col = c(1, 4), pch = 16, legend = c("Slide", "Digital"))
dev.off()

png("ASFigures/CombCon/Fig5_size_agreement_125.png")
with(size125, plot(Length, Agreement, pch = 16, main = "> 125", las = 1, ylab = "Agreement", xlab = expression(paste("Maximum diamter / ", mu, "m"))))
lines(tapply(size125$Length, size125$Agreement, max), names(tapply(size125$Length, size125$Agreement, max)), pch = 16)
dev.off()

png("ASFigures/CombCon/Fig5_size_agreement_150.png")
with(size150, plot(Length, Agreement, pch = 16, main = "> 150", las = 1, ylab = "Agreement", xlab = expression(paste("Maximum diamter / ", mu, "m"))))
lines(tapply(size150$Length, size150$Agreement, max), names(tapply(size150$Length, size150$Agreement, max)), pch = 16)
dev.off()

# 13j. Comparison with ForCenS -----------------------------------------------
# see how the datasets compare with the ForCenS values for that latitude
ForCenSred <- as.data.frame(read_excel("Data/Siccha_ForCenS.xlsx", sheet = "ForCenSred", na = "N/A"))
ForCenSred[, 22:62][is.na(ForCenSred[, 22:62])] <- 0
str(ForCenSred)

png("ASFigures/CombCon/Div_cf_ForCenSred.png", 550, 600)
par(mfrow = c(2,2))
# species richness
tmp.rich <- specnumber(ForCenSred[, 22:62]) # richness(ForCenSred[,22:62])
plot(ForCenSred$Latitude, tmp.rich, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Richness")
points(rep(30.2, nrow(divTemp)), divTemp$IF_Richness, col = "blue", pch = 16)
points(rep(30.2, 4), divTemp$IF_Richness[divTemp$Person == "consensus"], col = "red", pch = 16)

# ShannonWiener
tmp.sw <- diversity(ForCenSred[,22:62])
plot(ForCenSred$Latitude, tmp.sw, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Shannon Wiener")
points(rep(30.2, nrow(divTemp)), divTemp$IF_ShannonWiener, col = "blue", pch = 16)
points(rep(30.2, 4), divTemp$IF_ShannonWiener[divTemp$Person == "consensus"], col = "red", pch = 16)

# Dominance
tmp.dom <- (1 - diversity(ForCenSred[,22:62], index = "simpson"))
plot(ForCenSred$Latitude, tmp.dom, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Dominance")
points(rep(30.2, nrow(divTemp)), divTemp$IF_Dominance, col = "blue", pch = 16)
points(rep(30.2, 4), divTemp$IF_Dominance[divTemp$Person == "consensus"], col = "red", pch = 16)


# Evenness
tmp.eve <- (exp(diversity(ForCenSred[,22:62])) / specnumber(ForCenSred[,22:62]))
plot(ForCenSred$Latitude, tmp.eve, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Evenness")
points(rep(30.2, nrow(divTemp)), divTemp$IF_Evenness, col = "blue", pch = 16)
points(rep(30.2, 4), divTemp$IF_Evenness[divTemp$Person == "consensus"], col = "red", pch = 16)
par(mfrow = c(1,1))
dev.off()

png("ASFigures/CombCon/Div_cf_ForCenSred_Atl.png", 550, 600)
par(mfrow = c(2,2))
# species richness
with(ForCenSred[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11, ], plot(Latitude, tmp.rich[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], pch = 16, col = "grey50", xlab = "Latitude", ylab = "Richness"))
points(rep(30.2, nrow(divTemp)), divTemp$IF_Richness, col = "blue", pch = 16)
points(rep(30.2, 4), divTemp$IF_Richness[divTemp$Person == "consensus"], col = "red", pch = 16)


# ShannonWiener
plot(ForCenSred$Latitude[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], tmp.sw[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], pch = 16, col = "grey50", xlab = "Latitude", ylab = "Shannon Wiener")
points(rep(30.2, nrow(divTemp)), divTemp$IF_ShannonWiener, col = "blue", pch = 16)
points(rep(30.2, 4), divTemp$IF_ShannonWiener[divTemp$Person == "consensus"], col = "red", pch = 16)

# Dominance
plot(ForCenSred$Latitude[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], tmp.dom[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], pch = 16, col = "grey50", xlab = "Latitude", ylab = "Dominance")
points(rep(30.2, nrow(divTemp)), divTemp$IF_Dominance, col = "blue", pch = 16)
points(rep(30.2, 4), divTemp$IF_Dominance[divTemp$Person == "consensus"], col = "red", pch = 16)


# Evenness
plot(ForCenSred$Latitude[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], tmp.eve[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], pch = 16, col = "grey50", xlab = "Latitude", ylab = "Evenness")
points(rep(30.2, nrow(divTemp)), divTemp$IF_Evenness, col = "blue", pch = 16)
points(rep(30.2, 4), divTemp$IF_Evenness[divTemp$Person == "consensus"], col = "red", pch = 16)
par(mfrow = c(1,1))
dev.off()

rm(tmp.sw, tmp.dom, tmp.eve, tmp.rich, ForCenSred)


# 14. Save the data -------------------------------------------------------
save.image("Outputs/Reanalysis_NA_IF.RData")
# output the full species names
tmp.125<- full.125
tmp.150 <- full.150
for(i in 1:nrow(sp.abb)) {
  tmp.125 <- as.data.frame(apply(tmp.125, 2, function(x) gsub(paste("^", sp.abb$Abbreviation[i], "$", sep = ""), sp.abb$Species[i], x)))
  tmp.150 <- as.data.frame(apply(tmp.150, 2, function(x) gsub(paste("^", sp.abb$Abbreviation[i], "$", sep = ""), sp.abb$Species[i], x)))
}
write.csv(tmp.125, "Outputs/PersonIDs_125_IF.csv", row.names = FALSE)
write.csv(tmp.150, "Outputs/PersonIDs_150_IF.csv", row.names = FALSE)
rm(tmp.125, tmp.150, i)
