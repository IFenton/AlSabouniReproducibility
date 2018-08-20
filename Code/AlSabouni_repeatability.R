# Code for the repeatability analysis 
# Al-Sabouni, Fenton, Telford & Kucera
# Isabel Fenton
# Project: IF reanalysis
# Date created: 5 / 3 / 2018
# Date last edited: 17 / 8 / 2018
# 
# All the code needed to run the Al-Sabouni et al repeatability paper. Lab book repeatability
# 
# Previous file: Reanalysis_NA_IF.R
# Next file:
# 

rm(list = ls())

# Inputs ------------------------------------------------------------------
# these files are assumed to be in a folder called Data
# AlSabouni_PersonIDs.xlsx - the raw data and the species abbreviations
# AlSabouni_PeopleMetadata.xlsx - the person metadata for the participants
# AlSabouni_SpecimenSize.xlsx - the size data for the specimens
# AlSabouni_DiversityTemp.xlsx - the temperature data for the sites
# Siccha_ForCenS.xlsx - the ForCenS data for comparison of diversity ranges

# Outputs -----------------------------------------------------------------
# The assumed file struture is:
# Figures are written to:
# /ASFigures, and also /ASFigures/DigitalSlide and /ASFigures/Time
# Other outputs are written to:
# /ASOutputs


# Source files / libraries ------------------------------------------------
library(readxl) # reading in xlsx files
library(caret) # for the confusion matrix
library(colorRamps) # colours
library(RColorBrewer) # colours
library(viridis) # colours for confusion matrix
library(stringr) # for confusion matrix axes
library(vegan) # for the diversity metrics
library(cluster) # for the NMDS analysis

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
sp.abb <- as.data.frame(read_excel("Data/AlSabouni_PersonIDs.xlsx", sheet = "Abbreviations", na = "NA"))

# people metadata
people <- as.data.frame(read_excel("Data/AlSabouni_PeopleMetadata.xlsx", na = "NA"))

# specimen size
size125 <- as.data.frame(read_excel("Data/AlSabouni_SpecimenSize.xlsx", sheet = "Size125"))
size150 <- as.data.frame(read_excel("Data/AlSabouni_SpecimenSize.xlsx", sheet = "Size150"))

# diversity / temperature 
divTemp <- as.data.frame(read_excel("Data/AlSabouni_DiversityTemp.xlsx", na = "NA", sheet = "SST"))

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

# create plots for slide / digital separately
# figures showing the cumulative maximum consensus for the IDs
cum.sum <- list()

cum.sum$s125 <- cumsum(rev(table(full.125$sMaxCon)))/300*100
cum.sum$s150 <- cumsum(rev(table(full.150$sMaxCon)))/300*100
cum.sum$d125 <- cumsum(rev(table(full.125$dMaxCon)))/300*100
cum.sum$d150 <- cumsum(rev(table(full.150$dMaxCon)))/300*100

png("ASFigures/Cumulative_slide.png")
plot(names(cum.sum$s125), cum.sum$s125, type = "n", xlab = "Number of participants", ylab = "Fraction of specimens / %", las = 1, ylim = c(0, 100), lty = 2, main = "Slide", xaxt = "n", cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.1)
axis(1, 3:17)
legend("topright", lty = c(2, 1), legend = c(expression(paste(">125 ", mu, "m")), expression(paste(">150 ", mu, "m"))), cex = 1.3, lwd = 2)

# adding in the consensus lines
tmp.s150 <- rev(cum.sum$s150)[as.numeric(names(rev(cum.sum$s150))) > c50_cutoff$slide][1]
lines(c(names(tmp.s150), names(tmp.s150)), c(-5, tmp.s150), col = "grey")
lines(c(2, as.numeric(names(tmp.s150))), c(tmp.s150, tmp.s150), col = "grey")

tmp.s125 <- rev(cum.sum$s125)[as.numeric(names(rev(cum.sum$s125))) > c50_cutoff$slide][1]
lines(c(names(tmp.s125), names(tmp.s125)), c(-5, tmp.s125), col = "grey", lty = 2)
lines(c(2, names(tmp.s125)), c(tmp.s125, tmp.s125), col = "grey", lty = 2)

# plotting the cumulative curves
points(names(cum.sum$s125), cum.sum$s125, type = "s", lty = 2)
points(names(cum.sum$s150), cum.sum$s150, type = "s")
dev.off()

png("ASFigures/Cumulative_digital.png", 380, 480)
plot(names(cum.sum$d125), cum.sum$d125, type = "n", xlab = "Number of participants", ylab = "Fraction of specimens / %", las = 1, ylim = c(0, 100), lty = 2, main = "Digital", xaxt = "n", cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.1)
axis(1, 2:9)
legend("topright", lty = c(2, 1), legend = c(expression(paste(">125 ", mu, "m")), expression(paste(">150 ", mu, "m"))), cex = 1.3, lwd = 2, col = "blue")

# adding in the consensus lines
tmp.d150 <- rev(cum.sum$d150)[as.numeric(names(rev(cum.sum$d150))) > c50_cutoff$digital][1]
lines(c(names(tmp.d150), names(tmp.d150)), c(-5, tmp.d150), col = "grey")
lines(c(1, names(tmp.d150)), c(tmp.d150, tmp.d150), col = "grey")

tmp.d125 <- rev(cum.sum$d125)[as.numeric(names(rev(cum.sum$d125))) > c50_cutoff$digital][1]
lines(c(names(tmp.d125), names(tmp.d125)), c(-5, tmp.d125), col = "grey", lty = 2)
lines(c(1, names(tmp.d125)), c(tmp.d125, tmp.d125), col = "grey", lty = 2)

# plotting the cumulative curves
points(names(cum.sum$d125), cum.sum$d125, type = "s", lty = 2, col = "blue")
points(names(cum.sum$d150), cum.sum$d150, type = "s", col = "blue")
par(mfrow = c(1,1))
dev.off()

rm(tmp.s125, tmp.s150, tmp.d125, tmp.d150)

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
text(25, 90, expression(paste(">125 ", mu, "m")), cex = 1.3)
text(1, 40, "Slide", cex = 1.3)
text(21.75, 40, "Digital", cex = 1.3, col = "blue")

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
text(25, 90, expression(paste(">150 ", mu, "m")), cex = 1.3)
text(1, 40, "Slide", cex = 1.3)
text(21.75, 40, "Digital", cex = 1.3, col = "blue")

par(mfrow = c(1,1))
dev.off()

# Compare digital / slide accuracy to the separate consensus' (n.b. order is approx. experience)
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
text(25, 90, expression(paste(">125 ", mu, "m")), cex = 1.3)
text(1, 40, "Slide", cex = 1.3)
text(21.75, 40, "Digital", cex = 1.3, col = "blue")


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
text(25, 90, expression(paste(">150 ", mu, "m")), cex = 1.3)
text(1, 40, "Slide", cex = 1.3)
text(21.75, 40, "Digital", cex = 1.3, col = "blue")


par(mfrow = c(1,1))
dev.off()



# 4. Sensitivity to alphabetical order -----------------------------------
# Does the way in which the consensus split is chosen affect the results?
# try taking the last choice rather than the first

# 4a. Recalculate the Consensus ID ----------------------------------------
# for the combined consensus
full.125$cCIDr <- apply(full.125[, col.nam$c125], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1])
full.150$cCIDr <- apply(full.150[, col.nam$c150], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1])
# and split by slide / digital
full.125$sCIDr <- apply(full.125[, col.nam$s125], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1])
full.150$sCIDr <- apply(full.150[, col.nam$s150], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1])
full.125$dCIDr <- apply(full.125[, col.nam$d125], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1])
full.150$dCIDr <- apply(full.150[, col.nam$d150], 1, function (x) rev(names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))])[1])

# 4b. Consider the effect on agreement -------------------------------------
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

# 4c. Influence on the agreement plots --------------------------------------------------
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
text(25, 90, expression(paste(">125 ", mu, "m")), cex = 1.3)
text(1, 40, "Slide", cex = 1.3)
text(21.75, 40, "Digital", cex = 1.3, col = "blue")


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
text(25, 90, expression(paste(">150 ", mu, "m")), cex = 1.3)
text(1, 40, "Slide", cex = 1.3)
text(21.75, 40, "Digital", cex = 1.3, col = "blue")


par(mfrow = c(1,1))
dev.off()


# 4d. Excluding points with no consensus ----------------------------------
# for the combined consensus
full.125$cCIDnc <- apply(full.125[, col.nam$c125], 1, function (x) ifelse(sum(table(x[x != 'na']) == max(table(x[x != 'na']))) == 1, names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))], "nca"))
full.150$cCIDnc <- apply(full.150[, col.nam$c150], 1, function (x) ifelse(sum(table(x[x != 'na']) == max(table(x[x != 'na']))) == 1, names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))], "nca"))
# and split by slide / digital
full.125$sCIDnc <- apply(full.125[, col.nam$s125], 1, function (x) ifelse(sum(table(x[x != 'na']) == max(table(x[x != 'na']))) == 1, names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))], "nca"))
full.150$sCIDnc <- apply(full.150[, col.nam$s150], 1, function (x) ifelse(sum(table(x[x != 'na']) == max(table(x[x != 'na']))) == 1, names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))], "nca"))
full.125$dCIDnc <- apply(full.125[, col.nam$d125], 1, function (x) ifelse(sum(table(x[x != 'na']) == max(table(x[x != 'na']))) == 1, names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))], "nca"))
full.150$dCIDnc <- apply(full.150[, col.nam$d150], 1, function (x) ifelse(sum(table(x[x != 'na']) == max(table(x[x != 'na']))) == 1, names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))], "nca"))

# based on the consensus for the combined
accuracy$cPtA125nc <- apply(full.125[,col.nam$c125], 2, function(x) sum(x == full.125$cCIDnc) / 300 * 100)
accuracy$cPtA150nc <- apply(full.150[,col.nam$c150], 2, function(x) sum(x == full.150$cCIDnc) / 300 * 100)

# then split by slide / digital
accuracy$dPtA150nc <- accuracy$dPtA125nc <- accuracy$sPtA150nc <- accuracy$sPtA125nc <- NA 
accuracy$sPtA125nc[accuracy$Analysis == "Slide"] <- apply(full.125[,col.nam$s125], 2, function(x) sum(x == full.125$sCIDnc) / 300 * 100)
accuracy$sPtA150nc[accuracy$Analysis == "Slide"] <- apply(full.150[,col.nam$s150], 2, function(x) sum(x == full.150$sCIDnc) / 300 * 100)
accuracy$dPtA125nc[accuracy$Analysis == "Digital"] <- apply(full.125[,col.nam$d125], 2, function(x) sum(x == full.125$dCIDnc) / 300 * 100)
accuracy$dPtA150nc[accuracy$Analysis == "Digital"] <- apply(full.150[,col.nam$d150], 2, function(x) sum(x == full.150$dCIDnc) / 300 * 100)

# calculate the mean using the combined consensus
CID_mn$csPtA125nc <- mean(accuracy$cPtA125nc[accuracy$Analysis == "Slide"])
CID_mn$csPtA150nc <- mean(accuracy$cPtA150nc[accuracy$Analysis == "Slide"])
CID_mn$cdPtA125nc <- mean(accuracy$cPtA125nc[accuracy$Analysis == "Digital"])
CID_mn$cdPtA150nc <- mean(accuracy$cPtA150nc[accuracy$Analysis == "Digital"])
# and with the separate consensus'
CID_mn$ssPtA125nc <- mean(accuracy$sPtA125nc, na.rm = TRUE)
CID_mn$ssPtA150nc <- mean(accuracy$sPtA150nc, na.rm = TRUE)
CID_mn$ddPtA125nc <- mean(accuracy$dPtA125nc, na.rm = TRUE)
CID_mn$ddPtA150nc <- mean(accuracy$dPtA150nc, na.rm = TRUE)

# calculate the sd percentage agreement for each of the four analyses, using initially the combined consensus
CID_sd$csPtA125nc <- sd(accuracy$cPtA125nc[accuracy$Analysis == "Slide"])
CID_sd$csPtA150nc <- sd(accuracy$cPtA150nc[accuracy$Analysis == "Slide"])
CID_sd$cdPtA125nc <- sd(accuracy$cPtA125nc[accuracy$Analysis == "Digital"])
CID_sd$cdPtA150nc <- sd(accuracy$cPtA150nc[accuracy$Analysis == "Digital"])
# and with the separate consensus'
CID_sd$ssPtA125nc <- sd(accuracy$sPtA125nc, na.rm = TRUE)
CID_sd$ssPtA150nc <- sd(accuracy$sPtA150nc, na.rm = TRUE)
CID_sd$ddPtA125nc <- sd(accuracy$dPtA125nc, na.rm = TRUE)
CID_sd$ddPtA150nc <- sd(accuracy$dPtA150nc, na.rm = TRUE)

# Influence on the agreement plots 
png("ASFigures/CombCon_agreement_fullID_nc.png", 800, 1000)
par(mfrow = c(2, 1))
with(accuracy, plot(cPtA125nc[match(ord.div, PersonID)], ylim = c(40, 90), type = "n", xaxt = "n", ylab = "Percentage agreement", xlab = "Participant", cex.lab = 1.5, las = 2, cex.axis = 1.1))
axis(1, at = 1:26, labels = accuracy$PersonID[match(ord.div, accuracy$PersonID)], cex.axis = 1.1)
# adding mean values
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA125nc, 2), lty = 1)
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA125nc - CID_sd$csPtA125nc, 2), lty = 4)
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA125nc + CID_sd$csPtA125nc, 2), lty = 4)
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA125nc, 2), lty = 1, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA125nc - CID_sd$cdPtA125nc, 2), lty = 4, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA125nc + CID_sd$cdPtA125nc, 2), lty = 4, col = "blue")
abline(v = 20.5, col = "grey 50")
# points
with(accuracy, points(1:26, cPtA125nc[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
with(accuracy, arrows(1, cPtA125nc[PersonID == "1a"], 2, cPtA125nc[PersonID == "1b"], length = 0.14))
with(accuracy, arrows(3, cPtA125nc[PersonID == "2a"], 4, cPtA125nc[PersonID == "2b"], length = 0.14))
text(25, 90, expression(paste(">125 ", mu, "m")), cex = 1.3)
text(1, 40, "Slide", cex = 1.3)
text(21.75, 40, "Digital", cex = 1.3, col = "blue")


with(accuracy, plot(cPtA150nc[match(ord.div, PersonID)], ylim = c(40, 90), type = "n", xaxt = "n", ylab = "Percentage agreement", xlab = "Participant", cex.lab = 1.5, las = 2, cex.axis = 1.1))
axis(1, at = 1:26, labels = accuracy$PersonID[match(ord.div, accuracy$PersonID)], cex.axis = 1.1)
# adding mean values
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA150nc, 2), lty = 1)
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA150nc - CID_sd$csPtA150nc, 2), lty = 4)
lines(x = c(0, 20.5), y = rep(CID_mn$csPtA150nc + CID_sd$csPtA150nc, 2), lty = 4)
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA150nc, 2), lty = 1, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA150nc - CID_sd$cdPtA150nc, 2), lty = 4, col = "blue")
lines(x = c(20.5, 27), y = rep(CID_mn$cdPtA150nc + CID_sd$cdPtA150nc, 2), lty = 4, col = "blue")
abline(v = 20.5, col = "grey 50")
# points
with(accuracy, points(1:26, cPtA150nc[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
with(accuracy, arrows(1, cPtA150nc[PersonID == "1a"], 2, cPtA150nc[PersonID == "1b"], length = 0.14))
with(accuracy, arrows(3, cPtA150nc[PersonID == "2a"], 4, cPtA150nc[PersonID == "2b"], length = 0.14))
text(25, 90, expression(paste(">150 ", mu, "m")), cex = 1.3)
text(1, 40, "Slide", cex = 1.3)
text(21.75, 40, "Digital", cex = 1.3, col = "blue")


par(mfrow = c(1,1))
dev.off()


# 5. Confusion matrix ----------------------------------------------------

# 5a. Confusion matrices for the main results -----------------------------
# initially with slide 125
# this requires the data to be in long format
long <- list()
long$s125 <- reshape(full.125[, -col.nam$d125], varying = list(names(full.125)[col.nam$s125]), direction = "long", times = names(full.125)[col.nam$s125], timevar = "Person")
rownames(long$s125) <- 1:nrow(long$s125)
long$s125 <- long$s125[, (names(long$s125) != "id")]
names(long$s125)[names(long$s125) == "1a"] <- "origID"
head(long$s125)
tail(long$s125)

conf.mat <- list()

# plot the confusion matrix
png("ASFigures/CombCon_conf_slide125.png", 1000, 700)
conf_mat(long$s125, "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext(expression(paste("Slide >125", mu, "m")), 1, cex = 2, adj = -1.2,  line = -6)
dev.off()

# work out the order
sp.abb$confOrder <- NA
sp.abb$confOrder[sp.abb$Abbreviation %in% c("nc", "na", "nca")] <- nrow(sp.abb) - 0:2
sp.abb$confOrder[is.na(sp.abb$confOrder)] <- 1:(nrow(sp.abb) - 3)
# calculate the confusion matrix
conf.mat$s125 <- confusionMatrix(factor(sp.abb$Species[match(long$s125$cCID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s125$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# same for the other datasets
# so slide 150
long$s150 <- reshape(full.150[, -col.nam$d150], varying = list(names(full.150)[col.nam$s150]), direction = "long", times = names(full.150)[col.nam$s150], timevar = "Person")
rownames(long$s150) <- 1:nrow(long$s150)
long$s150 <- long$s150[, (names(long$s150) != "id")]
names(long$s150)[names(long$s150) == "1a"] <- "origID"
head(long$s150)
tail(long$s150)

png("ASFigures/CombCon_conf_slide150.png", 1000, 610)
conf_mat(long$s150, "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext(expression(paste("Slide >150", mu, "m")), 1, cex = 2, adj = -1.2,  line = -6)
dev.off()
conf.mat$s150 <- confusionMatrix(factor(sp.abb$Species[match(long$s150$cCID, sp.abb$Abbreviation)], levels = sp.abb$Species), factor(sp.abb$Species[match(long$s150$origID, sp.abb$Abbreviation)], levels = sp.abb$Species))

# calculate the confusion matrix
conf.mat$s150 <- confusionMatrix(factor(sp.abb$Species[match(long$s150$cCID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s150$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# digital 125
long$d125 <- reshape(full.125[, -col.nam$s125], varying = list(names(full.125)[col.nam$d125]), direction = "long", times = names(full.125)[col.nam$d125], timevar = "Person")
rownames(long$d125) <- 1:nrow(long$d125)
long$d125 <- long$d125[, (names(long$d125) != "id")]
names(long$d125)[names(long$d125) == "A"] <- "origID"
head(long$d125)
tail(long$d125)

png("ASFigures/CombCon_conf_digital125.png", 1000, 700)
conf_mat(long$d125, "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext(expression(paste("Digital >125", mu, "m")), 1, cex = 2, adj = -1.2,  line = -6)
dev.off()

# calculate the confusion matrix
conf.mat$d125 <- confusionMatrix(factor(sp.abb$Species[match(long$d125$cCID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$d125$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# digital 150
long$d150 <- reshape(full.150[, -col.nam$s125], varying = list(names(full.150)[col.nam$d150]), direction = "long", times = names(full.150)[col.nam$d150], timevar = "Person")
rownames(long$d150) <- 1:nrow(long$d150)
long$d150 <- long$d150[, (names(long$d150) != "id")]
names(long$d150)[names(long$d150) == "A"] <- "origID"
head(long$d150)
tail(long$d150)

png("ASFigures/CombCon_conf_digital150.png", 1000, 610)
conf_mat(long$d150, "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext(expression(paste("Digital >150", mu, "m")), 1, cex = 2, adj = -1.2,  line = -6)
dev.off()

# calculate the confusion matrix
conf.mat$d150 <- confusionMatrix(factor(sp.abb$Species[match(long$d150$cCID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$d150$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# What about for the separate consensus values
png("ASFigures/SepCon_conf_slide125.png", 1000, 700)
conf_mat(long$s125, "origID", "sCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Sep Slide 125", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

# calculate the confusion matrix
conf.mat$s125SC <- confusionMatrix(factor(sp.abb$Species[match(long$s125$sCID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s125$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/SepCon_conf_slide150.png", 1000, 610)
conf_mat(long$s150, "origID", "sCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Sep Slide 150", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

# calculate the confusion matrix
conf.mat$s150SC <- confusionMatrix(factor(sp.abb$Species[match(long$s150$sCID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s150$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/SepCon_conf_digital125.png", 1000, 700)
conf_mat(long$d125, "origID", "dCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Sep Digital 125", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

# calculate the confusion matrix
conf.mat$d125SC <- confusionMatrix(factor(sp.abb$Species[match(long$d125$dCID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$d125$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/SepCon_conf_digital150.png", 1000, 610)
conf_mat(long$d150, "origID", "dCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Sep Digital 150", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

conf.mat$d125SC <- confusionMatrix(factor(sp.abb$Species[match(long$d125$sCID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$d125$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# 5b. Influence of reversed order on the confusion matrices ---------------------------------
# influence of the alphabetical ordering
png("ASFigures/CombCon_conf_slide125r.png", 1000, 700)
conf_mat(long$s125, "origID", "cCIDr", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Slide 125", 1, cex = 2, adj = -0.8,  line = -5)
dev.off()

conf.mat$s125r <- confusionMatrix(factor(sp.abb$Species[match(long$s125$cCIDr, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s125$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/CombCon_conf_slide150r.png", 1000, 610)
conf_mat(long$s150, "origID", "cCIDr", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Slide 150", 1, cex = 2, adj = -.8,  line = -5)
dev.off()

conf.mat$s150r <- confusionMatrix(factor(sp.abb$Species[match(long$s150$cCIDr, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s150$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/CombCon_conf_digital125r.png", 1000, 700)
conf_mat(long$d125, "origID", "cCIDr", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Digital 125", 1, cex = 2, adj = -.8,  line = -5)
dev.off()

conf.mat$d125r <- confusionMatrix(factor(sp.abb$Species[match(long$d125$cCIDr, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$d125$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/CombCon_conf_digital150r.png", 1000, 610)
conf_mat(long$d150, "origID", "cCIDr", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Digital 150", 1, cex = 2, adj = -.8,  line = -5)
dev.off()

conf.mat$d150r <- confusionMatrix(factor(sp.abb$Species[match(long$d150$cCIDr, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$d150$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# and for the separate consensus values
png("ASFigures/SepCon_conf_slide125r.png", 1000, 700)
conf_mat(long$s125, "origID", "sCIDr", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Sep Slide 125", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

conf.mat$s125scr <- confusionMatrix(factor(sp.abb$Species[match(long$s125$sCIDr, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s125$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/SepCon_conf_slide150r.png", 1000, 610)
conf_mat(long$s150, "origID", "sCIDr", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Sep Slide 150", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

conf.mat$s150scr <- confusionMatrix(factor(sp.abb$Species[match(long$s150$sCIDr, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s150$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/SepCon_conf_digital125r.png", 1000, 700)
conf_mat(long$d125, "origID", "dCIDr", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Sep Digital 125", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

conf.mat$d125scr <- confusionMatrix(factor(sp.abb$Species[match(long$d125$dCIDr, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$d125$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/SepCon_conf_digital150r.png", 1000, 610)
conf_mat(long$d150, "origID", "dCIDr", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Sep Digital 150", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

conf.mat$d150scr <- confusionMatrix(factor(sp.abb$Species[match(long$d150$dCIDr, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$d150$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# 5c. Influence of no consensus (nca) on the confusion matrices ---------------------------------
# influence of the alphabetical ordering
png("ASFigures/CombCon_conf_slide125nc.png", 1000, 700)
conf_mat(long$s125, "origID", "cCIDnc", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nc')), ], abb.end = c("na", "nca"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext(expression(paste("Slide >125", mu, "m")), 1, cex = 2, adj = -1.2,  line = -5)
dev.off()

conf.mat$s125nc <- confusionMatrix(factor(sp.abb$Species[match(long$s125$cCIDnc, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s125$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/CombCon_conf_slide150nc.png", 1000, 610)
conf_mat(long$s150, "origID", "cCIDnc", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nc')), ], abb.end = c("na", "nca"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext(expression(paste("Slide >150", mu, "m")), 1, cex = 2, adj = -1.2,  line = -5)
dev.off()

conf.mat$s150nc <- confusionMatrix(factor(sp.abb$Species[match(long$s150$cCIDnc, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s150$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/CombCon_conf_digital125nc.png", 1000, 700)
conf_mat(long$d125, "origID", "cCIDnc", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nc')), ], abb.end = c("na", "nca"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext(expression(paste("Digital >125", mu, "m")), 1, cex = 2, adj = -1.2,  line = -5)
dev.off()

conf.mat$d125nc <- confusionMatrix(factor(sp.abb$Species[match(long$d125$cCIDnc, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$d125$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/CombCon_conf_digital150nc.png", 1000, 610)
conf_mat(long$d150, "origID", "cCIDnc", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nc')), ], abb.end = c("na", "nca"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext(expression(paste("Digital >150", mu, "m")), 1, cex = 2, adj = -1.2,  line = -5)
dev.off()

conf.mat$d150nc <- confusionMatrix(factor(sp.abb$Species[match(long$d150$cCIDnc, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$d150$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# and for the separate consensus values
png("ASFigures/SepCon_conf_slide125nc.png", 1000, 700)
conf_mat(long$s125, "origID", "sCIDnc", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nc')), ], abb.end = c("na", "nca"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Sep Slide 125", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

conf.mat$s125scr <- confusionMatrix(factor(sp.abb$Species[match(long$s125$sCIDnc, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s125$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/SepCon_conf_slide150nc.png", 1000, 610)
conf_mat(long$s150, "origID", "sCIDnc", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nc')), ], abb.end = c("na", "nca"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Sep Slide 150", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

conf.mat$s150scr <- confusionMatrix(factor(sp.abb$Species[match(long$s150$sCIDnc, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s150$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/SepCon_conf_digital125nc.png", 1000, 700)
conf_mat(long$d125, "origID", "dCIDnc", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nc')), ], abb.end = c("na", "nca"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Sep Digital 125", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

conf.mat$d125scr <- confusionMatrix(factor(sp.abb$Species[match(long$d125$dCIDnc, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$d125$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/SepCon_conf_digital150nc.png", 1000, 610)
conf_mat(long$d150, "origID", "dCIDnc", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nc')), ], abb.end = c("na", "nca"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID", palette = "viridis", grid = TRUE)
mtext("Sep Digital 150", 1, cex = 2, adj = -1.3,  line = -5)
dev.off()

conf.mat$d150scr <- confusionMatrix(factor(sp.abb$Species[match(long$d150$dCIDnc, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$d150$origID, sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))


# 6. NMDS -----------------------------------------------------------------

# 6a. Generate a dataframe for choosing colours ---------------------------
# create a dataframe for the colours
mds.col <- data.frame(person = names(full.125)[!(names(full.125) %in% c("Specimen", "cMaxCon", "sMaxCon", "dMaxCon", "cCIDr", "sCIDr", "dCIDr"))])
# add in the school
mds.col$school <- people$School[match(mds.col$person, people$SlideID)]
mds.col$school[mds.col$person %in% people$DigitalID] <- people$School[order(people$DigitalID)[1:sum(!is.na(people$DigitalID))]]
mds.col$school[grep("1[a-z]", mds.col$person)] <- people$School[people$SlideID == "1"][1]
mds.col$school[grep("2[a-z]", mds.col$person)] <- people$School[people$SlideID == "2"][1]
mds.col$school[is.na(mds.col$school)] <- as.character(mds.col$person[is.na(mds.col$school)])
# create a colour for the school
mds.col$sch.col <- NA
mds.col$sch.col[grep("^[1-9]", mds.col$school)] <- mds.col$school[grep("^[1-9]", mds.col$school)]
mds.col$sch.col <- gsub("_.*", "", mds.col$sch.col)
mds.col$sch.col <- as.numeric(mds.col$sch.col)

# add in the paired analyses
mds.col$pair <- NA
for (i in which(!is.na(people$SlideID) & !is.na(people$DigitalID))) {
  mds.col$pair[mds.col$person %in% people$DigitalID[i]] <- i
  mds.col$pair[grep(paste("^",people$SlideID[i], sep = ""), mds.col$person)] <- i
}
mds.col$pair[mds.col$person == "2a"] <- NA
rm(i)

# 6b. NMDS for combined consensus 125 --------------------------------------------------------
# generate a transposed dataframe for the NMDS analysis
trsp <- list()
trsp$CC125f <- data.frame(t(full.125[, !(names(full.125) %in% c("Specimen", "cMaxCon", "sMaxCon", "dMaxCon", "cCIDr", "sCIDr", "dCIDr"))]))

# consider the stress of the NMDS
stress <- list()
stress$CC125f <- rep(NA, 10)
for (i in 1:10) {
  stress$CC125f[i] <- metaMDS(daisy(trsp$CC125f), k = i)$stress
}
plot(stress$CC125f, type = "b")
rm(i)
# looks like between 2 and 3 dimensions would be reasonable

# run the NMDS
nmds <- list()
nmds$CC125f <- metaMDS(daisy(trsp$CC125f), k = 2, trymax = 1000)
stress$CCop125f <- nmds$CC125f$stress

stressplot(nmds$CC125f)

# there is still a small amount of variation in this plot each time it is written
png("ASFigures/CombCon_NMDS_125.png", 600, 600)
par(mar = c(5.1, 5.1, 3.1, 2.1))
plot(nmds$CC125f, type = "n", display = "sites", cex = 1, xlab = "Axis 1", ylab = "Axis 2", las = 1, cex.lab = 2, cex.main = 2.5, cex.axis = 1.4, main = expression(paste("NMDS >125 ", mu, "m")))
points(nmds$CC125f, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
text(nmds$CC125f$points[nchar(rownames(nmds$CC125f$points)) <3, ], labels = rownames(nmds$CC125f$points)[nchar(rownames(nmds$CC125f$points)) <3])
points(nmds$CC125f$points[rownames(nmds$CC125f$points) == "cCID" | rownames(nmds$CC125f$points) == "cSC50" , ], pch = "+")
text(nmds$CC125f$points[rownames(nmds$CC125f$points) == "cCID", 1]+0.03, nmds$CC125f$points[rownames(nmds$CC125f$points) == "cCID", 2]+0.00, labels = "CID", cex = 1.2)
text(nmds$CC125f$points[rownames(nmds$CC125f$points) == "cSC50", 1]+0.02, nmds$CC125f$points[rownames(nmds$CC125f$points) == "cSC50", 2]+0.00, labels = "sC", cex = 1.2)

legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16, pt.cex = 2, cex = 1.5)
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# 6c. NMDS for combined consensus 150 -------------------------------------
trsp$CC150f <- data.frame(t(full.150[, !(names(full.150) %in% c("Specimen", "cMaxCon", "sMaxCon", "dMaxCon", "cCIDr", "sCIDr", "dCIDr"))]))

# consider the stress of the NMDS
stress$CC150f <- rep(NA, 10)
for (i in 1:10) {
  stress$CC150f[i] <- metaMDS(daisy(trsp$CC150f), k = i)$stress
}
plot(stress$CC150f, type = "b")
rm(i)
# looks like between 2 and 3 dimensions would be reasonable, although the breakpoint is less obvious for 150 than it is for 125.

# run the NMDS
nmds$CC150f <- metaMDS(daisy(trsp$CC150f), k = 2, trymax = 1000)
stress$CCop150f <- nmds$CC150f$stress

stressplot(nmds$CC150f)

# the full plot
png("ASFigures/CombCon_NMDS_150.png", 600, 600)
par(mar = c(5.1, 5.1, 3.1, 2.1))
plot(nmds$CC150f, type = "n", display = "sites", cex = 1, xlab = "Axis 1", ylab = "Axis 2", las = 1, cex.lab = 2, cex.main = 2.5, cex.axis = 1.4, main = expression(paste("NMDS >150 ", mu, "m")))
points(nmds$CC150f, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
text(nmds$CC150f$points[nchar(rownames(nmds$CC150f$points)) <3, ], labels = rownames(nmds$CC150f$points)[nchar(rownames(nmds$CC150f$points)) <3])
points(nmds$CC150f$points[rownames(nmds$CC150f$points) == "cCID" | rownames(nmds$CC150f$points) == "cSC50" , ], pch = "+")
text(nmds$CC150f$points[rownames(nmds$CC150f$points) == "cCID", 1]+0.02, nmds$CC150f$points[rownames(nmds$CC150f$points) == "cCID", 2], labels = "CID", cex = 1.2)
text(nmds$CC150f$points[rownames(nmds$CC150f$points) == "cSC50", 1]+0.025, nmds$CC150f$points[rownames(nmds$CC150f$points) == "cSC50", 2], labels = "sC", cex = 1.2)
legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16, pt.cex = 2, cex = 1.5)
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# focussing on the main section. As noted above, it is better to run this as a new analysis rather than just zoom in, as the influence of the outliers means that the stability of the central points hasn't been tested

# consider the stress of the NMDS
stress$CC150z <- rep(NA, 10)
for (i in 1:10) {
  stress$CC150z[i] <- metaMDS(daisy(trsp$CC150f[!(rownames(trsp$CC150f) %in% c("3", "C", "E", "G")),]), k = i)$stress
}
plot(stress$CC150z, type = "b")
rm(i)
# again between 2 and 3 dimensions is probably reasonable

# run the NMDS
nmds$CC150z <- metaMDS(daisy(trsp$CC150f[!(rownames(trsp$CC150f) %in% c("3", "C", "E", "G")),]), k = 2, trymax = 1000)
stress$CCop150z <- nmds$CC150z$stress
stressplot(nmds$CC150z)

# the zoomed plot
png("ASFigures/CombCon_NMDS_150_no_out.png", 600, 600)
par(mar = c(5.1, 5.1, 3.1, 2.1))
plot(nmds$CC150z, type = "n", display = "sites", cex = 1, xlab = "Axis 1", ylab = "Axis 2", las = 1, cex.lab = 2, cex.main = 2.5, cex.axis = 1.4, main = expression(paste("NMDS >150 ", mu, "m: no outliers")), xlim = c(-0.25, 0.25))

points(nmds$CC150z, pch = 21, cex = 4, col = mds.col$pair[!(mds.col$person %in% c("3", "C", "E", "G"))], bg = brewer.pal(5, "Set2")[mds.col$sch.col[!(mds.col$person %in% c("3", "C", "E", "G"))]], lwd = 2)
text(nmds$CC150z$points[nchar(rownames(nmds$CC150z$points)) <3, ], labels = rownames(nmds$CC150z$points)[nchar(rownames(nmds$CC150z$points)) <3])
points(nmds$CC150z$points[rownames(nmds$CC150z$points) == "cCID" | rownames(nmds$CC150z$points) == "cSC50" , ], pch = "+")
text(nmds$CC150z$points[rownames(nmds$CC150z$points) == "cCID", 1]+0.012, nmds$CC150z$points[rownames(nmds$CC150z$points) == "cCID", 2]+0.008, labels = "CID", cex = 1.2)
text(nmds$CC150z$points[rownames(nmds$CC150z$points) == "cSC50", 1]+0.012, nmds$CC150z$points[rownames(nmds$CC150z$points) == "cSC50", 2]+0.008, labels = "sC", cex = 1.2)
legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pt.cex = 2, cex = 1.5, pch = 16)
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# output the scree plots
png("ASFigures/CombCon_Scree_plots.png")
par(mfrow = c(2, 2))
plot(stress$CC125f, type = "b", bty = "l", las = 1, xlab = "Dimensions", ylab = "Stress", main = "NMDS 125")
plot(stress$CC150f, type = "b", bty = "l", las = 1, xlab = "Dimensions", ylab = "Stress", main = "NMDS 150")
plot(stress$CC150z, type = "b", bty = "l", las = 1, xlab = "Dimensions", ylab = "Stress", main = "NMDS 150 zoomed")
par(mfrow = c(1,1))
dev.off()

# 7. Repeated analysis by workers -----------------------------------------
# for 1a / 1b slide 125
head(long$s125)

sum(full.125$`1a` == full.125$`1b`) # 183 or 61% similarity
sum(full.125$`1a` == full.125$cCID); sum(full.125$`1a` == full.125$sCID) # 192 or 64% accuracy (cf. 198)
sum(full.125$`1b` == full.125$cCID); sum(full.125$`1b` == full.125$sCID) # 227 or 76% accuracy (cf. 228)

png("ASFigures/Time/CombCon_conf_125_1aCon.png", 1000, 700)
conf_mat(long$s125[long$s125$Person == "1a", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "1a", ylab = "Consensus", palette = "viridis", grid = TRUE)
dev.off()

conf.mat$s125c1a <- confusionMatrix(factor(sp.abb$Species[match(long$s125$cCID[long$s125$Person == "1a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s125$origID[long$s125$Person == "1a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/Time/CombCon_conf_125_1bCon.png", 1000, 700)
conf_mat(long$s125[long$s125$Person == "1b", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "1b", ylab = "Consensus", palette = "viridis", grid = TRUE)
dev.off()

conf.mat$s125c1b <- confusionMatrix(factor(sp.abb$Species[match(long$s125$cCID[long$s125$Person == "1b"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s125$origID[long$s125$Person == "1b"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# and for 2a / 2b
sum(full.125$`2a` == full.125$`2b`) # 241 or 80% similarity
sum(full.125$`2a` == full.125$cCID); sum(full.125$`2a` == full.125$sCID) # 206 or 69% accuracy (cf. 208)
sum(full.125$`2b` == full.125$cCID); sum(full.125$`2b` == full.125$sCID) # 240 or 80% accuracy (cf. 237)

# again, do this as a confusion matrix
png("ASFigures/Time/comb_Con_conf_125_2aCon.png", 1000, 700)
conf_mat(long$s125[long$s125$Person == "2a", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "2a", ylab = "Consensus", palette = "viridis", grid = TRUE)
dev.off()

conf.mat$s125c2a <- confusionMatrix(factor(sp.abb$Species[match(long$s125$cCID[long$s125$Person == "2a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s125$origID[long$s125$Person == "2a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/Time/comb_Con_conf_125_2bCon.png", 1000, 700)
conf_mat(long$s125[long$s125$Person == "2b", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "2b", ylab = "Consensus", palette = "viridis", grid = TRUE)
dev.off()

conf.mat$s125c2b <- confusionMatrix(factor(sp.abb$Species[match(long$s125$cCID[long$s125$Person == "2b"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s125$origID[long$s125$Person == "2b"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# And for 150
sum(full.150$`1a` == full.150$`1b`) # 204 or 68% similarity
sum(full.150$`1a` == full.150$cCID); sum(full.150$`1a` == full.150$sCID) # 221 or 74% accuracy (cf. 221)
sum(full.150$`1b` == full.150$cCID); sum(full.150$`1b` == full.150$sCID) # 219 or 73% accuracy (cf. 223)

png("ASFigures/Time/comb_Con_conf_150_1aCon.png", 1000, 700)
conf_mat(long$s150[long$s150$Person == "1a", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "1a", ylab = "Consensus", palette = "viridis", grid = TRUE)
dev.off()

conf.mat$s150c1a <- confusionMatrix(factor(sp.abb$Species[match(long$s150$cCID[long$s150$Person == "1a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s150$origID[long$s150$Person == "1a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/Time/comb_Con_conf_150_1bCon.png", 1000, 700)
conf_mat(long$s150[long$s150$Person == "1b", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "1b", ylab = "Consensus", palette = "viridis", grid = TRUE)
dev.off()

conf.mat$s150c1b <- confusionMatrix(factor(sp.abb$Species[match(long$s150$cCID[long$s150$Person == "1b"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s150$origID[long$s150$Person == "1b"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# 2a/2b
sum(full.150$`2a` == full.150$`2b`) # 288 or 96% similarity
sum(full.150$`2a` == full.150$cCID); sum(full.150$`2a` == full.150$sCID) # 257 or 86% accuracy (cf. 255)
sum(full.150$`2b` == full.150$cCID); sum(full.150$`2b` == full.150$sCID) # 257 or 86% accuracy (cf. 255)

png("ASFigures/Time/comb_Con_conf_150_2aCon.png", 1000, 700)
conf_mat(long$s150[long$s150$Person == "2a", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "2a", ylab = "Consensus", palette = "viridis", grid = TRUE)
dev.off()

conf.mat$s150c2a <- confusionMatrix(factor(sp.abb$Species[match(long$s150$cCID[long$s150$Person == "2a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s150$origID[long$s150$Person == "2a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/Time/comb_Con_conf_150_2bCon.png", 1000, 700)
conf_mat(long$s150[long$s150$Person == "2b", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "2b", ylab = "Consensus", palette = "viridis", grid = TRUE)
dev.off()

conf.mat$s150c2b <- confusionMatrix(factor(sp.abb$Species[match(long$s150$cCID[long$s150$Person == "2b"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$s150$origID[long$s150$Person == "2b"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# 8. Digital vs. slides ---------------------------------------------------

# 8a. For 125 -------------------------------------------------------------
long$f125 <- rbind(long$s125, long$d125)

# 125 2b vs. A
sum(full.125$`2b` == full.125$A) # 171 or 57% similarity
sum(full.125$`2b` == full.125$cCID); sum(full.125$`2b` == full.125$sCID) # 240 or 80% (cf. 237)
sum(full.125$`A` == full.125$cCID); sum(full.125$`A` == full.125$dCID) # 175 or 58% accuracy (cf. 181)

png("ASFigures/DigitalSlide/CombCon_conf_125_2bA.png", 1000, 700)
conf_mat(long$f125, "origID", axis.col = "Person", axis1 = "2b", axis2 = "A", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, grid = TRUE, palette = "viridis")
dev.off()

conf.mat$c125A2b <- confusionMatrix(factor(sp.abb$Species[match(long$f125$origID[long$f125$Person == "A"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f125$origID[long$f125$Person == "2b"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/DigitalSlide/CombCon_conf_125_2bCon.png", 1000, 700)
conf_mat(long$f125[long$f125$Person == "2b", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "2b", ylab = "Consensus", grid = TRUE, palette = "viridis")
dev.off()

conf.mat$c125c2b <- confusionMatrix(factor(sp.abb$Species[match(long$f125$cCID[long$f125$Person == "2b"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f125$origID[long$f125$Person == "2b"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/DigitalSlide/CombCon_conf_125_ACon.png", 1000, 700)
conf_mat(long$f125[long$f125$Person == "A", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "A", ylab = "Consensus", grid = TRUE, palette = "viridis")
dev.off()

conf.mat$c125cA <- confusionMatrix(factor(sp.abb$Species[match(long$f125$cCID[long$f125$Person == "A"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f125$origID[long$f125$Person == "A"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# 125 6 vs. F
sum(full.125$`6` == full.125$F) # 187 or 62% similarity
sum(full.125$`6` == full.125$cCID); sum(full.125$`6` == full.125$sCID) # 212 or 71% (cf. 211)
sum(full.125$`F` == full.125$cCID); sum(full.125$`F` == full.125$dCID) # 229 or 76% accuracy (cf. 226)

png("ASFigures/DigitalSlide/CombCon_conf_125_6F.png", 1000, 700)
conf_mat(long$f125, "origID", axis.col = "Person", axis1 = "6", axis2 = "F", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, grid = TRUE, palette = "viridis")
dev.off()

conf.mat$c125F6 <- confusionMatrix(factor(sp.abb$Species[match(long$f125$origID[long$f125$Person == "F"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f125$origID[long$f125$Person == "6"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/DigitalSlide/CombCon_conf_125_FCon.png", 1000, 700)
conf_mat(long$f125[long$f125$Person == "F", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "F", ylab = "Consensus", grid = TRUE, palette = "viridis")
dev.off()

conf.mat$c125cF <- confusionMatrix(factor(sp.abb$Species[match(long$f125$cCID[long$f125$Person == "F"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f125$origID[long$f125$Person == "F"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/DigitalSlide/CombCon_conf_125_6Con.png", 1000, 700)
conf_mat(long$f125[long$f125$Person == "6", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "6", ylab = "Consensus", grid = TRUE, palette = "viridis")
dev.off()

conf.mat$c125c6 <- confusionMatrix(factor(sp.abb$Species[match(long$f125$cCID[long$f125$Person == "6"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f125$origID[long$f125$Person == "6"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# 125 9 vs. G
sum(full.125$`9` == full.125$G) # 232 or 77% similarity
sum(full.125$`9` == full.125$cCID); sum(full.125$`9` == full.125$sCID) # 212 or 71% accuracy (cf. 206)
sum(full.125$`G` == full.125$cCID); sum(full.125$`G` == full.125$dCID) # 160 or 53% accuracy (cf. 157)

png("ASFigures/DigitalSlide/CombCon_conf_125_9G.png", 1000, 700)
conf_mat(long$f125, "origID", axis.col = "Person", axis1 = "9", axis2 = "G", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, grid = TRUE, palette = "viridis")
dev.off()

conf.mat$c125G9 <- confusionMatrix(factor(sp.abb$Species[match(long$f125$origID[long$f125$Person == "G"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f125$origID[long$f125$Person == "9"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/DigitalSlide/CombCon_conf_125_GCon.png", 1000, 700)
conf_mat(long$f125[long$f125$Person == "G", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "G", ylab = "Consensus", grid = TRUE, palette = "viridis")
dev.off()

conf.mat$c125cG <- confusionMatrix(factor(sp.abb$Species[match(long$f125$cCID[long$f125$Person == "G"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f125$origID[long$f125$Person == "G"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/DigitalSlide/CombCon_conf_125_9Con.png", 1000, 700)
conf_mat(long$f125[long$f125$Person == "9", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "9", ylab = "Consensus", grid = TRUE, palette = "viridis")
dev.off()

conf.mat$c125c9 <- confusionMatrix(factor(sp.abb$Species[match(long$f125$cCID[long$f125$Person == "9"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f125$origID[long$f125$Person == "9"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# 8b. for 150 -----------------------------------------------------------------
long$f150 <- rbind(long$s150, long$d150)

# 150 2b vs. A
sum(full.150$`2b` == full.150$A) # 232 or 77% similarity
sum(full.150$`2b` == full.150$cCID); sum(full.150$`2b` == full.150$sCID) # 257 or 86% accuracy (cf. 255)
sum(full.150$`A` == full.150$cCID); sum(full.150$`A` == full.150$dCID) # 241 or 80% accuracy  (cf. 246)

png("ASFigures/DigitalSlide/CombCon_conf_150_2bA.png", 1000, 700)
conf_mat(long$f150, "origID", axis.col = "Person", axis1 = "2b", axis2 = "A", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, grid = TRUE, palette = "viridis")
dev.off()

conf.mat$c150A2b <- confusionMatrix(factor(sp.abb$Species[match(long$f150$origID[long$f150$Person == "A"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f150$origID[long$f150$Person == "2b"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/DigitalSlide/CombCon_conf_150_2bCon.png", 1000, 700)
conf_mat(long$f150[long$f150$Person == "2b", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "2b", ylab = "Consensus", palette = "viridis", grid = TRUE)
dev.off()

conf.mat$c150c2b <- confusionMatrix(factor(sp.abb$Species[match(long$f150$cCID[long$f150$Person == "2b"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f150$origID[long$f150$Person == "2b"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/DigitalSlide/CombCon_conf_150_ACon.png", 1000, 700)
conf_mat(long$f150[long$f150$Person == "A", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "A", ylab = "Consensus", palette = "viridis", grid = TRUE)
dev.off()

conf.mat$c150cA <- confusionMatrix(factor(sp.abb$Species[match(long$f150$cCID[long$f150$Person == "A"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f150$origID[long$f150$Person == "A"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# 150 6 vs. F
sum(full.150$`6` == full.150$F) # 216 or 72% similarity
sum(full.150$`6` == full.150$cCID); sum(full.150$`6` == full.150$sCID) # 246 or 82% accuracy (cf. 244)
sum(full.150$`F` == full.150$cCID); sum(full.150$`F` == full.150$dCID) # 244 or 81% accuracy (cf. 242)

png("ASFigures/DigitalSlide/CombCon_conf_150_6F.png", 1000, 700)
conf_mat(long$f150, "origID", axis.col = "Person", axis1 = "6", axis2 = "F", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, grid = TRUE, palette = "viridis")
dev.off()

conf.mat$c150F6 <- confusionMatrix(factor(sp.abb$Species[match(long$f150$origID[long$f150$Person == "F"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f150$origID[long$f150$Person == "6"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/DigitalSlide/CombCon_conf_150_FCon.png", 1000, 700)
conf_mat(long$f150[long$f150$Person == "F", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "F", ylab = "Consensus", palette = "viridis", grid = TRUE)
dev.off()

conf.mat$c150cF <- confusionMatrix(factor(sp.abb$Species[match(long$f150$cCID[long$f150$Person == "F"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f150$origID[long$f150$Person == "F"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/DigitalSlide/CombCon_conf_150_6Con.png", 1000, 700)
conf_mat(long$f150[long$f150$Person == "6", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "6", ylab = "Consensus", palette = "viridis", grid = TRUE)
dev.off()

conf.mat$c150c6 <- confusionMatrix(factor(sp.abb$Species[match(long$f150$cCID[long$f150$Person == "6"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f150$origID[long$f150$Person == "6"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))


# 150 9 vs. G
sum(full.150$`9` == full.150$G) # 212 or 71% similarity
sum(full.150$`9` == full.150$cCID); sum(full.150$`9` == full.150$sCID) # 216 or 72% accuracy (cf. 220)
sum(full.150$`G` == full.150$cCID); sum(full.150$`G` == full.150$dCID) # 148 or 49% accuracy (cf. 149)

png("ASFigures/DigitalSlide/CombCon_conf_150_9G.png", 1000, 700)
conf_mat(long$f150, "origID", axis.col = "Person", axis1 = "9", axis2 = "G", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, palette = "viridis", grid = TRUE)
dev.off()

conf.mat$c150G9 <- confusionMatrix(factor(sp.abb$Species[match(long$f150$origID[long$f150$Person == "G"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f150$origID[long$f150$Person == "9"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/DigitalSlide/CombCon_conf_150_GCon.png", 1000, 700)
conf_mat(long$f150[long$f150$Person == "G", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "G", ylab = "Consensus", palette = "viridis", grid = TRUE)
dev.off()

conf.mat$c150cG <- confusionMatrix(factor(sp.abb$Species[match(long$f150$cCID[long$f150$Person == "G"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f150$origID[long$f150$Person == "G"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/DigitalSlide/CombCon_conf_150_9Con.png", 1000, 700)
conf_mat(long$f150[long$f150$Person == "9", ], "origID", "cCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), axes.same = TRUE, xlab = "9", ylab = "Consensus", palette = "viridis", grid = TRUE)
dev.off()

conf.mat$c150c9 <- confusionMatrix(factor(sp.abb$Species[match(long$f150$cCID[long$f150$Person == "9"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f150$origID[long$f150$Person == "9"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# 8c. Consensus comparisons -----------------------------------------------

sum(full.125$sCID == full.125$dCID) # 234 or 78% accuracy
sum(full.150$sCID == full.150$dCID) # 248 or 83% accuracy

# plotting the consensus' against each other
png("ASFigures/DigitalSlide/conf_125_Con.png", 800, 700)
conf_mat(long$f125[long$f125$Person == "1a",], "dCID", "sCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), xlab = "Digital", ylab = "Slide", palette = "viridis", grid = TRUE)
text(-1, 1, "Consensus 125", cex = 1.5)
dev.off()

conf.mat$c125 <- confusionMatrix(factor(sp.abb$Species[match(long$f125$sCID[long$f125$Person == "1a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f125$dCID[long$f125$Person == "1a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/DigitalSlide/conf_150_Con.png", 800, 700)
conf_mat(long$f150[long$f150$Person == "1a",], "dCID", "sCID", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nca', 'nc')), ], abb.end = c("na"), xlab = "Digital", ylab = "Slide", palette = "viridis", grid = TRUE)
text(-1, 1, "Consensus 150", cex = 1.5)
dev.off()

conf.mat$c150 <- confusionMatrix(factor(sp.abb$Species[match(long$f150$sCID[long$f150$Person == "1a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f150$dCID[long$f150$Person == "1a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

# and with the tied IDs excluded
png("ASFigures/DigitalSlide/conf_125_Con_nc.png", 800, 700)
conf_mat(long$f125[long$f125$Person == "1a",], "dCIDnc", "sCIDnc", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nc')), ], abb.end = c("nca"), xlab = "Digital", ylab = "Slide", palette = "viridis", grid = TRUE)
text(-1, 1, "Consensus 125", cex = 1.5)
dev.off()

conf.mat$c125nc <- confusionMatrix(factor(sp.abb$Species[match(long$f125$sCIDnc[long$f125$Person == "1a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f125$dCIDnc[long$f125$Person == "1a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))

png("ASFigures/DigitalSlide/conf_150_Con_nc.png", 800, 700)
conf_mat(long$f150[long$f150$Person == "1a",], "dCIDnc", "sCIDnc", spec.abb = sp.abb[!(sp.abb$Abbreviation %in% c('nc')), ], abb.end = c("nca"), xlab = "Digital", ylab = "Slide", palette = "viridis", grid = TRUE)
text(-1, 1, "Consensus 150", cex = 1.5)
dev.off()

conf.mat$c150nc <- confusionMatrix(factor(sp.abb$Species[match(long$f150$sCIDnc[long$f150$Person == "1a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]), factor(sp.abb$Species[match(long$f150$dCIDnc[long$f150$Person == "1a"], sp.abb$Abbreviation)], levels = sp.abb$Species[order(sp.abb$confOrder)]))


# 9. SST ------------------------------------------------------------------

# 9a. Extract the species level data for the ANN --------------------------
full.ANNsp <- full.150sp

# add in the ANN species names
full.ANNsp$order <- sp.abb$ANNorder[match(full.150sp$species, sp.abb$Abbreviation)]
full.ANNsp$species <- sp.abb$ANNspecies[match(full.150sp$species, sp.abb$Abbreviation)]

# remove NAs (i.e. species not in ANN)
full.ANNsp <- na.omit(full.ANNsp)

# merge menardii and tumida
full.ANNsp[which(full.ANNsp$species == "Globorotalia menardii + tumida")[1], 2:(ncol(full.ANNsp) - 1)] <- colSums(full.ANNsp[full.ANNsp$species == "Globorotalia menardii + tumida", 2:(ncol(full.ANNsp) - 1)])
full.ANNsp <- full.ANNsp[-which(full.ANNsp$species == "Globorotalia menardii + tumida")[2],]

# reorder
full.ANNsp <- full.ANNsp[order(full.ANNsp$order), ]

# transpose
full.ANNsp.t <- as.data.frame(t(full.ANNsp[2:(ncol(full.ANNsp) - 1)])) # don't include order or the species names
colnames(full.ANNsp.t) <- full.ANNsp$species # add species names as column headings
head(full.ANNsp.t)

# calculate relative abundances
full.ANNsp.t <- full.ANNsp.t / rowSums(full.ANNsp.t) * 100

# output the data
write.csv(full.ANNsp.t, "ASOutputs/RelativeAbun_ANN.csv")

# check the percentage agreement when only these species are considered
tmp <- full.150[sp.abb$Abbreviation[!is.na(sp.abb$ANNspecies)] %in% full.150$cCID, ]
tmp2 <- apply(tmp[,col.nam$c150], 2, function(x) sum(x == tmp$cCID) / nrow(tmp) * 100)
summary(accuracy$cPtA150 - tmp2)
# on average the agreement is slightly lower for the species included in the ANN. 
rm(tmp, tmp2)

# 9b. Plot up the SST results ------------------------------------------------
# this was only done 150 size fraction as that is what ANN works on 
# specify the rownames
row.nam <- list()
row.nam$div <- which(nchar(divTemp$Person) < 3)
row.nam$s125 <- which(nchar(divTemp$Person) < 3 & divTemp$Analysis == "Slide" & divTemp$Size == 125)
row.nam$s150 <- which(nchar(divTemp$Person) < 3 & divTemp$Analysis == "Slide" & divTemp$Size == 150)
row.nam$d125 <- which(nchar(divTemp$Person) < 3 & divTemp$Analysis == "Digital" & divTemp$Size == 125)
row.nam$d150 <- which(nchar(divTemp$Person) < 3 & divTemp$Analysis == "Digital" & divTemp$Size == 150)
row.nam$s125c <- which(divTemp$Person == "consensus" & divTemp$Analysis == "Slide" & divTemp$Size == 125)
row.nam$s150c <- which(divTemp$Person == "consensus" & divTemp$Analysis == "Slide" & divTemp$Size == 150)
row.nam$d125c <- which(divTemp$Person == "consensus" & divTemp$Analysis == "Digital" & divTemp$Size == 125)
row.nam$d150c <- which(divTemp$Person == "consensus" & divTemp$Analysis == "Digital" & divTemp$Size == 150)
row.nam$c125c <- which(divTemp$Person == "consensus" & divTemp$Analysis == "Con" & divTemp$Size == 125)
row.nam$c150c <- which(divTemp$Person == "consensus" & divTemp$Analysis == "Con" & divTemp$Size == 150)
row.nam$sSST <- which(divTemp$Analysis == "Slide" & divTemp$Size == 150)
row.nam$dSST <- which(divTemp$Analysis == "Digital" & divTemp$Size == 150)

# generate an error bar function
err_bar <- function(mean, sd, xpos, length = 0.05, col = 1) {
  for(i in 1:length(mean)) {
    if (length(col) > 1) {
      arrows(xpos[i], mean[i] - sd[i], xpos[i], mean[i] + sd[i], angle = 90, code = 3, length = length, col = col[i])
    } else {
      arrows(xpos[i], mean[i] - sd[i], xpos[i], mean[i] + sd[i], angle = 90, code = 3, length = length, col = col)
    }
    
  }
}

# plot the data
png("ASFigures/CombCon_SST.png", 800, 500)
par(mar = c(5.1, 5.1, 4.1, 2.1))
# points
with(divTemp[divTemp$Size == 150,], plot(1:26, SST10m[match(ord.div, Person)], pch = 16, xaxt = "n", xlab = "Participant", ylab = expression(paste("SST / ", degree, "C")), col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(20.75, 24.35), cex.lab = 1.5, las = 1, cex.axis = 1.1, type = "n"))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
# consensus values
with(divTemp[row.nam$c150c,], lines(x = c(0, 27), y = rep(SST10m, each = 2), col = "green4"))
with(divTemp[row.nam$c150c,], rect(0, SST10m - SD, 27, SST10m + SD, col = rgb(0, 139, 0, 75, maxColorValue = 255), border = NA))
# mean values
lines(x = c(0, 20.5), y = rep(mean(divTemp$SST10m[row.nam$s150]), each = 2))
lines(x = c(20.5, 27), y = rep(mean(divTemp$SST10m[row.nam$d150]), each = 2), col = 4)
abline(v = 20.5, col = "grey 50")
# add points
with(divTemp[divTemp$Size == 150,], points(1:26, SST10m[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))

# error bars / actual  / legend
with(divTemp[divTemp$Size == 150,], err_bar(SST10m[match(ord.div, Person)], SD[match(ord.div, Person)], 1:26, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
abline(h = 21.76, col = "green4", lwd = 2)
text(23.25, 21.9, "WOA 1998", cex = 1.3, col = "green4")
legend("topleft", legend = c(expression(paste("Slide >150 ", mu, "m")), expression(paste("Digital >150 ", mu, "m"))), pch = 16, col = c(1, 4), pt.cex = 1.5, cex = 1.1)
text(1, 20.75, "Slide", cex = 1.3)
text(21.75, 20.75, "Digital", cex = 1.3, col = "blue")

par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# 10. Diversity ------------------------------------------------------------

# 10a. Calculating diversity -----------------------------------------------
head(full.125sp)

# calculate these without 'na' and 'nc'
# richness 
divTemp$Richness <- NA
divTemp$Richness[divTemp$Size == 125] <- specnumber(t(full.125sp[!grepl("nc|na", full.125sp$species), c("cCID", "sCID", names(full.125sp)[col.nam$s125], "dCID", names(full.125sp)[col.nam$d125])]))
divTemp$Richness[divTemp$Size == 150] <- specnumber(t(full.150sp[!grepl("nc|na", full.150sp$species), c("cCID", "sCID", names(full.150sp)[col.nam$s150], "dCID", names(full.150sp)[col.nam$d150])]))

# ShannonWiener
divTemp$ShannonWiener <- NA
divTemp$ShannonWiener[divTemp$Size == 125] <- diversity(t(full.125sp[!grepl("nc|na", full.125sp$species), c("cCID", "sCID", names(full.125sp)[col.nam$s125], "dCID", names(full.125sp)[col.nam$d125])]))
divTemp$ShannonWiener[divTemp$Size == 150] <- diversity(t(full.150sp[!grepl("nc|na", full.150sp$species), c("cCID", "sCID", names(full.150sp)[col.nam$s150], "dCID", names(full.150sp)[col.nam$d150])]))

# Dominance
divTemp$Dominance <- NA
divTemp$Dominance[divTemp$Size == 125] <- (1 - diversity(t(full.125sp[!grepl("nc|na", full.125sp$species), c("cCID", "sCID", names(full.125sp)[col.nam$s125], "dCID", names(full.125sp)[col.nam$d125])]), index = "simpson"))
divTemp$Dominance[divTemp$Size == 150] <- (1 - diversity(t(full.150sp[!grepl("nc|na", full.150sp$species), c("cCID", "sCID", names(full.150sp)[col.nam$s150], "dCID", names(full.150sp)[col.nam$d150])]), index = "simpson"))

# 10b. Plotting diversity combined consensus --------------------------------------------------
# Plotting diversity
png("ASFigures/div_CombCon_Richness.png", 800, 500)
par(mar = c(5.1, 5.1, 4.1, 2.1))
# points
with(divTemp[divTemp$Size == 125,], plot(1:26, Richness[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Richness", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(14, 30), cex.lab = 1.5, las = 1, cex.axis = 1.1, cex = 1.5))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
with(divTemp[divTemp$Size == 150,], points(1:26, Richness[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), cex = 1.5))
# consensus lines
abline(h = divTemp$Richness[row.nam$c125c], col = "green4", lty = 2)
abline(h = divTemp$Richness[row.nam$c150c], col = "green4")
# mean values
lines(x = c(20.5, 27), y = rep(mean(divTemp$Richness[row.nam$d125]), 2), col = "blue", lty = 2)
lines(x = c(0, 20.5), y = rep(mean(divTemp$Richness[row.nam$s125]), 2), lty = 2)
lines(x = c(20.5, 27), y = rep(mean(divTemp$Richness[row.nam$d150]), 2), col = "blue")
lines(x = c(0, 20.5), y = rep(mean(divTemp$Richness[row.nam$s150]), 2))
abline(v = 20.5, col = "grey 50")
# legend
legend("topleft", legend = c(expression(paste("Slide >125 ", mu, "m")), expression(paste("Slide >150 ", mu, "m")), expression(paste("Digital >125 ", mu, "m")), expression(paste("Digital >150 ", mu, "m"))), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4), pt.cex = 1.5, cex = 1.1)
text(1, 14, "Slide", cex = 1.3)
text(21.75, 14, "Digital", cex = 1.3, col = "blue")

par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("ASFigures/div_CombCon_Dominance.png", 800, 500)
# points
par(mar = c(5.1, 5.1, 4.1, 2.1))
with(divTemp[divTemp$Size == 125,], plot(1:26, Dominance[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Dominance", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(0.1, 0.22), cex.lab = 1.5, las = 1, cex.axis = 1.1, cex = 1.5))
with(divTemp[divTemp$Size == 150,], points(1:26, Dominance[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), cex = 1.5))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
# consensus lines
abline(h = divTemp$Dominance[row.nam$c125c], col = "green4", lty = 2)
abline(h = divTemp$Dominance[row.nam$c150c], col = "green4")
# mean values
lines(x = c(20.5, 27), y = rep(mean(divTemp$Dominance[row.nam$d125]), 2), col = "blue", lty = 2)
lines(x = c(0, 20.5), y = rep(mean(divTemp$Dominance[row.nam$s125]), 2), lty = 2)
lines(x = c(20.5, 27), y = rep(mean(divTemp$Dominance[row.nam$d150]), 2), col = "blue")
lines(x = c(0, 20.5), y = rep(mean(divTemp$Dominance[row.nam$s150]), 2))
abline(v = 20.5, col = "grey 50")
# legend
legend("topleft", legend = c(expression(paste("Slide >125 ", mu, "m")), expression(paste("Slide >150 ", mu, "m")), expression(paste("Digital >125 ", mu, "m")), expression(paste("Digital >150 ", mu, "m"))), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4), pt.cex = 1.5, cex = 1.1)
text(1, 0.1, "Slide", cex = 1.3)
text(21.75, 0.1, "Digital", cex = 1.3, col = "blue")

par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("ASFigures/div_CombCon_ShannonWiener.png", 800, 500)
# 125
par(mar = c(5.1, 5.1, 4.1, 2.1))
# points
with(divTemp[divTemp$Size == 125,], plot(1:26, ShannonWiener[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Shannon-Wiener Diversity", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(1.95, 2.6), cex.lab = 1.5, las = 1, cex.axis = 1.1, cex = 1.5))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
with(divTemp[divTemp$Size == 150,], points(1:26, ShannonWiener[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), cex = 1.5))
# consensus lines
abline(h = divTemp$ShannonWiener[row.nam$c125c], col = "green4", lty = 2)
abline(h = divTemp$ShannonWiener[row.nam$c150c], col = "green4")
# mean values
lines(x = c(20.5, 27), y = rep(mean(divTemp$ShannonWiener[row.nam$d125]), 2), col = "blue", lty = 2)
lines(x = c(0, 20.5), y = rep(mean(divTemp$ShannonWiener[row.nam$s125]), 2), lty = 2)
lines(x = c(20.5, 27), y = rep(mean(divTemp$ShannonWiener[row.nam$d150]), 2), col = "blue")
lines(x = c(0, 20.5), y = rep(mean(divTemp$ShannonWiener[row.nam$s150]), 2))
abline(v = 20.5, col = "grey 50")
# legend
legend("topleft", legend = c(expression(paste("Slide >125 ", mu, "m")), expression(paste("Slide >150 ", mu, "m")), expression(paste("Digital >125 ", mu, "m")), expression(paste("Digital >150 ", mu, "m"))), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4), pt.cex = 1.5, cex = 1.1)
text(1, 1.95, "Slide", cex = 1.3)
text(21.75, 1.95, "Digital", cex = 1.3, col = "blue")

par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("ASFigures/Div_pairs.png", 600, 600)
pairs(divTemp[, c("Richness", "ShannonWiener", "Dominance")], pch = c(1, 16)[(divTemp$Size == 150) + 1], col = c(2, 4, 1)[as.factor(divTemp$Analysis)], cex = 1.8, cex.axis = 1.6, las = 1)
dev.off()


# 10b. Plotting diversity separate consensus --------------------------------------------------
# Plotting diversity
png("ASFigures/div_Sep_Richness.png", 800, 500)
par(mar = c(5.1, 5.1, 4.1, 2.1))
# points
with(divTemp[divTemp$Size == 125,], plot(1:26, Richness[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Richness", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(14, 30), cex.lab = 1.5, las = 1, cex.axis = 1.1))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
with(divTemp[divTemp$Size == 150,], points(1:26, Richness[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
# consensus lines
lines(x = c(20.5, 27), y = rep(mean(divTemp$Richness[row.nam$d125c]), 2), col = "green4", lty = 2)
lines(x = c(0, 20.5), y = rep(mean(divTemp$Richness[row.nam$s125c]), 2), col = "green4", lty = 2)
lines(x = c(20.5, 27), y = rep(mean(divTemp$Richness[row.nam$d150c]), 2), col = "green4")
lines(x = c(0, 20.5), y = rep(mean(divTemp$Richness[row.nam$s150]), 2), col = "green4")

# mean values
lines(x = c(20.5, 27), y = rep(mean(divTemp$Richness[row.nam$d125]), 2), col = "blue", lty = 2)
lines(x = c(0, 20.5), y = rep(mean(divTemp$Richness[row.nam$s125]), 2), lty = 2)
lines(x = c(20.5, 27), y = rep(mean(divTemp$Richness[row.nam$d150]), 2), col = "blue")
lines(x = c(0, 20.5), y = rep(mean(divTemp$Richness[row.nam$s150]), 2))
abline(v = 20.5, col = "grey 50")
# legend
legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
text(1, 14, "Slide", cex = 1.3)
text(21.75, 14, "Digital", cex = 1.3, col = "blue")

par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("ASFigures/div_Sep_Dominance.png", 800, 500)
# points
par(mar = c(5.1, 5.1, 4.1, 2.1))
with(divTemp[divTemp$Size == 125,], plot(1:26, Dominance[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Dominance", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(0.1, 0.22), cex.lab = 1.5, las = 1, cex.axis = 1.1))
with(divTemp[divTemp$Size == 150,], points(1:26, Dominance[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
# consensus lines
lines(x = c(20.5, 27), y = rep(mean(divTemp$Dominance[row.nam$d125c]), 2), col = "green4", lty = 2)
lines(x = c(0, 20.5), y = rep(mean(divTemp$Dominance[row.nam$s125c]), 2), col = "green4", lty = 2)
lines(x = c(20.5, 27), y = rep(mean(divTemp$Dominance[row.nam$d150c]), 2), col = "green4")
lines(x = c(0, 20.5), y = rep(mean(divTemp$Dominance[row.nam$s150]), 2), col = "green4")
# mean values
lines(x = c(20.5, 27), y = rep(mean(divTemp$Dominance[row.nam$d125]), 2), col = "blue", lty = 2)
lines(x = c(0, 20.5), y = rep(mean(divTemp$Dominance[row.nam$s125]), 2), lty = 2)
lines(x = c(20.5, 27), y = rep(mean(divTemp$Dominance[row.nam$d150]), 2), col = "blue")
lines(x = c(0, 20.5), y = rep(mean(divTemp$Dominance[row.nam$s150]), 2))
abline(v = 20.5, col = "grey 50")
# legend
legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
text(1, 0.1, "Slide", cex = 1.3)
text(21.75, 0.1, "Digital", cex = 1.3, col = "blue")

par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("ASFigures/div_Sep_ShannonWiener.png", 800, 500)
# 125
par(mar = c(5.1, 5.1, 4.1, 2.1))
# points
with(divTemp[divTemp$Size == 125,], plot(1:26, ShannonWiener[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Shannon-Wiener Diversity", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(1.95, 2.6), cex.lab = 1.5, las = 1, cex.axis = 1.1))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
with(divTemp[divTemp$Size == 150,], points(1:26, ShannonWiener[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
# consensus lines
lines(x = c(20.5, 27), y = rep(mean(divTemp$ShannonWiener[row.nam$d125c]), 2), col = "green4", lty = 2)
lines(x = c(0, 20.5), y = rep(mean(divTemp$ShannonWiener[row.nam$s125c]), 2), col = "green4", lty = 2)
lines(x = c(20.5, 27), y = rep(mean(divTemp$ShannonWiener[row.nam$d150c]), 2), col = "green4")
lines(x = c(0, 20.5), y = rep(mean(divTemp$ShannonWiener[row.nam$s150]), 2), col = "green4")
# mean values
lines(x = c(20.5, 27), y = rep(mean(divTemp$ShannonWiener[row.nam$d125]), 2), col = "blue", lty = 2)
lines(x = c(0, 20.5), y = rep(mean(divTemp$ShannonWiener[row.nam$s125]), 2), lty = 2)
lines(x = c(20.5, 27), y = rep(mean(divTemp$ShannonWiener[row.nam$d150]), 2), col = "blue")
lines(x = c(0, 20.5), y = rep(mean(divTemp$ShannonWiener[row.nam$s150]), 2))
abline(v = 20.5, col = "grey 50")
# legend
legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
text(1, 1.95, "Slide", cex = 1.3)
text(21.75, 1.95, "Digital", cex = 1.3, col = "blue")

par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# 10c. Diversity comparisons by size -----------------------------------------
# the other variables show mixed results
# richness r2 = 0.1849, p = 0.0283*
summary(lm(divTemp$Richness[c(row.nam$s150, row.nam$d150)] ~ divTemp$Richness[c(row.nam$s125, row.nam$d125)]))
# dominance r2 = 0.0551, p = 0.248
summary(lm(divTemp$Dominance[c(row.nam$s150, row.nam$d150)] ~ divTemp$Dominance[c(row.nam$s125, row.nam$d125)]))
# Shannon-Wiener r2 = 0.1558, p = 0.04595*
summary(lm(divTemp$ShannonWiener[c(row.nam$s150, row.nam$d150)] ~ divTemp$ShannonWiener[c(row.nam$s125, row.nam$d125)]))


# 10d. Sensitivity to alphabetical order ---------------------------------
# checking diversity estimates from the three different consensus values
# 125 
divTemp$Richness[row.nam$c125c]; divTemp$Richness[row.nam$s125c]; divTemp$Richness[row.nam$d125c] # 22 (s - same, d - 23)
specnumber(table(full.125$cCIDr)); specnumber(table(full.125$sCIDr)); specnumber(table(full.125$dCIDr)) # 22 (s - same, d - 21)

# ShannonWiener
divTemp$ShannonWiener[row.nam$c125c]; divTemp$ShannonWiener[row.nam$s125c]; divTemp$ShannonWiener[row.nam$d125c] # 2.28 (s - same, d - 2.26)
diversity(table(full.125$cCIDr)); diversity(table(full.125$sCIDr)); diversity(table(full.125$dCIDr)) # 2.26 (s - 2.27, d - 2.23)

# Dominance
divTemp$Dominance[row.nam$c125c]; divTemp$Dominance[row.nam$s125c]; divTemp$Dominance[row.nam$d125c] # 0.134 (s - 0.133, d - 0.141)
1 - diversity(table(full.125$cCIDr), index = "simpson"); 1 - diversity(table(full.125$sCIDr), index = "simpson"); 1 - diversity(table(full.125$dCIDr), index = "simpson") # 0.137 (s - 0.138, d - 0.145)

# slide 150
divTemp$Richness[row.nam$c150c]; divTemp$Richness[row.nam$s150c]; divTemp$Richness[row.nam$d150c] # 20 (s - 18, d - 21)
specnumber(table(full.150$cCIDr)); specnumber(table(full.150$sCIDr)); specnumber(table(full.150$dCIDr)) # 20 (s - 19, d - 21)

# ShannonWiener
divTemp$ShannonWiener[row.nam$c150c]; divTemp$ShannonWiener[row.nam$s150c]; divTemp$ShannonWiener[row.nam$d150c] # 2.213 (s - 2.189, d - 2.303)
diversity(table(full.150$cCIDr)); diversity(table(full.150$sCIDr)); diversity(table(full.150$dCIDr)) # 2.213 (s - 2.206, d - 2.289)

# Dominance
divTemp$Dominance[row.nam$c150c]; divTemp$Dominance[row.nam$s150c]; divTemp$Dominance[row.nam$d150c] # 0.164 (s - 0.166, d - 0.151)
1 - diversity(table(full.150$cCIDr), index = "simpson"); 1 - diversity(table(full.150$sCIDr), index = "simpson"); 1 - diversity(table(full.150$dCIDr), index = "simpson") # 0.166 (s - 0.165, d - 0.156)


# 10e. Output the diversity file -------------------------------------------
tmp.div <- divTemp
tmp.div$ID <- paste(tmp.div$Analysis, tmp.div$Person)
head(tmp.div)
tmp.div <- reshape(tmp.div, direction = "wide", v.names = names(tmp.div)[!(names(tmp.div) %in% c("Analysis", "Person", "Size", "ID"))], timevar = "Size", idvar = "ID")
tmp.div <- tmp.div[, names(tmp.div) != "ID"]
tmp.div[, grepl("1", names(tmp.div)) & !grepl("Rich", names(tmp.div))] <- round(tmp.div[, grepl("1", names(tmp.div)) & !grepl("Rich", names(tmp.div))], 2)
write.csv(tmp.div, file = "ASOutputs/Diversity_temperature.csv", row.names = FALSE)
rm(tmp.div)

# 10f. Direction of change -------------------------------------------------
# add some columns to the diversity dataframe to investigate the direction of the change.
divTemp$Dir_D <- divTemp$Dir_SW <- divTemp$Dir_R <- divTemp$Dir_SST <- NA
divTemp[row.nam$s125, grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[row.nam$s125, c("SST10m", "Richness", "ShannonWiener", "Dominance")]), 2, as.numeric(divTemp[row.nam$c125c, c("SST10m", "Richness", "ShannonWiener", "Dominance")])) > 0, 1, -1)

divTemp[row.nam$s150, grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[row.nam$s150, c("SST10m", "Richness", "ShannonWiener", "Dominance")]), 2, as.numeric(divTemp[row.nam$c150c, c("SST10m", "Richness", "ShannonWiener", "Dominance")])) > 0, 1, -1)

divTemp[row.nam$d125, grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[row.nam$d125, c("SST10m", "Richness", "ShannonWiener", "Dominance")]), 2, as.numeric(divTemp[row.nam$c125c, c("SST10m", "Richness", "ShannonWiener", "Dominance")])) > 0, 1, -1)

divTemp[row.nam$d150, grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[row.nam$d150, c("SST10m", "Richness", "ShannonWiener", "Dominance")]), 2, as.numeric(divTemp[row.nam$c150c, c("SST10m", "Richness", "ShannonWiener", "Dominance")])) > 0, 1, -1)

# look at these summaries
table(divTemp$Dir_SST, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
table(divTemp$Dir_R, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
table(divTemp$Dir_SW, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
table(divTemp$Dir_D, paste(divTemp$Analysis, divTemp$Size, sep = "_"))

table(divTemp$Dir_R, divTemp$Dir_SW, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
table(divTemp$Dir_R, divTemp$Dir_D, paste(divTemp$Analysis, divTemp$Size, sep = "_"))


# 10g. Comparison with ForCenS -----------------------------------------------
# see how the datasets compare with the ForCenS values for that latitude
# only run with the >150 as that is what ForCenS is.
ForCenSred <- as.data.frame(read_excel("Data/Siccha_ForCenS.xlsx", sheet = "ForCenSred", na = "N/A"))
ForCenSred[, 22:62][is.na(ForCenSred[, 22:62])] <- 0
str(ForCenSred)

png("ASFigures/Div_cf_ForCenSred.png", 1000, 400)
par(mfrow = c(1,3), mar = c(5.1, 5.1, 3.1, 2.1))
# species richness
tmp.rich <- specnumber(ForCenSred[, 22:62]) # richness(ForCenSred[,22:62])
plot(ForCenSred$Latitude, tmp.rich, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Richness", cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(-75, 75))
points(rep(30.2, nrow(divTemp[divTemp$Size == 150,])), divTemp$Richness[divTemp$Size == 150], col = "blue", pch = 16, cex = 2)
points(30.2, divTemp$Richness[row.nam$c150c], col = "red", pch = 16, cex = 2)

# ShannonWiener
tmp.sw <- diversity(ForCenSred[,22:62])
plot(ForCenSred$Latitude, tmp.sw, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Shannon-Wiener", cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(-75, 75))
points(rep(30.2, nrow(divTemp[divTemp$Size == 150,])), divTemp$ShannonWiener[divTemp$Size == 150], col = "blue", pch = 16, cex = 1.8)
points(30.2, divTemp$ShannonWiener[row.nam$c150c], col = "red", pch = 16, cex = 1.8)

# Dominance
tmp.dom <- (1 - diversity(ForCenSred[,22:62], index = "simpson"))
plot(ForCenSred$Latitude, tmp.dom, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Dominance", cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(-75, 75))
points(rep(30.2, nrow(divTemp[divTemp$Size == 150,])), divTemp$Dominance[divTemp$Size == 150], col = "blue", pch = 16, cex = 1.8)
points(30.2, divTemp$Dominance[row.nam$c150c], col = "red", pch = 16, cex = 1.8)

par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("ASFigures/Div_cf_ForCenSred_Atl.png", 900, 360)
par(mfrow = c(1,3), mar = c(5.6, 5.1, 3.1, 2.1), mgp = c(3.5, 1, 0))
# species richness
plot(ForCenSred$Latitude[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], tmp.rich[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], pch = 16, col = "grey50", xlab = "Latitude", ylab = "Richness", cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(-65, 65))
points(rep(30.2, nrow(divTemp[divTemp$Size == 150,])), divTemp$Richness[divTemp$Size == 150], col = "blue", pch = 16, cex = 2)
points(30.2, divTemp$Richness[row.nam$c150c], col = "red", pch = 16, cex = 2)

# ShannonWiener
plot(ForCenSred$Latitude[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], tmp.sw[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], pch = 16, col = "grey50", xlab = "Latitude", ylab = "Shannon-Wiener Diversity", cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(-65, 65))
points(rep(30.2, nrow(divTemp[divTemp$Size == 150,])), divTemp$ShannonWiener[divTemp$Size == 150], col = "blue", pch = 16, cex = 1.8)
points(30.2, divTemp$ShannonWiener[row.nam$c150c], col = "red", pch = 16, cex = 1.8)

# Dominance
plot(ForCenSred$Latitude[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], tmp.dom[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], pch = 16, col = "grey50", xlab = "Latitude", ylab = "Dominance", cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(-65, 65))
points(rep(30.2, nrow(divTemp[divTemp$Size == 150,])), divTemp$Dominance[divTemp$Size == 150], col = "blue", pch = 16, cex = 1.8)
points(30.2, divTemp$Dominance[row.nam$c150c], col = "red", pch = 16, cex = 1.8)

par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))
dev.off()

# signal of diversity after richness
png("ASFigures/Div_cf_ForCenSred_Atl_rich.png", 900, 360)
par(mfrow = c(1,3), mar = c(5.6, 5.1, 3.1, 2.1), mgp = c(3.5, 1, 0))
# species richness
plot(ForCenSred$Latitude[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], tmp.rich[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], pch = 16, col = matlab.like(30)[tmp.rich[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11]], xlab = "Latitude", ylab = "Richness", cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(-65, 65))
points(rep(30.2, nrow(divTemp[divTemp$Size == 150,])), divTemp$Richness[divTemp$Size == 150], bg = matlab.like(30)[divTemp$Richness[divTemp$Size == 150]], pch = 21, cex = 2)
points(30.2, divTemp$Richness[row.nam$c150c], bg = matlab.like(30)[divTemp$Richness[row.nam$c150c]], pch = 21, cex = 3)

# ShannonWiener
plot(ForCenSred$Latitude[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], tmp.sw[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], pch = 16, col = matlab.like(30)[tmp.rich[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11]], xlab = "Latitude", ylab = "Shannon-Wiener Diversity", cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(-65, 65))
points(rep(30.2, nrow(divTemp[divTemp$Size == 150,])), divTemp$ShannonWiener[divTemp$Size == 150], bg = matlab.like(30)[divTemp$Richness[divTemp$Size == 150]], pch = 21, cex = 2)
points(30.2, divTemp$ShannonWiener[row.nam$c150c], bg = matlab.like(30)[divTemp$Richness[row.nam$c150c]], pch = 21, cex = 3)

# Dominance
plot(ForCenSred$Latitude[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], tmp.dom[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11], pch = 16, col = matlab.like(30)[tmp.rich[ForCenSred$Ocean == 7|ForCenSred$Ocean == 11]], xlab = "Latitude", ylab = "Dominance", cex.lab = 2, cex.axis = 1.5, las = 1, xlim = c(-65, 65), log = "y")
points(rep(30.2, nrow(divTemp[divTemp$Size == 150,])), divTemp$Dominance[divTemp$Size == 150], bg = matlab.like(30)[divTemp$Richness[divTemp$Size == 150]], pch = 21, cex = 2)
points(30.2, divTemp$Dominance[row.nam$c150c], bg = matlab.like(30)[divTemp$Richness[row.nam$c150c]], pch = 21, cex = 3)

par(mfrow = c(1,1), mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3, 1, 0))
dev.off()

rm(tmp.sw, tmp.dom, tmp.rich, ForCenSred)

# 10h. Simulations for studying changes -----------------------------------
tmp <- data.frame("a" = c(10, 10, 0, 20, 12, 10, 10, 5, 9, 10, 8, 8, 10, 11), "b" = c(10, 10, 10, 0, 10, 0, 10, 10, 10, 10, 12, 10, 10, 10), "c" = c(2, 2, 2, 2, 0, 12, 4, 2, 2, 1, 2, 4, 3, 1), "d" = c(2, 0, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 1, 2), "e" = c(0, 2, 10, 0, 0, 0, 0, 5, 1, 1, 0, 0, 0, 0))
rownames(tmp) <- c("Orig", "ChangeSpListR", "ChangeSpListC", "LumpCc", "LumpCr", "LumpRc", "LumpRr", "SplitCe", "SplitCu", "SplitR", "BoundShiftCc", "BoundShift2Rc", "BoundShift3Rr", "BoundShift4Cr")
tmp$Rich <- specnumber(tmp[, 1:5]) # richness
tmp$SW <- diversity(tmp[, 1:5]) # Shannon-Wiener
tmp$Dom <- (1 - diversity(tmp[, 1:5], index = "simpson")) # dominance
tmp
rm(tmp)


# 10i. Without incomplete data --------------------------------------------
# s125
summary(divTemp$Richness[row.nam$s125])
summary(divTemp$Richness[divTemp$Person %in% accuracy$PersonID[accuracy$ptID125 == 100] & divTemp$Analysis == "Slide"& divTemp$Size == 125])

# s150
summary(divTemp$Richness[row.nam$s150])
summary(divTemp$Richness[divTemp$Person %in% accuracy$PersonID[accuracy$ptID150 == 100] & divTemp$Analysis == "Slide"& divTemp$Size == 150])

# d125
summary(divTemp$Richness[row.nam$d125])
summary(divTemp$Richness[divTemp$Person %in% accuracy$PersonID[accuracy$ptID125 == 100] & divTemp$Analysis == "Digital"& divTemp$Size == 125])

# d150
summary(divTemp$Richness[row.nam$d150])
summary(divTemp$Richness[divTemp$Person %in% accuracy$PersonID[accuracy$ptID150 == 100] & divTemp$Analysis == "Digital"& divTemp$Size == 150])

# 11. Add person level characteristics to accuracy -------------------------

# 11a. How does agreement depend on whether workers routinely count specimens? --------
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

# 11b. Adding percentage of specimens identified ---------------------------
accuracy$ptID125 <- apply(full.125[,col.nam$c125], 2, function(x) sum(x != "na") / 300 * 100)
accuracy$ptID150 <- apply(full.150[,col.nam$c150], 2, function(x) sum(x != "na") / 300 * 100)[accuracy$PersonID]

# 11c. Adding Experience ---------------------------------------------------
accuracy$Experience <- people$ExperienceSlideA[match(accuracy$PersonID, people$SlideID)]
accuracy$Experience[accuracy$PersonID %in% c("1a", "2a")] <- people$ExperienceSlideA[1:2]
accuracy$Experience[accuracy$PersonID %in% c("1b", "2b")] <- people$ExperienceSlideB[1:2]
accuracy$Experience[accuracy$Analysis == "Digital"] <- people$ExperienceDigital[match(accuracy$PersonID[accuracy$Analysis == "Digital"], people$DigitalID)]

# 11d. Create a .csv file for the accuracy data ----------------------------
# create a subset of the data for this
tmp.sub <- accuracy[, c("PersonID", "Analysis", "Experience", "Routine", "cPtA125", "cPtA150", "sPtA125", "sPtA150", "dPtA125", "dPtA150", "ptID125", "ptID150")]
tmp.sub[, grep("125|150", names(tmp.sub))] <- round(tmp.sub[, grep("125|150", names(tmp.sub))], 2)
write.csv(tmp.sub, "ASOutputs/Accuracy.csv", row.names = FALSE)
rm(tmp.sub)

# 12. Species level accuracy ------------------------------------------------
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

# 13. Influence of size ----------------------------------------------------
# 13a. Sizes comparisons between >150 from >125 ------------------------------
size125[which(size125$Length < 150),]
summary(size125$cAgreement[which(size125$Length < 150)])
summary(size125$cAgreement[which(size125$Length > 150)])
summary(size150$cAgreement)

tmp <- data.frame(size = seq(130, 710, by = 10))
tmp$c125[tmp$size %in% names(table(round(size125$Length, -1)))] <- table(round(size125$Length, -1))
tmp$c150[tmp$size %in% names(table(round(size150$Length, -1)))] <- table(round(size150$Length, -1))

png("ASFigures/SizeComparison.png", 800, 500)
barplot(t(tmp[, 2:3]), beside = TRUE, legend.text = c("125", "150"), names.arg = tmp$size, las = 2)
dev.off()
# 13b. Size vs. maximum agreement ------------------------------------------
# when slide and digital consensus are calculated seperately
png("ASFigures/Size_agreement_sd125.png", 400, 480)
with(size125, plot(sAgreement, Length, pch = 16, main = "> 125", las = 1, xlab = "Agreement", ylab = expression(paste("Maximum diameter / ", mu, "m"))))
lines(names(tapply(size125$Length, size125$sAgreement, max)), tapply(size125$Length, size125$sAgreement, max), pch = 16)
with(size125, points(dAgreement, Length, pch = 16, col = "blue"))
lines(names(tapply(size125$Length, size125$dAgreement, max)), tapply(size125$Length, size125$dAgreement, max), pch = 16, col = 4)
legend("topleft", col = c(1, 4), pch = 16, legend = c("Slide", "Digital"))
dev.off()

png("ASFigures/Size_agreement_sd150.png", 400, 480)
with(size150, plot(sAgreement, Length, pch = 16, main = "> 150", las = 1, xlab = "Agreement", ylab = expression(paste("Maximum diameter / ", mu, "m"))))
lines(names(tapply(size150$Length, size150$sAgreement, max)), tapply(size150$Length, size150$sAgreement, max), pch = 16)
with(size150, points(dAgreement, Length, pch = 16, col = "blue"))
lines(names(tapply(size150$Length, size150$dAgreement, max)), tapply(size150$Length, size150$dAgreement, max), pch = 16, col = 4)
legend("topleft", col = c(1, 4), pch = 16, legend = c("Slide", "Digital"))
dev.off()

# with all the combined data
png("ASFigures/Size_agreement_c125.png")
with(size125, plot(Length, cAgreement, pch = 16, main = "> 125", las = 1, ylab = "Agreement", xlab = expression(paste("Maximum diameter / ", mu, "m"))))
lines(tapply(size125$Length, size125$cAgreement, max), names(tapply(size125$Length, size125$cAgreement, max)), pch = 16)
dev.off()

png("ASFigures/Size_agreement_c150.png")
with(size150, plot(Length, cAgreement, pch = 16, main = "> 150", las = 1, ylab = "Agreement", xlab = expression(paste("Maximum diameter / ", mu, "m"))))
lines(tapply(size150$Length, size150$cAgreement, max), names(tapply(size150$Length, size150$cAgreement, max)), pch = 16)
dev.off()

# the combined consensus, but slide / digital calculated separately
size125$scAgreement <- rowSums(full.125[,col.nam$s125] == full.125$cCID)/length(col.nam$s125)*100
size125$dcAgreement <- rowSums(full.125[,col.nam$d125] == full.125$cCID)/length(col.nam$d125)*100
size150$scAgreement <- rowSums(full.150[,col.nam$s150] == full.150$cCID)/length(col.nam$s150)*100
size150$dcAgreement <- rowSums(full.150[,col.nam$d150] == full.150$cCID)/length(col.nam$d150)*100

png("ASFigures/Size_agreement_sdc125.png", 400, 480)
with(size125, plot(Length, scAgreement, pch = 16, main = expression(paste(">125 ", mu, "m")), las = 1, ylab = "Percentage Agreement", xlab = expression(paste("Maximum diameter / ", mu, "m")), ylim = c(0, 100), cex.main = 1.5, cex.lab = 1.3))
lines(tapply(size125$Length, size125$scAgreement, max), names(tapply(size125$Length, size125$scAgreement, max)), pch = 16)
with(size125, points(Length, dcAgreement, pch = 16, col = "blue"))
lines(tapply(size125$Length, size125$dcAgreement, max), names(tapply(size125$Length, size125$dcAgreement, max)), pch = 16, col = 4)
legend("bottomright", col = c(1, 4), pch = 16, legend = c("Slide", "Digital"), cex = 1.2, pt.cex = 1.3)
dev.off()

png("ASFigures/Size_agreement_sdc150.png", 400, 480)
with(size150, plot(Length, scAgreement, pch = 16, main = expression(paste(">150 ", mu, "m")), las = 1, ylab = "Percentage Agreement", xlab = expression(paste("Maximum diameter / ", mu, "m")), ylim = c(0, 100), cex.main = 1.5, cex.lab = 1.3))
lines(tapply(size150$Length, size150$scAgreement, max), names(tapply(size150$Length, size150$scAgreement, max)), pch = 16)
with(size150, points(Length, dcAgreement, pch = 16, col = "blue"))
lines(tapply(size150$Length, size150$dcAgreement, max), names(tapply(size150$Length, size150$dcAgreement, max)), pch = 16, col = 4)
legend("bottomright", col = c(1, 4), pch = 16, legend = c("Slide", "Digital"), cex = 1.2, pt.cex = 1.3)
dev.off()


# 14. Outliers -------------------------------------------------------------
# Generate a dataframe of outliers.
# 14a. Outliers for each analysis ------------------------------------------
outliers <- data.frame(PersonID = accuracy$PersonID[1:26], Analysis = accuracy$Analysis[1:26])

# rank for MDS 125
outliers$MDS125 <- NA
tmp.pt <- nmds$CC125f$points["cCID", ]
tmp.tab <- sqrt((nmds$CC125f$points[, 1] - tmp.pt[1])^2 + (nmds$CC125f$points[, 2] - tmp.pt[2])^2)
outliers$MDS125[outliers$Analysis == "Slide"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))]))
outliers$MDS125[outliers$Analysis == "Digital"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))]))
rm(tmp.pt, tmp.tab)

# and mds 150
outliers$MDS150 <- NA
tmp.pt <- nmds$CC150f$points["cCID", ]
tmp.tab <- sqrt((nmds$CC150f$points[,1] - tmp.pt[1])^2 + (nmds$CC150f$points[,2] - tmp.pt[2])^2)
outliers$MDS150[outliers$Analysis == "Slide"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))]))
outliers$MDS150[outliers$Analysis == "Digital"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))]))
rm(tmp.pt, tmp.tab)

# for percentage accuracy
outliers$PtA125 <- NA
outliers$PtA125[outliers$Analysis == "Slide"][order(100-accuracy$cPtA125[match(outliers$PersonID[outliers$Analysis == "Slide"], accuracy$PersonID)])] <- sort(rank(100-accuracy$cPtA125[accuracy$Analysis == "Slide"]))
outliers$PtA125[outliers$Analysis == "Digital"][order(100-accuracy$cPtA125[match(outliers$PersonID[outliers$Analysis == "Digital"], accuracy$PersonID)])] <- sort(rank(100-accuracy$cPtA125[accuracy$Analysis == "Digital"]))

outliers$PtA150 <- NA
outliers$PtA150[outliers$Analysis == "Slide"][order(100-accuracy$cPtA150[match(outliers$PersonID[outliers$Analysis == "Slide"], accuracy$PersonID)])] <- sort(rank(100-accuracy$cPtA150[accuracy$Analysis == "Slide"]))
outliers$PtA150[outliers$Analysis == "Digital"][order(100-accuracy$cPtA150[match(outliers$PersonID[outliers$Analysis == "Digital"], accuracy$PersonID)])] <- sort(rank(100-accuracy$cPtA150[accuracy$Analysis == "Digital"]))

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
outliers$Richness125 <- NA
tmp.pt <- abs(divTemp$Richness[row.nam$s125] - divTemp$Richness[row.nam$c125c])
names(tmp.pt) <- divTemp$Person[row.nam$s125]
outliers$Richness125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$Richness[row.nam$d125] - divTemp$Richness[row.nam$c125c])
names(tmp.pt) <- divTemp$Person[row.nam$d125]
outliers$Richness125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$Richness150 <- NA
tmp.pt <- abs(divTemp$Richness[row.nam$s150] - divTemp$Richness[row.nam$c150c])
names(tmp.pt) <- divTemp$Person[row.nam$s150]
outliers$Richness150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$Richness[row.nam$d150] - divTemp$Richness[row.nam$c150c])
names(tmp.pt) <- divTemp$Person[row.nam$d150]
outliers$Richness150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# Dominance
outliers$Dominance125 <- NA
tmp.pt <- abs(divTemp$Dominance[row.nam$s125] - divTemp$Dominance[row.nam$c125c])
names(tmp.pt) <- divTemp$Person[row.nam$s125]
outliers$Dominance125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$Dominance[row.nam$d125] - divTemp$Dominance[row.nam$c125c])
names(tmp.pt) <- divTemp$Person[row.nam$d125]
outliers$Dominance125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$Dominance150 <- NA
tmp.pt <- abs(divTemp$Dominance[row.nam$s150] - divTemp$Dominance[row.nam$c150c])
names(tmp.pt) <- divTemp$Person[row.nam$s150]
outliers$Dominance150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$Dominance[row.nam$d150] - divTemp$Dominance[row.nam$c150c])
names(tmp.pt) <- divTemp$Person[row.nam$d150]
outliers$Dominance150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# ShannonWiener
outliers$ShannonWiener125 <- NA
tmp.pt <- abs(divTemp$ShannonWiener[row.nam$s125] - divTemp$ShannonWiener[row.nam$c125c])
names(tmp.pt) <- divTemp$Person[row.nam$s125]
outliers$ShannonWiener125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$ShannonWiener[row.nam$d125] - divTemp$ShannonWiener[row.nam$c125c])
names(tmp.pt) <- divTemp$Person[row.nam$d125]
outliers$ShannonWiener125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$ShannonWiener150 <- NA
tmp.pt <- abs(divTemp$ShannonWiener[row.nam$s150] - divTemp$ShannonWiener[row.nam$c150c])
names(tmp.pt) <- divTemp$Person[row.nam$s150]
outliers$ShannonWiener150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$ShannonWiener[row.nam$d150] - divTemp$ShannonWiener[row.nam$c150c])
names(tmp.pt) <- divTemp$Person[row.nam$d150]
outliers$ShannonWiener150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# weighted sums
outliers$fullSum <- rowSums(outliers[, 3:13])

outliers$wtSum <- rowSums(outliers[, grep("^MDS", names(outliers))])/2 + rowSums(outliers[, grep("^Pt", names(outliers))])/2 + outliers$SST + rowSums(outliers[, grep("^Richness", names(outliers))])/6 + rowSums(outliers[, grep("^Dominance", names(outliers))])/6 + rowSums(outliers[, grep("^ShannonWiener", names(outliers))])/6

outliers$wtSumPt <- outliers$wtSum / c(rep(4*17, 17), rep(4*9, 9)) * 100


# 14b. outliers for digital and slide combined -----------------------------
# rank for MDS 125
outliers$cMDS125 <- NA
tmp.pt <- nmds$CC125f$points["cCID", ]
tmp.tab <- sqrt((nmds$CC125f$points[row.names(nmds$CC125f$points),1] - tmp.pt[1])^2 + (nmds$CC125f$points[,2] - tmp.pt[2])^2)
outliers$cMDS125[order(tmp.tab[match(outliers$PersonID, names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID, names(tmp.tab))]))
rm(tmp.pt, tmp.tab)

# and mds 150
outliers$cMDS150 <- NA
tmp.pt <- nmds$CC150f$points["cCID", ]
tmp.tab <- sqrt((nmds$CC150f$points[,1] - tmp.pt[1])^2 + (nmds$CC150f$points[,2] - tmp.pt[2])^2)
outliers$cMDS150[order(tmp.tab[match(outliers$PersonID, names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID, names(tmp.tab))]))
rm(tmp.pt, tmp.tab)

# for percentage accuracy
outliers$cPtA125 <- NA
outliers$cPtA125[order(100-accuracy$cPtA125[match(outliers$PersonID, accuracy$PersonID)])] <- sort(rank(100-accuracy$cPtA125))

outliers$cPtA150 <- NA
outliers$cPtA150[order(100-accuracy$cPtA150[match(outliers$PersonID, accuracy$PersonID)])] <- sort(rank(100-accuracy$cPtA150))

# SST
outliers$cSST <- NA
tmp.pt <- abs(divTemp$SST10m[c(row.nam$s150, row.nam$d150)] - 21.76)
names(tmp.pt) <- divTemp$Person[c(row.nam$s150, row.nam$d150)]
outliers$cSST[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# richness
outliers$cRichness125 <- NA
tmp.pt <- abs(divTemp$Richness[c(row.nam$s125, row.nam$d125)] - divTemp$Richness[row.nam$c125c])
names(tmp.pt) <- divTemp$Person[c(row.nam$s125, row.nam$d125)]
outliers$cRichness125[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$cRichness150 <- NA
tmp.pt <- abs(divTemp$Richness[c(row.nam$s150, row.nam$d150)] - divTemp$Richness[row.nam$c150c])
names(tmp.pt) <- divTemp$Person[c(row.nam$s150, row.nam$d150)]
outliers$cRichness150[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# Dominance
outliers$cDominance125 <- NA
tmp.pt <- abs(divTemp$Dominance[c(row.nam$s125, row.nam$d125)] - divTemp$Dominance[row.nam$c125c])
names(tmp.pt) <- divTemp$Person[c(row.nam$s125, row.nam$d125)]
outliers$cDominance125[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$cDominance150 <- NA
tmp.pt <- abs(divTemp$Dominance[c(row.nam$s150, row.nam$d150)] - divTemp$Dominance[row.nam$c150c])
names(tmp.pt) <- divTemp$Person[c(row.nam$s150, row.nam$d150)]
outliers$cDominance150[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# ShannonWiener
outliers$cShannonWiener125 <- NA
tmp.pt <- abs(divTemp$ShannonWiener[c(row.nam$s125, row.nam$d125)] - divTemp$ShannonWiener[row.nam$c125c])
names(tmp.pt) <- divTemp$Person[c(row.nam$s125, row.nam$d125)]
outliers$cShannonWiener125[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$cShannonWiener150 <- NA
tmp.pt <- abs(divTemp$ShannonWiener[c(row.nam$s150, row.nam$d150)] - divTemp$ShannonWiener[row.nam$c150c])
names(tmp.pt) <- divTemp$Person[c(row.nam$s150, row.nam$d150)]
outliers$cShannonWiener150[order(tmp.pt[match(outliers$PersonID, names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# weighted sums
outliers$cfullSum <- rowSums(outliers[, grep("^c", names(outliers))])

outliers$cwtSum <- rowSums(outliers[, grep("cMDS", names(outliers))])/2 + rowSums(outliers[, grep("^c.*Pt", names(outliers))])/2 + outliers$cSST + rowSums(outliers[, grep("cRichness", names(outliers))])/6 + rowSums(outliers[, grep("cDominance", names(outliers))])/6 + rowSums(outliers[, grep("cShannonWiener", names(outliers))])/6

outliers$cwtSumPt <- outliers$cwtSum / (4*26) * 100


# 14c. Consider how to interpret these -------------------------------------
# looking for correlations
pairs(outliers[, grep("^c", names(outliers))])
pairs(outliers[, grep("^c.*125", names(outliers))])
pairs(outliers[, grep("^c.*150", names(outliers))])

# compare the outliers with person level information
# create a new dataframe that combines the person level information with the outlier information
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

# look at the plots for the slide / digital analyses separately
#par(ask = TRUE)
par(mfrow = c(4, 4))
par(mar = c(2,2,1,1))
for (i in c("Residence", "School", "Experience", "Routine", "Region", "ptID125", "ptID150")) {
  for (j in names(tmp.out)[3:16]) 
    plot(factor(tmp.out[tmp.out$Analysis == "Slide",i]), tmp.out[tmp.out$Analysis == "Slide",j], main = i, col = "red")
  plot(1,1, type = "n", axes = FALSE)
  plot(1,1, type = "n", axes = FALSE)
  for (j in names(tmp.out)[3:16])
    plot(factor(tmp.out[tmp.out$Analysis == "Digital",i]), tmp.out[tmp.out$Analysis == "Digital",j], main = j, col = "blue")
  plot(1,1, type = "n", axes = FALSE)
  plot(1,1, type = "n", axes = FALSE)
}
rm(i, j)
#par(ask = FALSE)

# look at the plots for the combined outliers (both slide and digital)
#par(ask = TRUE)
par(mfrow = c(4, 4))
par(mar = c(2,2,1,1))
for (i in c("Residence", "School", "Experience", "Routine", "Region", "ptID125", "ptID150")) {
  for (j in grep("^c", names(tmp.out), value = TRUE))
    plot(factor(tmp.out[,i]), tmp.out[,j], main = j, col = "red")
  plot(1,1, type = "n", axes = FALSE)
  plot(1,1, type = "n", axes = FALSE)
}
rm(i, j)
#par(ask = FALSE)
par(mfrow = c(1,1))
par(mar = c(5.1, 4.1, 4.1, 2.1))

plot(tmp.out$Analysis, tmp.out$cwtSumPt)

rm(tmp.s, tmp.out, tmp.d, tmp)

# 15. Save the data -------------------------------------------------------
save.image("ASOutputs/Reanalysis_NA_IF.RData")
# output the full species names
tmp.125<- full.125
tmp.150 <- full.150
for(i in 1:nrow(sp.abb)) {
  tmp.125 <- as.data.frame(apply(tmp.125, 2, function(x) gsub(paste("^", sp.abb$Abbreviation[i], "$", sep = ""), sp.abb$Species[i], x)))
  tmp.150 <- as.data.frame(apply(tmp.150, 2, function(x) gsub(paste("^", sp.abb$Abbreviation[i], "$", sep = ""), sp.abb$Species[i], x)))
}
write.csv(tmp.125, "ASOutputs/PersonIDs_125.csv", row.names = FALSE)
write.csv(tmp.150, "ASOutputs/PersonIDs_150.csv", row.names = FALSE)
rm(tmp.125, tmp.150, i)

# write out the confusion matrices
lapply(seq_along(conf.mat), function(i) write.csv(conf.mat[[i]]$table, paste("ASOutputs/Confusion matrices/", names(conf.mat)[i], ".csv", sep = "")))

