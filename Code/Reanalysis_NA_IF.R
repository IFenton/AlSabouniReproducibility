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
library(readxl)
library(caret) # for the confusion matrix
library(colorRamps) # colours
library(stringr) # for confusion matrix axes
library(vegan)
library(cluster) # for the dendrogram
library(RColorBrewer)
source("../../../../Code/Confusion_matrix.R")


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
people <- as.data.frame(read_excel("Data/PeopleMetadata.xlsx", na = "NA"))

# specimen size
size125 <- as.data.frame(read_excel("Data/SpecimenSize.xlsx", sheet = "Size125"))
size150 <- as.data.frame(read_excel("Data/SpecimenSize.xlsx", sheet = "Size150"))

# diversity / temperature 
divTemp <- as.data.frame(read_excel("Data/DiversityTemp.xlsx", na = "NA"))

# 1b. Sanity check --------------------------------------------------------
str(slide125)

# check that all the levels are correct and match to those in the abbreviations
for (i in datasets) {
  tmp <- sort(unique(unlist(sapply(i[, 2:ncol(i)], unique))))
  print(tmp %in% sp.abb$Abbreviation) 
}
rm(i, datasets, tmp)
# all of those are true, so all the abbreviations match to real species

# get the digital pair and the slide pair are the same size
dim(digital125) == dim(digital150)
dim(slide125) == dim(slide150)

# create a list of the columns containing data
col.nam <- list(s125 = which(nchar(names(slide125)) < 5), s150 = which(nchar(names(slide125)) < 5), d125 = which(nchar(names(digital125)) < 5), d150 = which(nchar(names(digital125)) < 5))

# 2. Calculate the consensus ----------------------------------------------

# 2a. Consensus 50 --------------------------------------------------------
# calculate the cutoff for consensus 50
# as both are odd numbres, by rounding down, the cutoff anything greater than this value.  
c50_cutoff <- list()
c50_cutoff$slide <- round(ncol(slide125[, col.nam$s125])/2, 0)
c50_cutoff$digital <- round(ncol(digital125[, col.nam$d125])/2, 0)

# calculate the consensus values
slide125$IFc50 <- apply(slide125[, col.nam$s125], 1, function (x) ifelse(length(which(table(x) > c50_cutoff$slide)) > 0, names(table(x))[which(table(x) > c50_cutoff$slide)], "nc")) 
# if one name has the majority (i.e. the length of the table > 8 is 1), then return that name
# if the there are no names that have more than 8 (given there are 17 IDs, then consensus-50 needs at least 9 to match), then return "nc"
slide150$IFc50 <- apply(slide150[, col.nam$s150], 1, function (x) ifelse(length(which(table(x) > c50_cutoff$slide)) > 0, names(table(x))[which(table(x) > c50_cutoff$slide)], "nc")) 
digital125$IFc50 <- apply(digital125[, col.nam$d125], 1, function (x) ifelse(length(which(table(x) > c50_cutoff$digital)) > 0, names(table(x))[which(table(x) > c50_cutoff$digital)], "nc")) 
digital150$IFc50 <- apply(digital150[, col.nam$d150], 1, function (x) ifelse(length(which(table(x) > c50_cutoff$digital)) > 0, names(table(x))[which(table(x) > c50_cutoff$digital)], "nc")) 


# see if this matches Nadia's results
write.csv(slide125[which(slide125$consensus50 != slide125$IFc50), ], "Outputs/Slide125Cmismatch.csv", row.names = FALSE)
slide150[which(slide150$consensus50 != slide150$IFc50), ] # errors are in Nadia's
digital125[which(digital125$consensus50 != digital125$IFc50), ] # 2 errors in Nadia's results
digital150[which(digital150$consensus50 != digital150$IFc50), ] # matches for everything

# does removing those that are NAs help?
tmp <- apply(slide125[, col.nam$s125], 1, function (x) ifelse(length(which(table(x) > (c50_cutoff$slide - sum(x == "na")))) > 0, names(table(x))[which(table(x) > (c50_cutoff$slide - sum(x == "na")))], "nc")) 
write.csv(slide125[which(slide125$consensus50 != tmp), ], "Outputs/Slide125CmismatchNA.csv", row.names = FALSE)
# no, it makes it worse
rm(tmp)

# 2b. Consensus 20 --------------------------------------------------------
# check whether 20% is actually the maximum for the consensus (so what is the highest count per specimen) ignoring na's
slide125$IFmaxCon <- apply(slide125[, col.nam$s125], 1, function (x)  max(table(x[x != 'na']))) 
slide150$IFmaxCon <- apply(slide150[, col.nam$s150], 1, function (x) max(table(x[x != 'na']))) 
digital125$IFmaxCon <- apply(digital125[, col.nam$d125], 1, function (x) max(table(x[x != 'na']))) 
digital150$IFmaxCon <- apply(digital150[, col.nam$d150], 1, function (x) max(table(x[x != 'na']))) 

# look at the summaries of these
table(slide125$IFmaxCon) # so actually for the slides, it should be consensus-18 (3/17)
table(slide150$IFmaxCon)
table(digital125$IFmaxCon) # and here it would be consensus-22 (2/9)
table(digital150$IFmaxCon)

# calculate consensus-20 
# I'm actually calculating this a the minimum consensus. So take most frequent name - if there are multiple take the first alphabetically.
# where the maximum is 'na', Nadia has ignored that
# I'll see if that is what Nadia has done.
slide125$IFcMin <- apply(slide125[, col.nam$s125], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][1]) 
slide150$IFcMin <- apply(slide150[, col.nam$s150], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][1])
digital125$IFcMin <- apply(digital125[, col.nam$d125], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][1])
digital150$IFcMin <- apply(digital150[, col.nam$d150], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][1])

# see if this matches Nadia's results
slide125[which(slide125$consensus20 != slide125$IFcMin), ] # I think they are errors
slide150[which(slide150$consensus20 != slide150$IFcMin), ] # I think they are errors
digital125[which(digital125$consensus20 != digital125$IFcMin), ] # pachS vs rub
digital150[which(digital150$consensus20 != digital150$IFcMin), ] # I think these were both errors (assuming ruber pink comes after ruber)

# by alphabetical, was that based on abbreviations or on the full names. 
slide150[order(slide150$IFmaxCon),][1,]
apply(slide150[146, col.nam$s150], 1, function (x) table(x))
sp.abb[sp.abb$Abbreviation %in% c("cal", "falc"),] # if it were on full names, then consensus name would be falconensis, but it is calida
# It was based on the abbreviations (which is what I am doing)

# was the alphabet only used when there wasn't a single maximum or was it used for all those greater than 20%
tmp <- apply(slide125[, col.nam$s125], 1, function (x) names(table(x))[which(table(x) > 2)][1])
sum(tmp != slide125$consensus20) # 71 mismatches
sum(slide125$consensus20 != slide125$IFcMin) # 14 mismatches
rm(tmp)
# therefore I think the alphabet was only used when there wasn't a single maximum.

# 2c. Generate a species total dataframe ----------------------------------
# create a blank data frame
slide125sp <- data.frame(species = c(sp.abb$Abbreviation[1:2], sort(sp.abb$Abbreviation[3:nrow(sp.abb)]))) 
slide125sp[, names(slide125)[2:ncol(slide125)]] <- NA
head(slide125sp)
# fill that dataframe
for (i in 2:ncol(slide125)) {
  # table each column
  tmp <- table(slide125[,i]) 
  # add them in in the right order
  slide125sp[, names(slide125)[i] == names(slide125sp)] <- tmp[match(slide125sp$species, names(tmp))]
}
slide125sp[is.na(slide125sp)] <- 0
rm(i, tmp)

# repeat for the other datasets
# species 150
slide150sp <- data.frame(species = c(sp.abb$Abbreviation[1:2], sort(sp.abb$Abbreviation[3:nrow(sp.abb)])))  
slide150sp[, names(slide150)[2:ncol(slide150)]] <- NA
# fill that dataframe
for (i in 2:ncol(slide150)) {
  tmp <- table(slide150[,i]) 
  slide150sp[, names(slide150)[i] == names(slide150sp)] <- tmp[match(slide150sp$species, names(tmp))]
}
slide150sp[is.na(slide150sp)] <- 0
rm(i, tmp)
# digital 125
digital125sp <- data.frame(species = c(sp.abb$Abbreviation[1:2], sort(sp.abb$Abbreviation[3:nrow(sp.abb)]))) 
digital125sp[, names(digital125)[2:ncol(digital125)]] <- NA
# fill that dataframe
for (i in 2:ncol(digital125)) {
  tmp <- table(digital125[,i]) 
  digital125sp[, names(digital125)[i] == names(digital125sp)] <- tmp[match(digital125sp$species, names(tmp))]
}
digital125sp[is.na(digital125sp)] <- 0
rm(i, tmp)
# digital 150
digital150sp <- data.frame(species = c(sp.abb$Abbreviation[1:2], sort(sp.abb$Abbreviation[3:nrow(sp.abb)]))) 
digital150sp[, names(digital150)[2:ncol(digital150)]] <- NA
# fill that dataframe
for (i in 2:ncol(digital150)) {
  tmp <- table(digital150[,i]) 
  digital150sp[, names(digital150)[i] == names(digital150sp)] <- tmp[match(digital150sp$species, names(tmp))]
}
digital150sp[is.na(digital150sp)] <- 0
rm(i, tmp)

# write these out to compare with Nadia's data (supplementary table 1)
write.csv(slide125sp, file = "Outputs/Supp1slide125.csv", row.names = FALSE)
write.csv(slide150sp, file = "Outputs/Supp1slide150.csv", row.names = FALSE)
write.csv(digital125sp, file = "Outputs/Supp1digital125.csv", row.names = FALSE)
write.csv(digital150sp, file = "Outputs/Supp1digital150.csv", row.names = FALSE)
# they all agree (with her data), though that is not really suprising given her data uses formulas from the same raw data. 


# 2d. Generate comparison tables ------------------------------------------
# Supp Table 3, add these to the size data
head(size125)

size125[, c("slideCon1", "slideCon2", "slideCon3", "slideCon4", "slideAgreement", "digitalCon1", "digitalCon2", "digitalCon3", "digitalCon4", "digitalAgreement")] <- NA
size150[, c("slideCon1", "slideCon2", "slideAgreement", "digitalCon1", "digitalCon2", "digitalCon3", "digitalAgreement")] <- NA


for (i in 1:4) {
  # add consensus values for 125
  size125[, grep("slideCon", names(size125))[i]] <- apply(slide125[, col.nam$s125], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][i]) 
  size125[, grep("digitalCon", names(size125))[i]] <- apply(digital125[, col.nam$d125], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][i])
  # and the consensus values for 150
  if (length(grep("slideCon", names(size150))) >= i)  # there aren't as many choices for the 150 consensus
    size150[, grep("slideCon", names(size150))[i]] <- apply(slide150[, col.nam$s150], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][i])
  if (length(grep("digitalCon", names(size150))) >= i)
    size150[, grep("digitalCon", names(size150))[i]] <- apply(digital150[, col.nam$d150], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][i])
}
rm(i)

# add in the specimen agreement values
size125$slideAgreement <- slide125$IFmaxCon/17*100
size150$slideAgreement <- slide150$IFmaxCon/17*100
size125$digitalAgreement <- digital125$IFmaxCon/9*100
size150$digitalAgreement <- digital150$IFmaxCon/9*100

# what fraction of specimens don't have a single consensus
sum(!is.na(size125$slideCon2)) # 16
sum(!is.na(size125$digitalCon2)) # 34
sum(!is.na(size150$slideCon2)) # 11
sum(!is.na(size150$digitalCon2)) # 18

# have three possibilities
sum(!is.na(size125$slideCon3)) # 1
sum(!is.na(size125$digitalCon3)) # 4
# slide 150 0
sum(!is.na(size150$digitalCon3)) # 0

# have four possibilities
sum(!is.na(size125$slideCon4)) # 1
sum(!is.na(size125$digitalCon4)) # 1
# slide 150 0
# digital 150 0

# save it to compare with Nadia's
write.csv(size125, file = "Outputs/Supp3_125.csv", row.names = FALSE)
write.csv(size150, file = "Outputs/Supp3_150.csv", row.names = FALSE)


# 3. Agreement between workers (pairwise comparisons) ---------------------

# 3a. C20_score_participant -----------------------------------------------
accuracySlide <- data.frame(PersonID = names(slide125)[col.nam$s125][c(1, 3, 2, 4:length(names(slide125)[col.nam$s125]))], stringsAsFactors = FALSE)
accuracyDigital <- data.frame(PersonID = names(digital125)[col.nam$d125], stringsAsFactors = FALSE)

# what is the percentage accuracy (i.e. how good is the match to the consensus)
# based on Nadia's consensus
accuracySlide$NA_PtAc125 <- apply(slide125[,col.nam$s125], 2, function(x) sum(x == slide125$consensus20) / 300 * 100)[accuracySlide$PersonID]
accuracySlide$NA_PtAc150 <- apply(slide150[,col.nam$s150], 2, function(x) sum(x == slide150$consensus20) / 300 * 100)[accuracySlide$PersonID]
accuracyDigital$NA_PtAc125 <- apply(digital125[,col.nam$d125], 2, function(x) sum(x == digital125$consensus20) / 300 * 100)
accuracyDigital$NA_PtAc150 <- apply(digital150[,col.nam$d150], 2, function(x) sum(x == digital150$consensus20) / 300 * 100)
# n.b. these values agree with Nadia's (with the odd rounding error), which is good

# based on my consensus
accuracySlide$IF_PtAc125 <- apply(slide125[,col.nam$s125], 2, function(x) sum(x == slide125$IFcMin) / 300 * 100)[accuracySlide$PersonID]
accuracySlide$IF_PtAc150 <- apply(slide150[,col.nam$s150], 2, function(x) sum(x == slide150$IFcMin) / 300 * 100)[accuracySlide$PersonID]
accuracyDigital$IF_PtAc125 <- apply(digital125[,col.nam$d125], 2, function(x) sum(x == digital125$IFcMin) / 300 * 100)
accuracyDigital$IF_PtAc150 <- apply(digital150[,col.nam$d150], 2, function(x) sum(x == digital150$IFcMin) / 300 * 100)

# 3b. C20_score_group -----------------------------------------------------
# the mean value for the group
c20_mn <- list()
# Nadia's value
c20_mn$NA_PtAc125s <- mean(accuracySlide$NA_PtAc125)
c20_mn$NA_PtAc150s <- mean(accuracySlide$NA_PtAc150)
c20_mn$NA_PtAc125d <- mean(accuracyDigital$NA_PtAc125)
c20_mn$NA_PtAc150d <- mean(accuracyDigital$NA_PtAc150)
# My value
c20_mn$IF_PtAc125s <- mean(accuracySlide$IF_PtAc125)
c20_mn$IF_PtAc150s <- mean(accuracySlide$IF_PtAc150)
c20_mn$IF_PtAc125d <- mean(accuracyDigital$IF_PtAc125)
c20_mn$IF_PtAc150d <- mean(accuracyDigital$IF_PtAc150)

# and sd
c20_sd <- list()
c20_sd$NA_PtAc125s <- sd(accuracySlide$NA_PtAc125)
c20_sd$NA_PtAc150s <- sd(accuracySlide$NA_PtAc150)
c20_sd$NA_PtAc125d <- sd(accuracyDigital$NA_PtAc125)
c20_sd$NA_PtAc150d <- sd(accuracyDigital$NA_PtAc150)
# My value
c20_sd$IF_PtAc125s <- sd(accuracySlide$IF_PtAc125)
c20_sd$IF_PtAc150s <- sd(accuracySlide$IF_PtAc150)
c20_sd$IF_PtAc125d <- sd(accuracyDigital$IF_PtAc125)
c20_sd$IF_PtAc150d <- sd(accuracyDigital$IF_PtAc150)

# mean consensus value based on doing a full pairwise comparison - these are the same as calculated above, so no need to do them.
sum(slide125$IFcMin == slide125[,col.nam$s125]) / (300*(length(col.nam$s125))) * 100
c20_mn$IF_PtAc125s
sum(slide150$IFcMin == slide150[,col.nam$s150]) / (300*(length(col.nam$s150))) * 100
c20_mn$IF_PtAc150s
sum(digital125$IFcMin == digital125[,col.nam$d125]) / (300*(length(col.nam$d125))) * 100
c20_mn$IF_PtAc125d
sum(digital150$IFcMin == digital150[,col.nam$d150]) / (300*(length(col.nam$d150))) * 100
c20_mn$IF_PtAc150d

# 3c. Average pairwise agreement scores -----------------------------------
# mean
# this is the mean similarity between all the columns
# the -300 and -1 come because the column shouldn't be compared with itself (it will be 100% similar)
accuracySlide$mnPA125 <- apply(slide125[,col.nam$s125], 2, function(x) (sum(x == slide125[,col.nam$s125]) - 300) / (300*(length(col.nam$s125) - 1)) * 100)[accuracySlide$PersonID]
accuracySlide$mnPA150 <- apply(slide150[,col.nam$s150], 2, function(x) (sum(x == slide150[,col.nam$s150]) - 300) / (300*(length(col.nam$s150) - 1)) * 100)[accuracySlide$PersonID]
accuracyDigital$mnPA125 <- apply(digital125[,col.nam$d125], 2, function(x) (sum(x == digital125[,col.nam$d125]) - 300) / (300*(length(col.nam$d125) - 1)) * 100)
accuracyDigital$mnPA150 <- apply(digital150[,col.nam$d150], 2, function(x) (sum(x == digital150[,col.nam$d150]) - 300) / (300*(length(col.nam$d150) - 1)) * 100)
# quite a lot of differences here, though fewer in the digital ones

# SD
accuracySlide$sdPA125 <- NA
accuracySlide$sdPA150 <- NA

for (i in accuracySlide$PersonID) {
  # the standard deviation on the scale of the mean (i.e. as a percentage)
  accuracySlide$sdPA125[accuracySlide$PersonID == i] <- sd(apply(slide125[,col.nam$s125][, i] == slide125[,col.nam$s125], 2, sum)[which(names(slide125[,col.nam$s125]) != i)])/300*100
  accuracySlide$sdPA150[accuracySlide$PersonID == i] <- sd(apply(slide150[,col.nam$s150][, i] == slide150[,col.nam$s150], 2, sum)[which(names(slide150[,col.nam$s150]) != i)])/300*100
}
rm(i)

accuracyDigital$sdPA125 <- NA
accuracyDigital$sdPA150 <- NA

for (i in accuracyDigital$PersonID) {
  accuracyDigital$sdPA125[accuracyDigital$PersonID == i] <- sd(apply(digital125[,col.nam$d125][, i] == digital125[,col.nam$d125], 2, sum)[which(names(digital125[,col.nam$d125]) != i)])/300*100
  accuracyDigital$sdPA150[accuracyDigital$PersonID == i] <- sd(apply(digital150[,col.nam$d150][, i] == digital150[,col.nam$d150], 2, sum)[which(names(digital150[,col.nam$d150]) != i)])/300*100
}
rm(i)

# How does this depend on whether workers routinely count specimens? 
# add in routine to the accuracy info
accuracySlide$Routine <- people$Routine[match(accuracySlide$PersonID, people$SlideID)]
accuracySlide$Routine[accuracySlide$PersonID %in% c("1a", "2a")] <- people$Routine[1:2]
accuracySlide$Routine[accuracySlide$PersonID %in% c("1b", "2b")] <- people$Routine[1:2]
accuracyDigital$Routine <- people$Routine[match(accuracyDigital$PersonID, people$DigitalID)]

tapply(accuracySlide$IF_PtAc125, accuracySlide$Routine, summary)
tapply(accuracySlide$IF_PtAc150, accuracySlide$Routine, summary)

tapply(accuracyDigital$IF_PtAc125, accuracyDigital$Routine, summary)
tapply(accuracyDigital$IF_PtAc150, accuracyDigital$Routine, summary)

png("Figures/Routine_accuracy.png")
par(mfrow = c(2,2))
boxplot(accuracySlide$IF_PtAc125 ~ accuracySlide$Routine, main = "Slide 125")
boxplot(accuracySlide$IF_PtAc150 ~ accuracySlide$Routine, main = "Slide 150")
boxplot(accuracyDigital$IF_PtAc125 ~ accuracyDigital$Routine, main = "Digital 125")
boxplot(accuracyDigital$IF_PtAc150 ~ accuracyDigital$Routine, main = "Digital 150")
par(mfrow = c(1,1))
dev.off()

# 3d. Look at this accuracy as plots ----------------------------
# compare these as plots
par(mfrow = c(2, 2))
# how does Nadia's results compare with mine
plot(factor(accuracySlide$PersonID), accuracySlide$NA_PtAc125, main = "Percentage Accuracy 125", ylim = c(45, 90))
points(factor(accuracySlide$PersonID), accuracySlide$IF_PtAc125, col = "blue", pch = 16)
abline(h = c20_mn$NA_PtAc125s)
abline(h = c(c20_mn$NA_PtAc125s - c20_sd$NA_PtAc125s, c20_mn$NA_PtAc125s + c20_sd$NA_PtAc125s), lty = 2)
abline(h = c20_mn$IF_PtAc125s, col = "blue")
abline(h = c(c20_mn$IF_PtAc125s - c20_sd$IF_PtAc125s, c20_mn$IF_PtAc125s + c20_sd$IF_PtAc125s), lty = 2, col = "blue")

plot(factor(accuracySlide$PersonID), accuracySlide$NA_PtAc150, main = "Percentage Accuracy 150", ylim = c(45, 90))
points(factor(accuracySlide$PersonID), accuracySlide$IF_PtAc150, col = "blue", pch = 16)
abline(h = c20_mn$NA_PtAc150s)
abline(h = c(c20_mn$NA_PtAc150s - c20_sd$NA_PtAc150s, c20_mn$NA_PtAc150s + c20_sd$NA_PtAc150s), lty = 2)
abline(h = c20_mn$IF_PtAc150s, col = "blue")
abline(h = c(c20_mn$IF_PtAc150s - c20_sd$IF_PtAc150s, c20_mn$IF_PtAc150s + c20_sd$IF_PtAc150s), lty = 2, col = "blue")

plot(factor(accuracyDigital$PersonID), accuracyDigital$NA_PtAc125, main = "Percentage Accuracy 125", ylim = c(45, 90))
points(factor(accuracyDigital$PersonID), accuracyDigital$IF_PtAc125, col = "blue", pch = 16)
abline(h = c20_mn$NA_PtAc125d)
abline(h = c(c20_mn$NA_PtAc125d - c20_sd$A_PA125d, c20_mn$NA_PtAc125d + c20_sd$NA_PtAc125d), lty = 2)
abline(h = c20_mn$IF_PtAc125d, col = "blue")
abline(h = c(c20_mn$IF_PtAc125d - c20_sd$IF_PtAc125d, c20_mn$IF_PtAc125d + c20_sd$IF_PtAc125d), lty = 2, col = "blue")

plot(factor(accuracyDigital$PersonID), accuracyDigital$NA_PtAc150, main = "Percentage Accuracy 150", ylim = c(45, 90))
points(factor(accuracyDigital$PersonID), accuracyDigital$IF_PtAc150, col = "blue", pch = 16)
abline(h = c20_mn$NA_PtAc150d)
abline(h = c(c20_mn$NA_PtAc150d - c20_sd$NA_PtAc150d, c20_mn$NA_PtAc150d + c20_sd$NA_PtAc150d), lty = 2)
abline(h = c20_mn$IF_PtAc150d, col = "blue")
abline(h = c(c20_mn$IF_PtAc150d - c20_sd$IF_PtAc150d, c20_mn$IF_PtAc150d + c20_sd$IF_PtAc150d), lty = 2, col = "blue")
par(mfrow = c(1,1))

# recreate figure 3
err_bar <- function(mean, sd, xpos, length = 0.05, col = 1) {
  for(i in 1:length(mean)) {
    if (length(col) > 1) {
      arrows(xpos[i], mean[i] - sd[i], xpos[i], mean[i] + sd[i], angle = 90, code = 3, length = length, col = col[i])
    } else {
      arrows(xpos[i], mean[i] - sd[i], xpos[i], mean[i] + sd[i], angle = 90, code = 3, length = length, col = col)
    }
    
  }
}

accuracySlide$Experience <- people$ExperienceSlideA[match(accuracySlide$PersonID, people$SlideID)]
accuracySlide$Experience[accuracySlide$PersonID %in% c("1a", "2a")] <- people$ExperienceSlideA[1:2]
accuracySlide$Experience[accuracySlide$PersonID %in% c("1b", "2b")] <- people$ExperienceSlideB[1:2]
accuracyDigital$Experience <- people$ExperienceDigital[match(accuracyDigital$PersonID, people$DigitalID)]

png("Figures/Fig3_Pairwise_agreement.png", 500, 800)
par(mfrow = c(4, 1), mar = c(2.5, 4.1, .5, 1))
with(accuracySlide, plot(Experience, mnPA125, ylim = c(30, 90), type = "n", xlim = c(0, 40)))
rect(-5, c20_mn$IF_PtAc125s - c20_sd$IF_PtAc125s, 45, c20_mn$IF_PtAc125s + c20_sd$IF_PtAc125s, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_mn$IF_PtAc125s, lty = 4)
with(accuracySlide, points(Experience, mnPA125, pch = 16))
with(accuracySlide, err_bar(mnPA125, sdPA125, Experience))
with(accuracySlide[1, ], text(Experience - 0.5, mnPA125, labels = PersonID, cex = 0.7))
with(accuracySlide[2:nrow(accuracySlide), ], text(Experience + 0.5, mnPA125, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(mnPA125[PersonID == "1a"], mnPA125[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(mnPA125[PersonID == "2a"], mnPA125[PersonID == "2b"])))
text(38, 35, "Slide 125", cex = 1.5)

with(accuracySlide, plot(Experience, mnPA150, ylim = c(30, 90), type = "n", xlim = c(0, 40)))
rect(-5, c20_mn$IF_PtAc150s - c20_sd$IF_PtAc150s, 45, c20_mn$IF_PtAc150s + c20_sd$IF_PtAc150s, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_mn$IF_PtAc150s, lty = 4)
with(accuracySlide, points(Experience, mnPA150, pch = 16))
with(accuracySlide, err_bar(mnPA150, sdPA150, Experience))
with(accuracySlide[1, ], text(Experience - 0.5, mnPA150, labels = PersonID, cex = 0.7))
with(accuracySlide[2:nrow(accuracySlide), ], text(Experience + 0.5, mnPA150, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(mnPA150[PersonID == "1a"], mnPA150[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(mnPA150[PersonID == "2a"], mnPA150[PersonID == "2b"])))
text(38, 35, "Slide 150", cex = 1.5)

with(accuracyDigital, plot(Experience, mnPA125, ylim = c(30, 90), type = "n", xlim = c(0, 40)))
rect(-5, c20_mn$IF_PtAc125d - c20_sd$IF_PtAc125d, 45, c20_mn$IF_PtAc125d + c20_sd$IF_PtAc125d, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_mn$IF_PtAc125d, lty = 4)
with(accuracyDigital, points(Experience, mnPA125, pch = 16))
with(accuracyDigital, err_bar(mnPA125, sdPA125, Experience))
with(accuracyDigital[1, ], text(Experience - 0.5, mnPA125, labels = PersonID, cex = 0.7))
with(accuracyDigital[2:nrow(accuracyDigital), ], text(Experience + 0.5, mnPA125, labels = PersonID, cex = 0.7))
text(38, 35, "Digital 125", cex = 1.5)

with(accuracyDigital, plot(Experience, mnPA150, ylim = c(30, 90), type = "n", xlim = c(0, 40)))
rect(-5, c20_mn$IF_PtAc150d - c20_sd$IF_PtAc150d, 45, c20_mn$IF_PtAc150d + c20_sd$IF_PtAc150d, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_mn$IF_PtAc150d, lty = 4)
with(accuracyDigital, points(Experience, mnPA150, pch = 16))
with(accuracyDigital, err_bar(mnPA150, sdPA150, Experience))
with(accuracyDigital[1, ], text(Experience - 0.5, mnPA150, labels = PersonID, cex = 0.7))
with(accuracyDigital[2:nrow(accuracyDigital), ], text(Experience + 0.5, mnPA150, labels = PersonID, cex = 0.7))
text(38, 35, "Digital 150", cex = 1.5)
par(mfrow = c(1,1))
dev.off()

# is there a more useful way to plot this? What does it actually mean?
# compare consensus PA and MPA
png("Figures/cPtAc_mnPA.png")
plot(1, xlim = c(40, 90), ylim = c(38, 74), type = "n", xlab = "Consensus Percentage Accuracy", ylab = "Mean pairwise agreement")
points(accuracySlide$IF_PtAc125, accuracySlide$mnPA125, pch = 4)
abline(lm(accuracySlide$mnPA125 ~ accuracySlide$IF_PtAc125))

points(accuracySlide$NA_PtAc125, accuracySlide$mnPA125, pch = 3)
abline(lm(accuracySlide$mnPA125 ~ accuracySlide$NA_PtAc125), lty = 2)

points(accuracySlide$IF_PtAc150, accuracySlide$mnPA150, pch = 4, col = "red")
abline(lm(accuracySlide$mnPA150 ~ accuracySlide$IF_PtAc150), col = "red")

points(accuracySlide$NA_PtAc150, accuracySlide$mnPA150, pch = 3, col = "red")
abline(lm(accuracySlide$mnPA150 ~ accuracySlide$NA_PtAc150), lty = 2, col = "red")

points(accuracyDigital$IF_PtAc125, accuracyDigital$mnPA125, pch = 4, col = "blue")
abline(lm(accuracyDigital$mnPA125 ~ accuracyDigital$IF_PtAc125), col = "blue")

points(accuracyDigital$NA_PtAc125, accuracyDigital$mnPA125, pch = 3, col = "blue")
abline(lm(accuracyDigital$mnPA125 ~ accuracyDigital$NA_PtAc125), lty = 2, col = "blue")

points(accuracyDigital$IF_PtAc150, accuracyDigital$mnPA150, pch = 4, col = "purple")
abline(lm(accuracyDigital$mnPA150 ~ accuracyDigital$IF_PtAc150), col = "purple")

points(accuracyDigital$NA_PtAc150, accuracyDigital$mnPA150, pch = 3, col = "purple")
abline(lm(accuracyDigital$mnPA150 ~ accuracyDigital$NA_PtAc150), lty = 2, col = "purple")

legend("topleft", legend = c("Slide125", "Slide150", "Digital125", "Digital150", "Isabel (x)", "Nadia (+)"), lty = c(rep(1, 4), 1, 2), col = c("black", "red", "blue", "purple", "black", "black"))
dev.off()

# is it more useful to plot each person's agreement relative to the consensus?
png("Figures/Fig3_Consensus_agreement.png", 500, 800)
par(mfrow = c(4, 1), mar = c(2.5, 4.1, .5, 1))
with(accuracySlide, plot(Experience, IF_PtAc125, ylim = c(45, 90), type = "n", xlim = c(0, 40), ylab = "Percentage Accuracy"))
rect(-5, c20_mn$IF_PtAc125s - c20_sd$IF_PtAc125s, 45, c20_mn$IF_PtAc125s + c20_sd$IF_PtAc125s, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_mn$IF_PtAc125s, lty = 4)
with(accuracySlide, points(Experience, IF_PtAc125, pch = 16))
with(accuracySlide, text(Experience + 0.5, IF_PtAc125, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(IF_PtAc125[PersonID == "1a"], IF_PtAc125[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(IF_PtAc125[PersonID == "2a"], IF_PtAc125[PersonID == "2b"])))
text(38, 48, "Slide 125", cex = 1.5)

with(accuracySlide, plot(Experience, IF_PtAc150, ylim = c(45, 90), type = "n", xlim = c(0, 40), ylab = "Percentage Accuracy"))
rect(-5, c20_mn$IF_PtAc150s - c20_sd$IF_PtAc150s, 45, c20_mn$IF_PtAc150s + c20_sd$IF_PtAc150s, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_mn$IF_PtAc150s, lty = 4)
with(accuracySlide, points(Experience, IF_PtAc150, pch = 16))
with(accuracySlide, text(Experience + 0.5, IF_PtAc150, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(IF_PtAc150[PersonID == "1a"], IF_PtAc150[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(IF_PtAc150[PersonID == "2a"], IF_PtAc150[PersonID == "2b"])))
text(38, 48, "Slide 150", cex = 1.5)

with(accuracyDigital, plot(Experience, IF_PtAc125, ylim = c(45, 90), type = "n", xlim = c(0, 40), ylab = "Percentage Accuracy"))
rect(-5, c20_mn$IF_PtAc125d - c20_sd$IF_PtAc125d, 45, c20_mn$IF_PtAc125d + c20_sd$IF_PtAc125d, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_mn$IF_PtAc125d, lty = 4)
with(accuracyDigital, points(Experience, IF_PtAc125, pch = 16))
with(accuracyDigital, text(Experience + 0.5, IF_PtAc125, labels = PersonID, cex = 0.7))
text(38, 48, "Digital 125", cex = 1.5)

with(accuracyDigital, plot(Experience, IF_PtAc150, ylim = c(45, 90), type = "n", xlim = c(0, 40), ylab = "Percentage Accuracy"))
rect(-5, c20_mn$IF_PtAc150d - c20_sd$IF_PtAc150d, 45, c20_mn$IF_PtAc150d + c20_sd$IF_PtAc150d, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_mn$IF_PtAc150d, lty = 4)
with(accuracyDigital, points(Experience, IF_PtAc150, pch = 16))
with(accuracyDigital, text(Experience + 0.5, IF_PtAc150, labels = PersonID, cex = 0.7))
text(38, 48, "Digital 150", cex = 1.5)
par(mfrow = c(1,1))
dev.off()

# or plotted against Person ID, not experience
png("Figures/Fig3_Consensus_agreement_ID.png", 500, 800)
par(mfrow = c(4, 1), mar = c(2.5, 4.1, .5, 1))
with(accuracySlide, plot(1:17, IF_PtAc125, ylim = c(45, 90), type = "n", xaxt = "n", xlab = "", ylab = "Percentage Accuracy"))
axis(1, at = 1:17, labels = accuracySlide$PersonID)
rect(0, c20_mn$IF_PtAc125s - c20_sd$IF_PtAc125s, 18, c20_mn$IF_PtAc125s + c20_sd$IF_PtAc125s, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_mn$IF_PtAc125s, lty = 4)
with(accuracySlide, points(1:17, IF_PtAc125, pch = 16))
with(accuracySlide, lines(c(1:2), c(IF_PtAc125[PersonID == "1a"], IF_PtAc125[PersonID == "1b"])))
with(accuracySlide, lines(c(3:4), c(IF_PtAc125[PersonID == "2a"], IF_PtAc125[PersonID == "2b"])))
text(2, 48, "Slide 125", cex = 1.5)

with(accuracySlide, plot(1:17, IF_PtAc150, ylim = c(45, 90), type = "n", xaxt = "n", xlab = "", ylab = "Percentage Accuracy"))
axis(1, at = 1:17, labels = accuracySlide$PersonID)
rect(-5, c20_mn$IF_PtAc150s - c20_sd$IF_PtAc150s, 45, c20_mn$IF_PtAc150s + c20_sd$IF_PtAc150s, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_mn$IF_PtAc150s, lty = 4)
with(accuracySlide, points(1:17, IF_PtAc150, pch = 16))
with(accuracySlide, lines(c(1:2), c(IF_PtAc150[PersonID == "1a"], IF_PtAc150[PersonID == "1b"])))
with(accuracySlide, lines(c(3:4), c(IF_PtAc150[PersonID == "2a"], IF_PtAc150[PersonID == "2b"])))
text(2, 48, "Slide 150", cex = 1.5)

with(accuracyDigital, plot(1:9, IF_PtAc125, ylim = c(45, 90), type = "n", xaxt = "n", xlab = "", ylab = "Percentage Accuracy"))
axis(1, at = 1:9, labels = accuracyDigital$PersonID)
rect(-5, c20_mn$IF_PtAc125d - c20_sd$IF_PtAc125d, 45, c20_mn$IF_PtAc125d + c20_sd$IF_PtAc125d, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_mn$IF_PtAc125d, lty = 4)
with(accuracyDigital, points(1:9, IF_PtAc125, pch = 16))
text(1.5, 48, "Digital 125", cex = 1.5)

with(accuracyDigital, plot(1:9, IF_PtAc150, ylim = c(45, 90), type = "n", xaxt = "n", xlab = "", ylab = "Percentage Accuracy"))
axis(1, at = 1:9, labels = accuracyDigital$PersonID)
rect(-5, c20_mn$IF_PtAc150d - c20_sd$IF_PtAc150d, 45, c20_mn$IF_PtAc150d + c20_sd$IF_PtAc150d, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_mn$IF_PtAc150d, lty = 4)
with(accuracyDigital, points(1:9, IF_PtAc150, pch = 16))
text(1.5, 48, "Digital 150", cex = 1.5)
par(mfrow = c(1,1))
dev.off()

# Digital and slide on one plot
ord.div <- c("1a", "1b", "2a", "2b", "A", 3:6, "F", 7:9, "G", 10:15, LETTERS[2:5], LETTERS[8:9])
accuracyFull <- rbind(accuracySlide, accuracyDigital)
accuracyFull$Analysis <- c(rep("Slide", 17), rep("Digital", 9))

png("Figures/Fig3_Consensus_agreement_fullID.png", 800, 1000)
par(mfrow = c(2, 1))
with(accuracyFull, plot(IF_PtAc125[match(ord.div, PersonID)], ylim = c(43, 90), type = "n", xaxt = "n", ylab = "Percentage accuracy", xlab = "Participant", cex.lab = 1.5, las = 2, cex.axis = 1.1))
axis(1, at = 1:26, labels = accuracyFull$PersonID[match(ord.div, accuracyFull$PersonID)], cex.axis = 1.1)
abline(h = c20_mn$IF_PtAc125s, lty = 1)
abline(h = c(c20_mn$IF_PtAc125s - c20_sd$IF_PtAc125s, c20_mn$IF_PtAc125s + c20_sd$IF_PtAc125s), lty = 4)
abline(h = c20_mn$IF_PtAc125d, lty = 1, col = "blue")
abline(h = c(c20_mn$IF_PtAc125d - c20_sd$IF_PtAc125d, c20_mn$IF_PtAc125d + c20_sd$IF_PtAc125d), lty = 4, col = "blue")
with(accuracyFull, points(1:26, IF_PtAc125[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
with(accuracyFull, arrows(1, IF_PtAc125[PersonID == "1a"], 2, IF_PtAc125[PersonID == "1b"], length = 0.14))
with(accuracyFull, arrows(3, IF_PtAc125[PersonID == "2a"], 4, IF_PtAc125[PersonID == "2b"], length = 0.14))
text(26, 90, "125", cex = 1.5)

with(accuracyFull, plot(IF_PtAc150[match(ord.div, PersonID)], ylim = c(43, 90), type = "n", xaxt = "n", ylab = "Percentage accuracy", xlab = "Participant", cex.lab = 1.5, las = 2, cex.axis = 1.1))
axis(1, at = 1:26, labels = accuracyFull$PersonID[match(ord.div, accuracyFull$PersonID)], cex.axis = 1.1)
abline(h = c20_mn$IF_PtAc150s, lty = 1)
abline(h = c(c20_mn$IF_PtAc150s - c20_sd$IF_PtAc150s, c20_mn$IF_PtAc150s + c20_sd$IF_PtAc150s), lty = 4)
abline(h = c20_mn$IF_PtAc150d, lty = 1, col = "blue")
abline(h = c(c20_mn$IF_PtAc150d - c20_sd$IF_PtAc150d, c20_mn$IF_PtAc150d + c20_sd$IF_PtAc150d), lty = 4, col = "blue")
with(accuracyFull, points(1:26, IF_PtAc150[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
with(accuracyFull, arrows(1, IF_PtAc150[PersonID == "1a"], 2, IF_PtAc150[PersonID == "1b"], length = 0.14))
with(accuracyFull, arrows(3, IF_PtAc150[PersonID == "2a"], 4, IF_PtAc150[PersonID == "2b"], length = 0.14))
text(26, 90, "150", cex = 1.5)

par(mfrow = c(1,1))
dev.off()

# looking at the summary of this data
tapply(accuracyFull$IF_PtAc125, accuracyFull$Analysis, summary)
tapply(accuracyFull$IF_PtAc150, accuracyFull$Analysis, summary)


# 3e. Lumped pairwise agreement scores ------------------------------------
# Table 4
# what happens to agreements if morphologically similar species are grouped. So:
lump <- list()
lump$s125 <- slide125
lump$s150 <- slide150
lump$d125 <- digital125
lump$d150 <- digital150

lump$s125[lump$s125 == "falc"] <- "bull" # bulloides with falconensis
lump$s125[lump$s125 == "cal"] <- "siph"# siphonifera with calida
lump$s125[lump$s125 == "tri"] <- "sac"# sacculifer with trilobus

lump$s150[lump$s150 == "falc"] <- "bull" # bulloides with falconensis
lump$s150[lump$s150 == "cal"] <- "siph"# siphonifera with calida
lump$s150[lump$s150 == "tri"] <- "sac"# sacculifer with trilobus

lump$d125[lump$d125 == "falc"] <- "bull" # bulloides with falconensis
lump$d125[lump$d125 == "cal"] <- "siph"# siphonifera with calida
lump$d125[lump$d125 == "tri"] <- "sac"# sacculifer with trilobus

lump$d150[lump$d150 == "falc"] <- "bull" # bulloides with falconensis
lump$d150[lump$d150 == "cal"] <- "siph"# siphonifera with calida
lump$d150[lump$d150 == "tri"] <- "sac"# sacculifer with trilobus

table(lump$s125$IFcMin)
table(lump$s150$IFcMin)

# now rerun the previous analyses
# based on Nadia's consensus
accuracySlide$l.NA_PtAc125 <- apply(lump$s125[,col.nam$s125], 2, function(x) sum(x == lump$s125$consensus20) / 300 * 100)[accuracySlide$PersonID]
accuracySlide$l.NA_PtAc150 <- apply(lump$s150[,col.nam$s150], 2, function(x) sum(x == lump$s150$consensus20) / 300 * 100)[accuracySlide$PersonID]
accuracyDigital$l.NA_PtAc125 <- apply(lump$d125[,col.nam$d125], 2, function(x) sum(x == lump$d125$consensus20) / 300 * 100)
accuracyDigital$l.NA_PtAc150 <- apply(lump$d150[,col.nam$d150], 2, function(x) sum(x == lump$d150$consensus20) / 300 * 100)
# n.b. these values agree with Nadia's (with the odd rounding error), which is good

# based on my consensus
accuracySlide$l.IF_PtAc125 <- apply(lump$s125[,col.nam$s125], 2, function(x) sum(x == lump$s125$IFcMin) / 300 * 100)[accuracySlide$PersonID]
accuracySlide$l.IF_PtAc150 <- apply(lump$s150[,col.nam$s150], 2, function(x) sum(x == lump$s150$IFcMin) / 300 * 100)[accuracySlide$PersonID]
accuracyDigital$l.IF_PtAc125 <- apply(lump$d125[,col.nam$d125], 2, function(x) sum(x == lump$d125$IFcMin) / 300 * 100)
accuracyDigital$l.IF_PtAc150 <- apply(lump$d150[,col.nam$d150], 2, function(x) sum(x == lump$d150$IFcMin) / 300 * 100)

# C20_score_group
# the mean value for the group
c20_lmn <- list()
# Nadia's value
c20_lmn$NA_PtAc125s <- mean(accuracySlide$l.NA_PtAc125)
c20_lmn$NA_PtAc150s <- mean(accuracySlide$l.NA_PtAc150)
c20_lmn$NA_PtAc125d <- mean(accuracyDigital$l.NA_PtAc125)
c20_lmn$NA_PtAc150d <- mean(accuracyDigital$l.NA_PtAc150)
# My value
c20_lmn$IF_PtAc125s <- mean(accuracySlide$l.IF_PtAc125)
c20_lmn$IF_PtAc150s <- mean(accuracySlide$l.IF_PtAc150)
c20_lmn$IF_PtAc125d <- mean(accuracyDigital$l.IF_PtAc125)
c20_lmn$IF_PtAc150d <- mean(accuracyDigital$l.IF_PtAc150)

# and sd
c20_lsd <- list()
c20_lsd$NA_PtAc125s <- sd(accuracySlide$l.NA_PtAc125)
c20_lsd$NA_PtAc150s <- sd(accuracySlide$l.NA_PtAc150)
c20_lsd$NA_PtAc125d <- sd(accuracyDigital$l.NA_PtAc125)
c20_lsd$NA_PtAc150d <- sd(accuracyDigital$l.NA_PtAc150)
# My value
c20_lsd$IF_PtAc125s <- sd(accuracySlide$l.IF_PtAc125)
c20_lsd$IF_PtAc150s <- sd(accuracySlide$l.IF_PtAc150)
c20_lsd$IF_PtAc125d <- sd(accuracyDigital$l.IF_PtAc125)
c20_lsd$IF_PtAc150d <- sd(accuracyDigital$l.IF_PtAc150)

# mean consensus value based on doing a full pairwise comparison - these are the same as calculated above, so no need to do them.
sum(lump$s125$IFcMin == lump$s125[,col.nam$s125]) / (300*(length(col.nam$s125))) * 100
c20_lmn$IF_PtAc125s
sum(lump$s150$IFcMin == lump$s150[,col.nam$s150]) / (300*(length(col.nam$s150))) * 100
c20_lmn$IF_PtAc150s
sum(lump$d125$IFcMin == lump$d125[,col.nam$d125]) / (300*(length(col.nam$d125))) * 100
c20_lmn$IF_PtAc125d
sum(lump$d150$IFcMin == lump$d150[,col.nam$d150]) / (300*(length(col.nam$d150))) * 100
c20_lmn$IF_PtAc150d


# Average pairwise agreement scores
# mean
accuracySlide$l.mnPA125 <- apply(lump$s125[,col.nam$s125], 2, function(x) (sum(x == lump$s125[,col.nam$s125]) - 300) / (300*(length(col.nam$s125) - 1)) * 100)[accuracySlide$PersonID]
accuracySlide$l.mnPA150 <- apply(lump$s150[,col.nam$s150], 2, function(x) (sum(x == lump$s150[,col.nam$s150]) - 300) / (300*(length(col.nam$s150) - 1)) * 100)[accuracySlide$PersonID]
accuracyDigital$l.mnPA125 <- apply(lump$d125[,col.nam$d125], 2, function(x) (sum(x == lump$d125[,col.nam$d125]) - 300) / (300*(length(col.nam$d125) - 1)) * 100)
accuracyDigital$l.mnPA150 <- apply(lump$d150[,col.nam$d150], 2, function(x) (sum(x == lump$d150[,col.nam$d150]) - 300) / (300*(length(col.nam$d150) - 1)) * 100)
# quite a lot of differences here, though fewer in the digital ones

# SD
accuracySlide$l.sdPA125 <- NA
accuracySlide$l.sdPA150 <- NA

for (i in accuracySlide$PersonID) {
  accuracySlide$l.sdPA125[accuracySlide$PersonID == i] <- sd(apply(lump$s125[,col.nam$s125][, i] == lump$s125[,col.nam$s125], 2, sum)[which(names(lump$s125[,col.nam$s125]) != i)])/300*100
  accuracySlide$l.sdPA150[accuracySlide$PersonID == i] <- sd(apply(lump$s150[,col.nam$s150][, i] == lump$s150[,col.nam$s150], 2, sum)[which(names(lump$s150[,col.nam$s150]) != i)])/300*100
}
rm(i)

accuracyDigital$l.sdPA125 <- NA
accuracyDigital$l.sdPA150 <- NA

for (i in accuracyDigital$PersonID) {
  accuracyDigital$l.sdPA125[accuracyDigital$PersonID == i] <- sd(apply(lump$d125[,col.nam$d125][, i] == lump$d125[,col.nam$d125], 2, sum)[which(names(lump$d125[,col.nam$d125]) != i)])/300*100
  accuracyDigital$l.sdPA150[accuracyDigital$PersonID == i] <- sd(apply(lump$d150[,col.nam$d150][, i] == lump$d150[,col.nam$d150], 2, sum)[which(names(lump$d150[,col.nam$d150]) != i)])/300*100
}
rm(i)

# Look at this percentage accuracy as plots
# col.MPAre these as plots
par(mfrow = c(2, 2))
plot(factor(accuracySlide$PersonID), accuracySlide$l.NA_PtAc125, main = "Percentage Accuracy 125", ylim = c(45, 90))
points(factor(accuracySlide$PersonID), accuracySlide$l.IF_PtAc125, col = "blue", pch = 16)
abline(h = c20_lmn$NA_PtAc125s)
abline(h = c(c20_lmn$NA_PtAc125s - c20_lsd$NA_PtAc125s, c20_lmn$NA_PtAc125s + c20_lsd$NA_PtAc125s), lty = 2)
abline(h = c20_lmn$IF_PtAc125s, col = "blue")
abline(h = c(c20_lmn$IF_PtAc125s - c20_lsd$IF_PtAc125s, c20_lmn$IF_PtAc125s + c20_lsd$IF_PtAc125s), lty = 2, col = "blue")

plot(factor(accuracySlide$PersonID), accuracySlide$l.NA_PtAc150, main = "Percentage Accuracy 150", ylim = c(45, 90))
points(factor(accuracySlide$PersonID), accuracySlide$l.IF_PtAc150, col = "blue", pch = 16)
abline(h = c20_lmn$NA_PtAc150s)
abline(h = c(c20_lmn$NA_PtAc150s - c20_lsd$NA_PtAc150s, c20_lmn$NA_PtAc150s + c20_lsd$NA_PtAc150s), lty = 2)
abline(h = c20_lmn$IF_PtAc150s, col = "blue")
abline(h = c(c20_lmn$IF_PtAc150s - c20_lsd$IF_PtAc150s, c20_lmn$IF_PtAc150s + c20_lsd$IF_PtAc150s), lty = 2, col = "blue")

plot(factor(accuracyDigital$PersonID), accuracyDigital$l.NA_PtAc125, main = "Percentage Accuracy 125", ylim = c(45, 90))
points(factor(accuracyDigital$PersonID), accuracyDigital$l.IF_PtAc125, col = "blue", pch = 16)
abline(h = c20_lmn$NA_PtAc125d)
abline(h = c(c20_lmn$NA_PtAc125d - c20_lsd$NA_PtAc125d, c20_lmn$NA_PtAc125d + c20_lsd$NA_PtAc125d), lty = 2)
abline(h = c20_lmn$IF_PtAc125d, col = "blue")
abline(h = c(c20_lmn$IF_PtAc125d - c20_lsd$IF_PtAc125d, c20_lmn$IF_PtAc125d + c20_lsd$IF_PtAc125d), lty = 2, col = "blue")

plot(factor(accuracyDigital$PersonID), accuracyDigital$l.NA_PtAc150, main = "Percentage Accuracy 150", ylim = c(45, 90))
points(factor(accuracyDigital$PersonID), accuracyDigital$l.IF_PtAc150, col = "blue", pch = 16)
abline(h = c20_lmn$NA_PtAc150d)
abline(h = c(c20_lmn$NA_PtAc150d - c20_lsd$NA_PtAc150d, c20_lmn$NA_PtAc150d + c20_lsd$NA_PtAc150d), lty = 2)
abline(h = c20_lmn$IF_PtAc150d, col = "blue")
abline(h = c(c20_lmn$IF_PtAc150d - c20_lsd$IF_PtAc150d, c20_lmn$IF_PtAc150d + c20_lsd$IF_PtAc150d), lty = 2, col = "blue")
par(mfrow = c(1,1))

# recreate figure 3
png("Figures/Fig3_Pairwise_agreement_lump.png", 500, 800)
par(mfrow = c(4, 1), mar = c(2.5, 4.1, .5, 1))
with(accuracySlide, plot(Experience, l.mnPA125, ylim = c(30, 90), type = "n", xlim = c(0, 40), ylab = "Lumped Pairwise Agreement"))
rect(-5, c20_lmn$IF_PtAc125s - c20_lsd$IF_PtAc125s, 45, c20_lmn$IF_PtAc125s + c20_lsd$IF_PtAc125s, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_lmn$IF_PtAc125s, lty = 4)
with(accuracySlide, points(Experience, l.mnPA125, pch = 16))
with(accuracySlide, err_bar(l.mnPA125, l.sdPA125, Experience))
with(accuracySlide, text(Experience + 0.5, l.mnPA125, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(l.mnPA125[PersonID == "1a"], l.mnPA125[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(l.mnPA125[PersonID == "2a"], l.mnPA125[PersonID == "2b"])))
text(38, 35, "Slide 125", cex = 1.5)

with(accuracySlide, plot(Experience, l.mnPA150, ylim = c(30, 90), type = "n", xlim = c(0, 40), ylab = "Lumped Pairwise Agreement"))
rect(-5, c20_lmn$IF_PtAc150s - c20_lsd$IF_PtAc150s, 45, c20_lmn$IF_PtAc150s + c20_lsd$IF_PtAc150s, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_lmn$IF_PtAc150s, lty = 4)
with(accuracySlide, points(Experience, l.mnPA150, pch = 16))
with(accuracySlide, err_bar(l.mnPA150, l.sdPA150, Experience))
with(accuracySlide, text(Experience + 0.5, l.mnPA150, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(l.mnPA150[PersonID == "1a"], l.mnPA150[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(l.mnPA150[PersonID == "2a"], l.mnPA150[PersonID == "2b"])))
text(38, 35, "Slide 150", cex = 1.5)

with(accuracyDigital, plot(Experience, l.mnPA125, ylim = c(30, 90), type = "n", xlim = c(0, 40), ylab = "Lumped Pairwise Agreement"))
rect(-5, c20_lmn$IF_PtAc125d - c20_lsd$IF_PtAc125d, 45, c20_lmn$IF_PtAc125d + c20_lsd$IF_PtAc125d, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_lmn$IF_PtAc125d, lty = 4)
with(accuracyDigital, points(Experience, l.mnPA125, pch = 16))
with(accuracyDigital, err_bar(l.mnPA125, l.sdPA125, Experience))
with(accuracyDigital, text(Experience + 0.5, l.mnPA125, labels = PersonID, cex = 0.7))
text(38, 35, "Digital 125", cex = 1.5)

with(accuracyDigital, plot(Experience, l.mnPA150, ylim = c(30, 90), type = "n", xlim = c(0, 40), ylab = "Lumped Pairwise Agreement"))
rect(-5, c20_lmn$IF_PtAc150d - c20_lsd$IF_PtAc150d, 45, c20_lmn$IF_PtAc150d + c20_lsd$IF_PtAc150d, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_lmn$IF_PtAc150d, lty = 4)
with(accuracyDigital, points(Experience, l.mnPA150, pch = 16))
with(accuracyDigital, err_bar(l.mnPA150, l.sdPA150, Experience))
with(accuracyDigital, text(Experience + 0.5, l.mnPA150, labels = PersonID, cex = 0.7))
text(38, 35, "Digital 150", cex = 1.5)
par(mfrow = c(1,1))
dev.off()

# is there a more useful way to plot this? What does it actually mean?
# col.MPAre consensus PA and l.MPA
png("Figures/cPA_mnPA_lump.png")
plot(1, xlim = c(40, 90), ylim = c(38, 74), type = "n", xlab = "Percentage Accuracy", ylab = "Lumped Pairwise Agreement")
points(accuracySlide$l.IF_PtAc125, accuracySlide$l.mnPA125, pch = 4)
abline(lm(accuracySlide$l.mnPA125 ~ accuracySlide$l.IF_PtAc125))

points(accuracySlide$l.NA_PtAc125, accuracySlide$l.mnPA125, pch = 3)
abline(lm(accuracySlide$l.mnPA125 ~ accuracySlide$l.NA_PtAc125), lty = 2)

points(accuracySlide$l.IF_PtAc150, accuracySlide$l.mnPA150, pch = 4, col = "red")
abline(lm(accuracySlide$l.mnPA150 ~ accuracySlide$l.IF_PtAc150), col = "red")

points(accuracySlide$l.NA_PtAc150, accuracySlide$l.mnPA150, pch = 3, col = "red")
abline(lm(accuracySlide$l.mnPA150 ~ accuracySlide$l.NA_PtAc150), lty = 2, col = "red")

points(accuracyDigital$l.IF_PtAc125, accuracyDigital$l.mnPA125, pch = 4, col = "blue")
abline(lm(accuracyDigital$l.mnPA125 ~ accuracyDigital$l.IF_PtAc125), col = "blue")

points(accuracyDigital$l.NA_PtAc125, accuracyDigital$l.mnPA125, pch = 3, col = "blue")
abline(lm(accuracyDigital$l.mnPA125 ~ accuracyDigital$l.NA_PtAc125), lty = 2, col = "blue")

points(accuracyDigital$l.IF_PtAc150, accuracyDigital$l.mnPA150, pch = 4, col = "purple")
abline(lm(accuracyDigital$l.mnPA150 ~ accuracyDigital$l.IF_PtAc150), col = "purple")

points(accuracyDigital$l.NA_PtAc150, accuracyDigital$l.mnPA150, pch = 3, col = "purple")
abline(lm(accuracyDigital$l.mnPA150 ~ accuracyDigital$l.NA_PtAc150), lty = 2, col = "purple")

legend("topleft", legend = c("slide125", "slide150", "digital125", "digital150", "Isabel (x)", "Nadia (+)"), lty = c(rep(1, 4), 1, 2), col = c("black", "red", "blue", "purple", "black", "black"))
dev.off()

# is it more useful to plot each person relative to the consensus?
png("Figures/Fig3_Consensus_agreement_lump.png", 500, 800)
par(mfrow = c(4, 1), mar = c(2.5, 4.1, .5, 1))
with(accuracySlide, plot(Experience, l.IF_PtAc125, ylim = c(45, 95), type = "n", xlim = c(0, 40), ylab = "Lumped Percentage Accuracy"))
rect(-5, c20_lmn$IF_PtAc125s - c20_lsd$IF_PtAc125s, 45, c20_lmn$IF_PtAc125s + c20_lsd$IF_PtAc125s, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_lmn$IF_PtAc125s, lty = 4)
with(accuracySlide, points(Experience, l.IF_PtAc125, pch = 16))
with(accuracySlide, text(Experience + 0.5, l.IF_PtAc125, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(l.IF_PtAc125[PersonID == "1a"], l.IF_PtAc125[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(l.IF_PtAc125[PersonID == "2a"], l.IF_PtAc125[PersonID == "2b"])))
text(38, 50, "Slide 125", cex = 1.5)

with(accuracySlide, plot(Experience, l.IF_PtAc150, ylim = c(45, 95), type = "n", xlim = c(0, 40), ylab = "Lumped Percentage Accuracy"))
rect(-5, c20_lmn$IF_PtAc150s - c20_lsd$IF_PtAc150s, 45, c20_lmn$IF_PtAc150s + c20_lsd$IF_PtAc150s, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_lmn$IF_PtAc150s, lty = 4)
with(accuracySlide, points(Experience, l.IF_PtAc150, pch = 16))
with(accuracySlide, text(Experience + 0.5, l.IF_PtAc150, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(l.IF_PtAc150[PersonID == "1a"], l.IF_PtAc150[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(l.IF_PtAc150[PersonID == "2a"], l.IF_PtAc150[PersonID == "2b"])))
text(38, 50, "Slide 150", cex = 1.5)

with(accuracyDigital, plot(Experience, l.IF_PtAc125, ylim = c(45, 95), type = "n", xlim = c(0, 40), ylab = "Lumped Percentage Accuracy"))
rect(-5, c20_lmn$IF_PtAc125d - c20_lsd$IF_PtAc125d, 45, c20_lmn$IF_PtAc125d + c20_lsd$IF_PtAc125d, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_lmn$IF_PtAc125d, lty = 4)
with(accuracyDigital, points(Experience, l.IF_PtAc125, pch = 16))
with(accuracyDigital, text(Experience + 0.5, l.IF_PtAc125, labels = PersonID, cex = 0.7))
text(38, 50, "Digital 125", cex = 1.5)

with(accuracyDigital, plot(Experience, l.IF_PtAc150, ylim = c(45, 95), type = "n", xlim = c(0, 40), ylab = "Lumped Percentage Accuracy"))
rect(-5, c20_lmn$IF_PtAc150d - c20_lsd$IF_PtAc150d, 45, c20_lmn$IF_PtAc150d + c20_lsd$IF_PtAc150d, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_lmn$IF_PtAc150d, lty = 4)
with(accuracyDigital, points(Experience, l.IF_PtAc150, pch = 16))
with(accuracyDigital, text(Experience + 0.5, l.IF_PtAc150, labels = PersonID, cex = 0.7))
text(38, 50, "Digital 150", cex = 1.5)
par(mfrow = c(1,1))
dev.off()

# or plotted against Person ID, not experience
png("Figures/Fig3_Consensus_agreement_ID_lump.png", 500, 800)
par(mfrow = c(4, 1), mar = c(2.5, 4.1, .5, 1))
with(accuracySlide, plot(1:17, l.IF_PtAc125, ylim = c(45, 95), type = "n", xaxt = "n", xlab = "", ylab = "Lumped Percentage Accuracy"))
axis(1, at = 1:17, labels = accuracySlide$PersonID)
rect(0, c20_lmn$IF_PtAc125s - c20_lsd$IF_PtAc125s, 18, c20_lmn$IF_PtAc125s + c20_lsd$IF_PtAc125s, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_lmn$IF_PtAc125s, lty = 4)
with(accuracySlide, points(1:17, l.IF_PtAc125, pch = 16))
with(accuracySlide, lines(c(1:2), c(l.IF_PtAc125[PersonID == "1a"], l.IF_PtAc125[PersonID == "1b"])))
with(accuracySlide, lines(c(3:4), c(l.IF_PtAc125[PersonID == "2a"], l.IF_PtAc125[PersonID == "2b"])))
text(2, 48, "Slide 125", cex = 1.5)

with(accuracySlide, plot(1:17, l.IF_PtAc150, ylim = c(45, 95), type = "n", xaxt = "n", xlab = "", ylab = "Lumped Percentage Accuracy"))
axis(1, at = 1:17, labels = accuracySlide$PersonID)
rect(-5, c20_lmn$IF_PtAc150s - c20_lsd$IF_PtAc150s, 45, c20_lmn$IF_PtAc150s + c20_lsd$IF_PtAc150s, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_lmn$IF_PtAc150s, lty = 4)
with(accuracySlide, points(1:17, l.IF_PtAc150, pch = 16))
with(accuracySlide, lines(c(1:2), c(l.IF_PtAc150[PersonID == "1a"], l.IF_PtAc150[PersonID == "1b"])))
with(accuracySlide, lines(c(3:4), c(l.IF_PtAc150[PersonID == "2a"], l.IF_PtAc150[PersonID == "2b"])))
text(2, 48, "Slide 150", cex = 1.5)

with(accuracyDigital, plot(1:9, l.IF_PtAc125, ylim = c(45, 95), type = "n", xaxt = "n", xlab = "", ylab = "Lumped Percentage Accuracy"))
axis(1, at = 1:9, labels = accuracyDigital$PersonID)
rect(-5, c20_lmn$IF_PtAc125d - c20_lsd$IF_PtAc125d, 45, c20_lmn$IF_PtAc125d + c20_lsd$IF_PtAc125d, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_lmn$IF_PtAc125d, lty = 4)
with(accuracyDigital, points(1:9, l.IF_PtAc125, pch = 16))
text(1.5, 48, "Digital 125", cex = 1.5)

with(accuracyDigital, plot(1:9, l.IF_PtAc150, ylim = c(45, 95), type = "n", xaxt = "n", xlab = "", ylab = "Lumped Percentage Accuracy"))
axis(1, at = 1:9, labels = accuracyDigital$PersonID)
rect(-5, c20_lmn$IF_PtAc150d - c20_lsd$IF_PtAc150d, 45, c20_lmn$IF_PtAc150d + c20_lsd$IF_PtAc150d, col = rgb(0, .8, .2, alpha = .5))
abline(h = c20_lmn$IF_PtAc150d, lty = 4)
with(accuracyDigital, points(1:9, l.IF_PtAc150, pch = 16))
text(1.5, 48, "Digital 150", cex = 1.5)
par(mfrow = c(1,1))
dev.off()

# combined plots showing both digital and slide on the same scale
accuracyFull <- rbind(accuracySlide, accuracyDigital)
accuracyFull$Analysis <- c(rep("Slide", 17), rep("Digital", 9))

png("Figures/Fig3_Consensus_agreement_fullID_lump.png", 800, 1000)
par(mfrow = c(2, 1))
with(accuracyFull, plot(l.IF_PtAc125[match(ord.div, PersonID)], ylim = c(45, 92), type = "n", xaxt = "n", ylab = "Percentage accuracy", xlab = "Person"))
axis(1, at = 1:26, labels = accuracyFull$PersonID[match(ord.div, accuracyFull$PersonID)])
abline(h = c20_lmn$IF_PtAc125s, lty = 1)
abline(h = c(c20_lmn$IF_PtAc125s - c20_lsd$IF_PtAc125s, c20_lmn$IF_PtAc125s + c20_lsd$IF_PtAc125s), lty = 4)
abline(h = c20_lmn$IF_PtAc125d, lty = 1, col = "blue")
abline(h = c(c20_lmn$IF_PtAc125d - c20_lsd$IF_PtAc125d, c20_lmn$IF_PtAc125d + c20_lsd$IF_PtAc125d), lty = 4, col = "blue")
with(accuracyFull, points(1:26, l.IF_PtAc125[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
with(accuracyFull, lines(c(1:2), c(l.IF_PtAc125[PersonID == "1a"], l.IF_PtAc125[PersonID == "1b"])))
with(accuracyFull, lines(c(3:4), c(l.IF_PtAc125[PersonID == "2a"], l.IF_PtAc125[PersonID == "2b"])))
text(26, 92, "125", cex = 1.5)

with(accuracyFull, plot(l.IF_PtAc150[match(ord.div, PersonID)], ylim = c(45, 92), type = "n", xaxt = "n", ylab = "Percentage accuracy", xlab = "Person"))
axis(1, at = 1:26, labels = accuracyFull$PersonID[match(ord.div, accuracyFull$PersonID)])
abline(h = c20_lmn$IF_PtAc150s, lty = 1)
abline(h = c(c20_lmn$IF_PtAc150s - c20_lsd$IF_PtAc150s, c20_lmn$IF_PtAc150s + c20_lsd$IF_PtAc150s), lty = 4)
abline(h = c20_lmn$IF_PtAc150d, lty = 1, col = "blue")
abline(h = c(c20_lmn$IF_PtAc150d - c20_lsd$IF_PtAc150d, c20_lmn$IF_PtAc150d + c20_lsd$IF_PtAc150d), lty = 4, col = "blue")
with(accuracyFull, points(1:26, l.IF_PtAc150[match(ord.div, PersonID)], pch = 16, col = (Analysis[match(ord.div, PersonID)] == "Digital")*3 + 1))
with(accuracyFull, lines(c(1:2), c(l.IF_PtAc150[PersonID == "1a"], l.IF_PtAc150[PersonID == "1b"])))
with(accuracyFull, lines(c(3:4), c(l.IF_PtAc150[PersonID == "2a"], l.IF_PtAc150[PersonID == "2b"])))
text(26, 92, "150", cex = 1.5)

par(mfrow = c(1,1))
dev.off()



# 3f. Radial plots --------------------------------------------------------

# for Slide125
for(j in sp.abb$Abbreviation) {
  # for each species
  tmp.sp <- j
  # assuming it occurs in the consensus
  tmp.no <- sum(slide125$IFcMin == tmp.sp)
  if (tmp.no > 0) {
    # sorted list of IDs
    tmp.tab <- rev(sort(table(unlist(head(slide125[which(slide125$IFcMin == tmp.sp),col.nam$s125])))))
    # create a radial barplot
    png(paste("Figures/Radial/slide125_", j, ".png", sep = ""))
    plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    for (i in 1:length(tmp.tab)) {
      # each given name is a bar, width relating to how often it is used
      lines(c(0.5, 0.5 + 0.3*sin((i-1) * pi/(length(tmp.tab)/2))), c(0.5, 0.5 + 0.3*cos((i-1) * pi/(length(tmp.tab)/2))), lwd = tmp.tab[i]/tmp.tab[1]*20, col = "orange")
      text(0.5 + 0.35*sin((i-1) * pi/(length(tmp.tab)/2)), 0.5 + 0.35*cos((i-1) * pi/(length(tmp.tab)/2)), names(tmp.tab[i]))
    }
    dev.off()
  }
}
rm(i, j, tmp.sp, tmp.no, tmp.tab)

# slide150
for(j in sp.abb$Abbreviation) {
  # for each species
  tmp.sp <- j
  # assuming it occurs in the consensus
  tmp.no <- sum(slide150$IFcMin == tmp.sp)
  if (tmp.no > 0) {
    # sorted list of IDs
    tmp.tab <- rev(sort(table(unlist(head(slide150[which(slide150$IFcMin == tmp.sp),col.nam$s150])))))
    # create a radial barplot
    png(paste("Figures/Radial/slide150_", j, ".png", sep = ""))
    plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    for (i in 1:length(tmp.tab)) {
      # each given name is a bar, width relating to how often it is used
      lines(c(0.5, 0.5 + 0.3*sin((i-1) * pi/(length(tmp.tab)/2))), c(0.5, 0.5 + 0.3*cos((i-1) * pi/(length(tmp.tab)/2))), lwd = tmp.tab[i]/tmp.tab[1]*20, col = "orange")
      text(0.5 + 0.35*sin((i-1) * pi/(length(tmp.tab)/2)), 0.5 + 0.35*cos((i-1) * pi/(length(tmp.tab)/2)), names(tmp.tab[i]))
    }
    dev.off()
  }
}
rm(i, j, tmp.sp, tmp.no, tmp.tab)


# digital125
for(j in sp.abb$Abbreviation) {
  # for each species
  tmp.sp <- j
  # assuming it occurs in the consensus
  tmp.no <- sum(digital125$IFcMin == tmp.sp)
  if (tmp.no > 0) {
    # sorted list of IDs
    tmp.tab <- rev(sort(table(unlist(head(digital125[which(digital125$IFcMin == tmp.sp),col.nam$d125])))))
    # create a radial barplot
    png(paste("Figures/Radial/digital125_", j, ".png", sep = ""))
    plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    for (i in 1:length(tmp.tab)) {
      # each given name is a bar, width relating to how often it is used
      lines(c(0.5, 0.5 + 0.3*sin((i-1) * pi/(length(tmp.tab)/2))), c(0.5, 0.5 + 0.3*cos((i-1) * pi/(length(tmp.tab)/2))), lwd = tmp.tab[i]/tmp.tab[1]*20, col = "orange")
      text(0.5 + 0.35*sin((i-1) * pi/(length(tmp.tab)/2)), 0.5 + 0.35*cos((i-1) * pi/(length(tmp.tab)/2)), names(tmp.tab[i]))
    }
    dev.off()
  }
}
rm(i, j, tmp.sp, tmp.no, tmp.tab)


# digital150
for(j in sp.abb$Abbreviation) {
  # for each species
  tmp.sp <- j
  # assuming it occurs in the consensus
  tmp.no <- sum(digital150$IFcMin == tmp.sp)
  if (tmp.no > 0) {
    # sorted list of IDs
    tmp.tab <- rev(sort(table(unlist(head(digital150[which(digital150$IFcMin == tmp.sp),col.nam$d150])))))
    # create a radial barplot
    png(paste("Figures/Radial/digital150_", j, ".png", sep = ""))
    plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    for (i in 1:length(tmp.tab)) {
      # each given name is a bar, width relating to how often it is used
      lines(c(0.5, 0.5 + 0.3*sin((i-1) * pi/(length(tmp.tab)/2))), c(0.5, 0.5 + 0.3*cos((i-1) * pi/(length(tmp.tab)/2))), lwd = tmp.tab[i]/tmp.tab[1]*20, col = "orange")
      text(0.5 + 0.35*sin((i-1) * pi/(length(tmp.tab)/2)), 0.5 + 0.35*cos((i-1) * pi/(length(tmp.tab)/2)), names(tmp.tab[i]))
    }
    dev.off()
  }
}
rm(i, j, tmp.sp, tmp.no, tmp.tab)


# Figure 4

# 3g. Confusion matrix ----------------------------------------------------
# initially with slide 125
# this requires the data to be in long format
long <- list()
long$s125 <- reshape(slide125, varying = list(names(slide125)[col.nam$s125]), direction = "long", times = names(slide125)[col.nam$s125], timevar = "Person")
rownames(long$s125) <- 1:nrow(long$s125)
long$s125 <- long$s125[, (names(long$s125) != "id")]
names(long$s125)[names(long$s125) == "1a"] <- "origID"
head(long$s125)
tail(long$s125)

png("Figures/confusion_slide125.png", 1000, 700)
conf_mat(long$s125, "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Slide 125", 1, cex = 2, adj = -0.8,  line = -5)
dev.off() 

# same for the other datasets
# so slide 150
long$s150 <- reshape(slide150, varying = list(names(slide150)[col.nam$s150]), direction = "long", times = names(slide150)[col.nam$s150], timevar = "Person")
rownames(long$s150) <- 1:nrow(long$s150)
long$s150 <- long$s150[, (names(long$s150) != "id")]
names(long$s150)[names(long$s150) == "1a"] <- "origID"
head(long$s150)
tail(long$s150)

png("Figures/confusion_slide150.png", 1000, 600)
conf_mat(long$s150, "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Slide 150", 1, cex = 2, adj = -.8,  line = -5)
dev.off() 


# digital 125
long$d125 <- reshape(digital125, varying = list(names(digital125)[col.nam$d125]), direction = "long", times = names(digital125)[col.nam$d125], timevar = "Person")
rownames(long$d125) <- 1:nrow(long$d125)
long$d125 <- long$d125[, (names(long$d125) != "id")]
names(long$d125)[names(long$d125) == "A"] <- "origID"
head(long$d125)
tail(long$d125)

png("Figures/confusion_digital125.png", 1000, 700)
conf_mat(long$d125, "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Digital 125", 1, cex = 2, adj = -.8,  line = -5)
dev.off() 


# digital 150
long$d150 <- reshape(digital150, varying = list(names(digital150)[col.nam$d150]), direction = "long", times = names(digital150)[col.nam$d150], timevar = "Person")
rownames(long$d150) <- 1:nrow(long$d150)
long$d150 <- long$d150[, (names(long$d150) != "id")]
names(long$d150)[names(long$d150) == "A"] <- "origID"
head(long$d150)
tail(long$d150)

png("Figures/confusion_digital150.png", 1000, 600)
conf_mat(long$d150, "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = FALSE, sp.list = "full", xlab = "Participant ID", ylab = "Consensus ID")
mtext("Digital 150", 1, cex = 2, adj = -.8,  line = -5)
dev.off() 


# 4. NMDS -----------------------------------------------------------------

# 4a. Recreating Nadia's NMDS ---------------------------------------------
# I think this was originally run on the species counts (i.e. working at the community level)
# I attempted to do this and failed to produce the same results

# Figure 2
# for size fraction 125 (slide and digital are combined, with the consensus values)
# the data needs to be the other way round, so transpose it
sp.125 <- merge(slide125sp, digital125sp, by = "species")

trsp <- list()
trsp$NA125sp <- data.frame(t(sp.125[, !(names(sp.125) %in% c("species", "IFmaxCon.x", "IFmaxCon.y", "IFc50.x", "IFcMin.x", "IFc50.y", "IFcMin.y"))]))
names(trsp$NA125sp) <- sp.125$species
rownames(trsp$NA125sp)[nchar(rownames(trsp$NA125sp)) >= 5] <- c("Sc50", "Sc20", "Dc50", "Dc20")
# re-running this creates a certain amount of movement
nmds <- list()
nmds$NA125sp <- metaMDS(trsp$NA125sp)

# create a dataframe for the colours
mds.col <- data.frame(person = rownames(trsp$NA125sp))
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

# plot the NMDS
plot(nmds$NA125sp, type = "n", display = "sites")
points(nmds$NA125sp, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
text(nmds$NA125sp, labels = rownames(trsp$NA125sp))
legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16)

# for 150
sp.150 <- merge(slide150sp, digital150sp, by = "species")
trsp$NA150sp <- data.frame(t(sp.150[, !(names(sp.150) %in% c("species", "IFmaxCon.x", "IFmaxCon.y", "IFc50.x", "IFcMin.x", "IFc50.y", "IFcMin.y"))]))
names(trsp$NA150sp) <- sp.150$species
rownames(trsp$NA150sp)[nchar(rownames(trsp$NA150sp)) >= 5] <- c("Sc50", "Sc20", "Dc50", "Dc20")
nmds$NA150sp <- metaMDS(trsp$NA150sp)
plot(nmds$NA150sp, type = "n", display = "sites")
points(nmds$NA150sp, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
text(nmds$NA150sp, labels = rownames(trsp$NA150sp))
legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16)
# these are currently very different from Nadia's data

# 4b. Using raw data rather than species counts ---------------------------
# again this is run initially for Nadia's consensus values
# I think instead it is more sensible to use the raw data
# try using the original data instead (n.b. I can't seem to do this in PAST)

# for 125
full.125 <- merge(slide125, digital125, by = "Specimen")
trsp$NA125f <- data.frame(t(full.125[, !(names(full.125) %in% c("Specimen", "IFmaxCon.x", "IFmaxCon.y", "IFc50.x", "IFcMin.x", "IFc50.y", "IFcMin.y"))]))
rownames(trsp$NA125f)[nchar(rownames(trsp$NA125f)) >= 5] <- c("Sc50", "Sc20", "Dc50", "Dc20")
nmds$NA125f <- metaMDS(daisy(trsp$NA125f))

# consider the stress of the NMDS
stress <- list()
stress$NA125f <- rep(NA, 10)
for (i in 1:10) {
  stress$NA125f[i] <- metaMDS(daisy(trsp$NA125f), k = i)$stress
}
plot(stress$NA125f, type = "b")  
rm(i)
stressplot(nmds$NA125f)
# looks like between 2 and 3 dimensions would be reasonable

plot(nmds$NA125f, type = "n", display = "sites")
points(nmds$NA125f, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
text(nmds$NA125f, labels = rownames(trsp$NA125f))
legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16)


# and for 150
full.150 <- merge(slide150, digital150, by = "Specimen")
trsp$NA150f <- data.frame(t(full.150[, !(names(full.150) %in% c("Specimen", "IFmaxCon.x", "IFmaxCon.y", "IFc50.x", "IFcMin.x", "IFc50.y", "IFcMin.y"))]))
rownames(trsp$NA150f)[nchar(rownames(trsp$NA150f)) >= 5] <- c("Sc50", "Sc20", "Dc50", "Dc20")
nmds$NA150f <- metaMDS(daisy(trsp$NA150f))

# consider the stress of the NMDS
stress$NA150f <- rep(NA, 10)
for (i in 1:10) {
  stress$NA150f[i] <- metaMDS(daisy(trsp$NA150f), k = i)$stress
}
plot(stress$NA150f, type = "b")  
rm(i)
stressplot(nmds$NA150f)
# again between 2 and 3 dimensions is probably reasonable

# the full plot
plot(nmds$NA150f, type = "n", display = "sites", cex = 1)
points(nmds$NA150f, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
text(nmds$NA150f, labels = rownames(trsp$NA150f))
legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16)

# Given the influence of the outliers, a more informative relationship between these points can be obtained by running the analysis excluding the outliers. 
nmds$NA150z <- metaMDS(daisy(trsp$NA150f[!(rownames(trsp$NA150f) %in% c("3", "C", "E", "G")),]))

# consider the stress of the NMDS
stress$NA150z <- rep(NA, 10)
for (i in 1:10) {
  stress$NA150z[i] <- metaMDS(daisy(trsp$NA150f[!(rownames(trsp$NA150f) %in% c("3", "C", "E", "G")),]), k = i)$stress
}
plot(stress$NA150z, type = "b")
rm(i)
stressplot(nmds$NA150z)
# again between 2 and 3 dimensions is probably reasonable

# the full plot
plot(nmds$NA150z, type = "n", display = "sites", cex = 1)
points(nmds$NA150z, pch = 21, cex = 4, col = mds.col$pair[!(mds.col$person %in% c("3", "C", "E", "G"))], bg = brewer.pal(5, "Set2")[mds.col$sch.col[!(mds.col$person %in% c("3", "C", "E", "G"))]], lwd = 2)
text(nmds$NA150z, labels = rownames(trsp$NA150f[!(rownames(trsp$NA150f) %in% c("3", "C", "E", "G")),]))
legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16)

# 4c. plotting using my consensus estimates -------------------------------
# for 125
trsp$IF125f <- data.frame(t(full.125[, !(names(full.125) %in% c("Specimen", "consensus50.x", "consensus20.x", "consensus50.y", "consensus20.y", "IFmaxCon.x", "IFmaxCon.y"))]))
rownames(trsp$IF125f)[nchar(rownames(trsp$IF125f)) >= 5] <- c("SsC", "SCID", "DsC", "DCID")

# consider the stress of the NMDS
stress$IF125f <- rep(NA, 10)
for (i in 1:10) {
  stress$IF125f[i] <- metaMDS(daisy(trsp$IF125f), k = i)$stress
}
plot(stress$IF125f, type = "b")  
rm(i)
# looks like between 2 and 3 dimensions would be reasonable

# check for variation
tmp <- metaMDS(daisy(trsp$IF125f))
tmp$stress
plot(tmp, display = "sites", type = "t", cex = 1.5)
# some variation, so optimise

stress$IFop125f <- 1
# find the nmds plot with the lowest stress out of 20000 runs
# for (i in 1:1000) {
#   tmp <- metaMDS(daisy(trsp$IF125f))
#   if (stress$IFop125f > tmp$stress) {
#     nmds$IF125f <- tmp
#     stress$IFop125f <- tmp$stress
#   }
# }
# rm(i, tmp)
# save(stress, nmds, file = "Outputs/NMDS_125.RData")
load("Outputs/NMDS_125.RData")

stressplot(nmds$IF125f)

png("Figures/IF_NMDS_125.png", 600, 600)
plot(nmds$IF125f, type = "n", display = "sites", cex = 1, xlab = "Axis 1", ylab = "Axis 2", las = 1, cex.lab = 1.5, cex.axis = 1.2, main = "NMDS 125")
points(nmds$IF125f, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
text(nmds$IF125f$points[nchar(rownames(nmds$IF125f$points)) <3, ], labels = rownames(nmds$IF125f$points)[nchar(rownames(nmds$IF125f$points)) <3])
points(nmds$IF125f$points[nchar(rownames(nmds$IF125f$points)) >2, ], pch = "+")
text(sweep(data.matrix(nmds$IF125f$points[nchar(rownames(nmds$IF125f$points)) >2 & rownames(nmds$IF125f$points) != "SCID", ]), 2, c(-0.03, 0)), labels = rownames(nmds$IF125f$points)[nchar(rownames(nmds$IF125f$points)) >2 & rownames(nmds$IF125f$points) != "SCID"])
text(nmds$IF125f$points[rownames(nmds$IF125f$points) == "SCID", 1]+0.04, nmds$IF125f$points[rownames(nmds$IF125f$points) == "SCID", 2]+0.01, labels = "SCID")
legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16, cex = 1.3)
dev.off()

# for 150
full.150 <- merge(slide150, digital150, by = "Specimen")
trsp$IF150f <- data.frame(t(full.150[, !(names(full.150) %in% c("Specimen", "consensus50.x", "consensus20.x", "consensus50.y", "consensus20.y", "IFmaxCon.x", "IFmaxCon.y"))]))
rownames(trsp$IF150f)[nchar(rownames(trsp$IF150f)) >= 5] <- c("SsC", "SCID", "DsC", "DCID")

# consider the stress of the NMDS
stress$IF150f <- rep(NA, 10)
for (i in 1:10) {
  stress$IF150f[i] <- metaMDS(daisy(trsp$IF150f), k = i)$stress
}
plot(stress$IF150f, type = "b")
rm(i)
# looks like between 2 and 3 dimensions would be reasonable, although as with Nadia's data, the breakpoint is less obvious for 150 than it is for 125. 

# check for variation
tmp <- metaMDS(daisy(trsp$IF150f))
tmp$stress
plot(tmp, display = "sites", type = "t", cex = 1.5)
# some variation, so optimise

stress$IFop150f <- 1
# find the nmds plot with the lowest stress out of 1000 runs
# for (i in 1:1000) {
#   tmp <- metaMDS(daisy(trsp$IF150f))
#   if (stress$IFop150f > tmp$stress) {
#     nmds$IF150f <- tmp
#     stress$IFop150f <- tmp$stress
#   }
# }
# rm(i, tmp)
# save(stress, nmds, file = "Outputs/NMDS_150f.RData")
load("Outputs/NMDS_150f.RData")



stressplot(nmds$IF150f)

# the full plot
png("Figures/IF_NMDS_150.png", 600, 600)
plot(nmds$IF150f, type = "n", display = "sites", cex = 1, xlab = "Axis 1", ylab = "Axis 2", las = 1, cex.lab = 1.5, cex.axis = 1.2, main = "NMDS 150")
points(nmds$IF150f, pch = 21, cex = 4, col = mds.col$pair, bg = brewer.pal(5, "Set2")[mds.col$sch.col], lwd = 2)
text(nmds$IF150f$points[nchar(rownames(nmds$IF150f$points)) <3, ], labels = rownames(nmds$IF150f$points)[nchar(rownames(nmds$IF150f$points)) <3])
points(nmds$IF150f$points[nchar(rownames(nmds$IF150f$points)) >2, ], pch = "+")
text(sweep(data.matrix(nmds$IF150f$points[nchar(rownames(nmds$IF150f$points)) >2 & rownames(nmds$IF150f$points) != "DsC", ]), 2, c(-0.025, 0)), labels = rownames(nmds$IF150f$points)[nchar(rownames(nmds$IF150f$points)) >2 & rownames(nmds$IF150f$points) != "DsC"])
text(nmds$IF150f$points[rownames(nmds$IF150f$points) == "DsC", 1], nmds$IF150f$points[rownames(nmds$IF150f$points) == "DsC", 2]-0.02, labels = "DsC")
legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16, cex = 1.3)
dev.off()

# focussing on the main section. As noted above, it is better to run this as a new analysis rather than just zoom in, as the influence of the outliers means that the stability of the central points hasn't been tested
nmds$IF150z <- metaMDS(daisy(trsp$IF150f[!(rownames(trsp$IF150f) %in% c("3", "C", "E", "G")),]))

# consider the stress of the NMDS
stress$IF150z <- rep(NA, 10)
for (i in 1:10) {
  stress$IF150z[i] <- metaMDS(daisy(trsp$IF150f[!(rownames(trsp$IF150f) %in% c("3", "C", "E", "G")),]), k = i)$stress
}
plot(stress$IF150z, type = "b")
rm(i)
stressplot(nmds$IF150z)
# again between 2 and 3 dimensions is probably reasonable

# check for variation
tmp <- metaMDS(daisy(trsp$IF150f[!(rownames(trsp$IF150f) %in% c("3", "C", "E", "G")),]))
tmp$stress
plot(tmp, display = "sites", type = "t", cex = 1.5)
# some variation, so optimise

stress$IFop150z <- 1
# find the nmds plot with the lowest stress out of 1000 runs
# for (i in 1:1000) {
#   tmp <- metaMDS(daisy(trsp$IF150f[!(rownames(trsp$IF150f) %in% c("3", "C", "E", "G")),]))
#   if (stress$IFop150z > tmp$stress) {
#     nmds$IF150z <- tmp
#     stress$IFop150z <- tmp$stress
#   }
# }
# rm(i, tmp)
# save(stress, nmds, file = "Outputs/NMDS_150z.RData")
load("Outputs/NMDS_150z.RData")

stress$IFop150z
stressplot(nmds$IF150z)

# the zoomed plot
png("Figures/IF_NMDS_150_zoom.png", 600, 600)
plot(nmds$IF150z, type = "n", display = "sites", cex = 1, xlab = "Axis 1", ylab = "Axis 2", las = 1, cex.lab = 1.5, cex.axis = 1.2, main = "NMDS 150 zoomed")
points(nmds$IF150z, pch = 21, cex = 4, col = mds.col$pair[!(mds.col$person %in% c("3", "C", "E", "G"))], bg = brewer.pal(5, "Set2")[mds.col$sch.col[!(mds.col$person %in% c("3", "C", "E", "G"))]], lwd = 2)
text(nmds$IF150z$points[nchar(rownames(nmds$IF150z$points)) <3, ], labels = rownames(nmds$IF150z$points)[nchar(rownames(nmds$IF150z$points)) <3])
points(nmds$IF150z$points[nchar(rownames(nmds$IF150z$points)) >2, ], pch = "+")
text(sweep(data.matrix(nmds$IF150z$points[nchar(rownames(nmds$IF150z$points)) >2 & rownames(nmds$IF150z$points) != "SsC", ]), 2, c(-0.012, -0.008)), labels = rownames(nmds$IF150z$points)[nchar(rownames(nmds$IF150z$points)) >2 & rownames(nmds$IF150z$points) != "SsC"])
text(nmds$IF150z$points[rownames(nmds$IF150z$points) == "SsC", 1]-0.012, nmds$IF150z$points[rownames(nmds$IF150z$points) == "SsC", 2]-0.008, labels = "SsC")
legend("topright", legend = paste("School", 1:5), col = brewer.pal(5, "Set2")[1:5], pch = 16, cex = 1.3)
dev.off()

# output the scree plots
png("Figures/Scree plots.png")
par(mfrow = c(2, 2))
plot(stress$IF125f, type = "b", bty = "l", las = 1, xlab = "Dimensions", ylab = "Stress", main = "NMDS 125")
plot(stress$IF150f, type = "b", bty = "l", las = 1, xlab = "Dimensions", ylab = "Stress", main = "NMDS 150")
plot(stress$IF150z, type = "b", bty = "l", las = 1, xlab = "Dimensions", ylab = "Stress", main = "NMDS 150 zoomed")
par(mfrow = c(1,1))
dev.off()

# 4b. Dendrogram ----------------------------------------------------------
# an alternative way of plotting this is as a dendrogram
# based on the community
plot(hclust(daisy(data.frame(trsp$IF150f))))

# based on the original IDs
plot(hclust(daisy(data.frame(t(merge(slide125, digital125, by = "Specimen")[, !(names(sp.125) %in% c("species", "IFmaxCon.x", "IFmaxCon.y", "IFc50.x", "IFcMin.x", "IFc50.y", "IFcMin.y"))])))))

# but I don't think this is as helpful

# 5. Repeated analysis by workers -----------------------------------------
# ex figure 6
# again, do this as a confusion matrix

# for 1a / 1b slide 125
head(long$s125)

# I considered different ways of plotting this:
conf_mat(long$s125, "origID", axis.col = "Person", axis1 = "1b", axis2 = "1a", spec.abb = sp.abb, abb.end = c("na", "nc"), key = FALSE)
# or
conf_mat(long$s125, "origID", axis.col = "Person", axis1 = "1a", axis2 = "1b", spec.abb = sp.abb, abb.end = c("na", "nc"), key = FALSE)
# but these only highlight changes, not increasing accuracy, so instead I'm plotting them both against the consensus

sum(slide125$`1a` == slide125$'1b') # 183 or 61% similarity
sum(slide125$`1a` == slide125$IFcMin) # 198 or 66% accuracy
sum(slide125$`1b` == slide125$IFcMin) # 228 or 76% accuracy

png("Figures/Time/confusion_125_1aCon.png", 1000, 700)
conf_mat(long$s125[long$s125$Person == "1a", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "1a", ylab = "Consensus")
dev.off() 
png("Figures/Time/confusion_125_1bCon.png", 1000, 700)
conf_mat(long$s125[long$s125$Person == "1b", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "1b", ylab = "Consensus")
dev.off() 

# and for 2a / 2b
sum(slide125$`2a` == slide125$'2b') # 241 or 80% similarity
sum(slide125$`2a` == slide125$IFcMin) # 208 or 69% accuracy
sum(slide125$`2b` == slide125$IFcMin) # 237 or 79% accuracy

# again, do this as a confusion matrix
png("Figures/Time/confusion_125_2aCon.png", 1000, 700)
conf_mat(long$s125[long$s125$Person == "2a", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2a", ylab = "Consensus")
dev.off() 
png("Figures/Time/confusion_125_2bCon.png", 1000, 700)
conf_mat(long$s125[long$s125$Person == "2b", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2b", ylab = "Consensus")
dev.off() 

# And for 150
sum(slide150$`1a` == slide150$'1b') # 204 or 68% similarity
sum(slide150$`1a` == slide150$IFcMin) # 221 or 74% accuracy
sum(slide150$`1b` == slide150$IFcMin) # 223 or 74% accuracy

png("Figures/Time/confusion_150_1aCon.png", 1000, 700)
conf_mat(long$s150[long$s150$Person == "1a", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "1a", ylab = "Consensus")
dev.off() 
png("Figures/Time/confusion_150_1bCon.png", 1000, 700)
conf_mat(long$s150[long$s150$Person == "1b", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "1b", ylab = "Consensus")
dev.off() 

# 2a/2b
sum(slide150$`2a` == slide150$'2b') # 288 or 96% similarity
sum(slide150$`2a` == slide150$IFcMin) # 255 or 85% accuracy
sum(slide150$`2b` == slide150$IFcMin) # 255 or 85% accuracy

png("Figures/Time/confusion_150_2aCon.png", 1000, 700)
conf_mat(long$s150[long$s150$Person == "2a", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2a", ylab = "Consensus")
dev.off()
png("Figures/Time/confusion_150_2bCon.png", 1000, 700)
conf_mat(long$s150[long$s150$Person == "2b", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2b", ylab = "Consensus")
dev.off()



# 6. Digital vs. slides ---------------------------------------------------

# 6a. Individual comparisons ----------------------------------------------
# Table 6
long$f125 <- rbind(long$s125, long$d125)

# 125 2b vs. A
sum(full.125$`2b` == full.125$'A') # 171 or 57% similarity
sum(full.125$`2b` == full.125$IFcMin.x) # 237 or 79% accuracy
sum(full.125$`A` == full.125$IFcMin.y) # 181 or 60% accuracy

png("Figures/DigitalSlide/confusion_125_2bA.png", 1000, 700)
conf_mat(long$f125, "origID", axis.col = "Person", axis1 = "2b", axis2 = "A", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
dev.off() 

png("Figures/DigitalSlide/confusion_125_2bCon.png", 1000, 700)
conf_mat(long$f125[long$f125$Person == "2b", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2b", ylab = "Consensus")
dev.off() 
png("Figures/DigitalSlide/confusion_125_ACon.png", 1000, 700)
conf_mat(long$f125[long$f125$Person == "A", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "A", ylab = "Consensus")
dev.off() 

# 125 6 vs. F
sum(full.125$`6` == full.125$'F') # 187 or 62% similarity
sum(full.125$`6` == full.125$IFcMin.x) # 211 or 70% accuracy
sum(full.125$`F` == full.125$IFcMin.y) # 226 or 75% accuracy

png("Figures/DigitalSlide/confusion_125_6F.png", 1000, 700)
conf_mat(long$f125, "origID", axis.col = "Person", axis1 = "6", axis2 = "F", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
dev.off() 
png("Figures/DigitalSlide/confusion_125_FCon.png", 1000, 700)
conf_mat(long$f125[long$f125$Person == "F", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "F", ylab = "Consensus")
dev.off() 
png("Figures/DigitalSlide/confusion_125_6Con.png", 1000, 700)
conf_mat(long$f125[long$f125$Person == "6", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "6", ylab = "Consensus")
dev.off() 

# 125 9 vs. G
sum(full.125$`9` == full.125$'G') # 232 or 77% similarity
sum(full.125$`9` == full.125$IFcMin.x) # 206 or 69% accuracy
sum(full.125$`G` == full.125$IFcMin.y) # 157 or 52% accuracy

png("Figures/DigitalSlide/confusion_125_9G.png", 1000, 700)
conf_mat(long$f125, "origID", axis.col = "Person", axis1 = "9", axis2 = "G", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
dev.off() 

png("Figures/DigitalSlide/confusion_125_GCon.png", 1000, 700)
conf_mat(long$f125[long$f125$Person == "G", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "G", ylab = "Consensus")
dev.off() 
png("Figures/DigitalSlide/confusion_125_9Con.png", 1000, 700)
conf_mat(long$f125[long$f125$Person == "9", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "9", ylab = "Consensus")
dev.off() 

# for 150
long$f150 <- rbind(long$s150, long$d150)

# 150 2b vs. A
sum(full.150$`2b` == full.150$'A') # 232 or 77% similarity
sum(full.150$`2b` == full.150$IFcMin.x) # 255 or 85% accuracy
sum(full.150$`A` == full.150$IFcMin.y) # 246 or 82% accuracy

png("Figures/DigitalSlide/confusion_150_2bA.png", 1000, 700)
conf_mat(long$f150, "origID", axis.col = "Person", axis1 = "2b", axis2 = "A", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
dev.off() 

png("Figures/DigitalSlide/confusion_150_2bCon.png", 1000, 700)
conf_mat(long$f150[long$f150$Person == "2b", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "2b", ylab = "Consensus")
dev.off() 
png("Figures/DigitalSlide/confusion_150_ACon.png", 1000, 700)
conf_mat(long$f150[long$f150$Person == "A", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "A", ylab = "Consensus")
dev.off() 

# 150 6 vs. F
sum(full.150$`6` == full.150$'F') # 216 or 72% similarity
sum(full.150$`6` == full.150$IFcMin.x) # 244 or 81% accuracy
sum(full.150$`F` == full.150$IFcMin.y) # 242 or 81% accuracy

png("Figures/DigitalSlide/confusion_150_6F.png", 1000, 700)
conf_mat(long$f150, "origID", axis.col = "Person", axis1 = "6", axis2 = "F", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
dev.off() 

png("Figures/DigitalSlide/confusion_150_FCon.png", 1000, 700)
conf_mat(long$f150[long$f150$Person == "F", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "F", ylab = "Consensus")
dev.off() 
png("Figures/DigitalSlide/confusion_150_6Con.png", 1000, 700)
conf_mat(long$f150[long$f150$Person == "6", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "6", ylab = "Consensus")
dev.off() 

# 150 9 vs. G
sum(full.150$`9` == full.150$'G') # 212 or 71% similarity
sum(full.150$`9` == full.150$IFcMin.x) # 220 or 73% accuracy
sum(full.150$`G` == full.150$IFcMin.y) # 149 or 50% accuracy

png("Figures/DigitalSlide/confusion_150_9G.png", 1000, 700)
conf_mat(long$f150, "origID", axis.col = "Person", axis1 = "9", axis2 = "G", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE)
dev.off() 

png("Figures/DigitalSlide/confusion_150_GCon.png", 1000, 700)
conf_mat(long$f150[long$f150$Person == "G", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "G", ylab = "Consensus")
dev.off() 
png("Figures/DigitalSlide/confusion_150_9Con.png", 1000, 700)
conf_mat(long$f150[long$f150$Person == "9", ], "origID", "IFcMin", spec.abb = sp.abb, abb.end = c("na", "nc"), axes.same = TRUE, xlab = "9", ylab = "Consensus")
dev.off() 

# 6b. Consensus comparisons -----------------------------------------------

sum(full.125$IFcMin.x == full.125$IFcMin.y) # 234 or 78% accuracy
sum(full.150$IFcMin.x == full.150$IFcMin.y) # 248 or 83% accuracy

# plotting the consensus' against each other
png("Figures/DigitalSlide/confusion_125_Con.png", 1000, 700)
conf_mat(long$f125, "IFcMin", axis.col = "Person", axis1 = "A", axis2 = "1a", spec.abb = sp.abb, abb.end = c("na", "nc"), xlab = "Digital", ylab = "Slide")
dev.off()
png("Figures/DigitalSlide/confusion_150_Con.png", 1000, 700)
conf_mat(long$f150, "IFcMin", axis.col = "Person", axis1 = "A", axis2 = "1a", spec.abb = sp.abb, abb.end = c("na", "nc"), xlab = "Digital", ylab = "Slide")
dev.off()

# 7. SST ------------------------------------------------------------------
# this was only done 150 size fraction.
# it is also not perfect as it currently uses Nadia's consensus values not mine. But that's because I can't currently rerun the ANN analysis. 
# Figure 6
head(divTemp)
row.nam <- list()
row.nam$div <- which(nchar(divTemp$Person) < 5)
row.nam$s125 <- which(nchar(divTemp$Person) < 5 & divTemp$Analysis == "Slide" & divTemp$Size == 125)
row.nam$s150 <- which(nchar(divTemp$Person) < 5 & divTemp$Analysis == "Slide" & divTemp$Size == 150)
row.nam$d125 <- which(nchar(divTemp$Person) < 5 & divTemp$Analysis == "Digital" & divTemp$Size == 125)
row.nam$d150 <- which(nchar(divTemp$Person) < 5 & divTemp$Analysis == "Digital" & divTemp$Size == 150)
row.nam$s125c <- which(divTemp$Person == "consensus" & divTemp$Analysis == "Slide" & divTemp$Size == 125)
row.nam$s150c <- which(divTemp$Person == "consensus" & divTemp$Analysis == "Slide" & divTemp$Size == 150)
row.nam$d125c <- which(divTemp$Person == "consensus" & divTemp$Analysis == "Digital" & divTemp$Size == 125)
row.nam$d150c <- which(divTemp$Person == "consensus" & divTemp$Analysis == "Digital" & divTemp$Size == 150)
row.nam$sSST <- which(divTemp$Analysis == "Slide" & divTemp$Size == 150)
row.nam$dSST <- which(divTemp$Analysis == "Digital" & divTemp$Size == 150)

png("Figures/Fig6_SST.png", 900, 800)
par(mfrow = c(2, 1), mar = c(2.5, 4.1, .5, 1))
tmp <- divTemp[row.nam$sSST,]
plot(1:17, tmp$SST10m[tmp$Person != "consensus"], pch = 16, type = "n", xaxt = "n", xlab = "Person", ylab = expression(paste("SST / ", degree, "C")), ylim = c(20, 24))
axis(1, at = 1:17, labels = tmp$Person[tmp$Person != "consensus"])
with(tmp[tmp$Person == "consensus", ], rect(0, SST10m - SD, 18, SST10m + SD, col = rgb(0, .8, .2, alpha = .5)))
abline(h = tmp$SST10m[tmp$Person == "consensus"], lty = 4)
points(1:17, tmp$SST10m[match(accuracySlide$PersonID, tmp$Person)], pch = 16)
with(tmp[match(accuracySlide$PersonID, tmp$Person), ], err_bar(SST10m, SD, 1:17))
with(tmp, lines(c(1:2), c(SST10m[Person == "1a"], SST10m[Person == "1b"])))
with(tmp, lines(c(3:4), c(SST10m[Person == "2a"], SST10m[Person == "2b"])))
text(16.5, 20, "Slide 150", cex = 1.5)
abline(h = 21.76, col = 4)
text(16.5, 21.65, "WOA 1998", cex = 1.3, col = 4)

tmp <- divTemp[row.nam$dSST,]
plot(1:9, tmp$SST10m[tmp$Person != "consensus"], pch = 16, type = "n", xaxt = "n", xlab = "Person", ylab = expression(paste("SST / ", degree, "C")), ylim = c(20, 24))
axis(1, at = 1:9, labels = tmp$Person[tmp$Person != "consensus"])
with(tmp[tmp$Person == "consensus", ], rect(0, SST10m - SD, 18, SST10m + SD, col = rgb(0, .8, .2, alpha = .5)))
abline(h = tmp$SST10m[tmp$Person == "consensus"], lty = 4)
points(1:9, tmp$SST10m[tmp$Person != "consensus"], pch = 16)
with(tmp[tmp$Person != "consensus", ], err_bar(SST10m, SD, 1:9))
text(8.75, 20, "Digital 150", cex = 1.5)
abline(h = 21.76, col = 4)
text(8.75, 21.65, "WOA 1998", cex = 1.3, col = 4)

par(mfrow = c(1,1))
dev.off()
rm(tmp)

png("Figures/Fig6_SST_comb.png", 800, 500)
par(mar = c(5.1, 5.1, 4.1, 2.1))
with(divTemp[divTemp$Size == 150,], plot(1:26, SST10m[match(ord.div, Person)], pch = 16, xaxt = "n", xlab = "Participant", ylab = expression(paste("SST / ", degree, "C")), col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(20, 24), cex.lab = 1.5, las = 1, cex.axis = 1.1))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
with(divTemp[c(row.nam$s150c, row.nam$d150c),], abline(h = c(SST10m - SD, SST10m + SD), col = ((Analysis != "Slide")*3 + 1), lty = 4))
with(divTemp[c(row.nam$s150c, row.nam$d150c),], abline(h = SST10m, col = ((Analysis != "Slide")*3 + 1)))
with(divTemp[divTemp$Size == 150,], err_bar(SST10m[match(ord.div, Person)], SD[match(ord.div, Person)], 1:26, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
abline(h = 21.76, col = "green4")
text(25, 21.65, "WOA 1998", cex = 1.3, col = "green4")
legend("topleft", legend = c("Slide 150", "Digital 150"), pch = 16, col = c(1, 4))
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# 8. Diversity ------------------------------------------------------------

# 8a. Calculating diversity -----------------------------------------------
# Table 5
# compare the results I calculate with those that Nadia has
head(slide125sp)

# these all included 'na' (and 'nc') as a species (to get the match. )
# what do I get if I calculate these without those 
# richness Slide 125
divTemp$IF_Richness <- NA
divTemp$IF_Richness[row.nam$s125] <- specnumber(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))[1:17] # richness
divTemp$IF_Richness[row.nam$s125c] <- specnumber(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))["IFcMin"]

# ShannonWiener
divTemp$IF_ShannonWiener <- NA
divTemp$IF_ShannonWiener[row.nam$s125] <- diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))[1:17] 
divTemp$IF_ShannonWiener[row.nam$s125c] <- diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))["IFcMin"]

# Dominance
divTemp$IF_Dominance <- NA
divTemp$IF_Dominance[row.nam$s125] <- (1 - diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]), index = "simpson"))[1:17]  
divTemp$IF_Dominance[row.nam$s125c] <- (1 - diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]), index = "simpson"))["IFcMin"]  

# Evenness
divTemp$IF_Evenness <- NA
divTemp$IF_Evenness[row.nam$s125] <- (exp(diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))) / specnumber(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)])))[1:17]
divTemp$IF_Evenness[row.nam$s125c] <- (exp(diversity(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)]))) / specnumber(t(slide125sp[3:nrow(slide125sp),2:ncol(slide125sp)])))["IFcMin"]


# Slide 150
# richness 
divTemp$IF_Richness[row.nam$s150] <- specnumber(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]))[1:17] # richness
divTemp$IF_Richness[row.nam$s150c] <- specnumber(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]))["IFcMin"]

# ShannonWiener
divTemp$IF_ShannonWiener[row.nam$s150] <- diversity(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]))[1:17] 
divTemp$IF_ShannonWiener[row.nam$s150c] <- diversity(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]))["IFcMin"]

# Dominance
divTemp$IF_Dominance[row.nam$s150] <- (1 - diversity(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]), index = "simpson"))[1:17]  
divTemp$IF_Dominance[row.nam$s150c] <- (1 - diversity(t(slide150sp[3:nrow(slide150sp),2:ncol(slide150sp)]), index = "simpson"))["IFcMin"]  

# Evenness
divTemp$IF_Evenness[row.nam$s150] <- (exp(diversity(t(slide150sp[3:nrow(slide125sp),2:ncol(slide150sp)]))) / specnumber(t(slide150sp[3:nrow(slide125sp),2:ncol(slide150sp)])))[1:17]
divTemp$IF_Evenness[row.nam$s150c] <- (exp(diversity(t(slide150sp[3:nrow(slide125sp),2:ncol(slide150sp)]))) / specnumber(t(slide150sp[3:nrow(slide125sp),2:ncol(slide150sp)])))["IFcMin"]

# Digital 125
# richness 
divTemp$IF_Richness[row.nam$d125] <- specnumber(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]))[1:9] # richness
divTemp$IF_Richness[row.nam$d125c] <- specnumber(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]))["IFcMin"]

# ShannonWiener
divTemp$IF_ShannonWiener[row.nam$d125] <- diversity(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]))[1:9] 
divTemp$IF_ShannonWiener[row.nam$d125c] <- diversity(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]))["IFcMin"]

# Dominance
divTemp$IF_Dominance[row.nam$d125] <- (1 - diversity(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]), index = "simpson"))[1:9]  
divTemp$IF_Dominance[row.nam$d125c] <- (1 - diversity(t(digital125sp[3:nrow(digital125sp),2:ncol(digital125sp)]), index = "simpson"))["IFcMin"]  

# Evenness
divTemp$IF_Evenness[row.nam$d125] <- (exp(diversity(t(digital125sp[3:nrow(slide125sp),2:ncol(digital125sp)]))) / specnumber(t(digital125sp[3:nrow(slide125sp),2:ncol(digital125sp)])))[1:9]
divTemp$IF_Evenness[row.nam$d125c] <- (exp(diversity(t(digital125sp[3:nrow(slide125sp),2:ncol(digital125sp)]))) / specnumber(t(digital125sp[3:nrow(slide125sp),2:ncol(digital125sp)])))["IFcMin"]

# Digital 150
# richness 
divTemp$IF_Richness[row.nam$d150] <- specnumber(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]))[1:9] # richness
divTemp$IF_Richness[row.nam$d150c] <- specnumber(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]))["IFcMin"]

# ShannonWiener
divTemp$IF_ShannonWiener[row.nam$d150] <- diversity(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]))[1:9] 
divTemp$IF_ShannonWiener[row.nam$d150c] <- diversity(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]))["IFcMin"]

# Dominance
divTemp$IF_Dominance[row.nam$d150] <- (1 - diversity(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]), index = "simpson"))[1:9]  
divTemp$IF_Dominance[row.nam$d150c] <- (1 - diversity(t(digital150sp[3:nrow(digital150sp),2:ncol(digital150sp)]), index = "simpson"))["IFcMin"]  

# Evenness
divTemp$IF_Evenness[row.nam$d150] <- (exp(diversity(t(digital150sp[3:nrow(slide125sp),2:ncol(digital150sp)]))) / specnumber(t(digital150sp[3:nrow(slide125sp),2:ncol(digital150sp)])))[1:9]
divTemp$IF_Evenness[row.nam$d150c] <- (exp(diversity(t(digital150sp[3:nrow(slide125sp),2:ncol(digital150sp)]))) / specnumber(t(digital150sp[3:nrow(slide125sp),2:ncol(digital150sp)])))["IFcMin"]

# 8b. Plotting diversity --------------------------------------------------
ord.div <- c("1a", "1b", "2a", "2b", "A", 3:6, "F", 7:9, "G", 10:15, LETTERS[2:5], LETTERS[8:9])

# Figure 7
png("Figures/Fig7_richness.png", 800, 500)
# 125
par(mar = c(5.1, 5.1, 4.1, 2.1))
with(divTemp[divTemp$Size == 125,], plot(1:26, IF_Richness[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Richness", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(14, 30), cex.lab = 1.5, las = 1, cex.axis = 1.1))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
with(divTemp[c(row.nam$s125c, row.nam$d125c),], abline(h = IF_Richness, lty = 2, col = ((Analysis != "Slide")*3 + 1)))
# 150
with(divTemp[divTemp$Size == 150,], points(1:26, IF_Richness[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
with(divTemp[c(row.nam$s150c, row.nam$d150c),], abline(h = IF_Richness, col = ((Analysis != "Slide")*3 + 1)))
legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("Figures/Fig7_Dominance.png", 800, 500)
# 125
par(mar = c(5.1, 5.1, 4.1, 2.1))
with(divTemp[divTemp$Size == 125,], plot(1:26, IF_Dominance[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Dominance", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(0.1, 0.22), cex.lab = 1.5, las = 1, cex.axis = 1.1))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
with(divTemp[c(row.nam$s125c, row.nam$d125c),], abline(h = IF_Dominance, lty = 2, col = ((Analysis != "Slide")*3 + 1)))
# 150
with(divTemp[divTemp$Size == 150,], points(1:26, IF_Dominance[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
with(divTemp[c(row.nam$s150c, row.nam$d150c),], abline(h = IF_Dominance, col = ((Analysis != "Slide")*3 + 1)))
legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("Figures/Fig7_ShannonWiener.png", 800, 500)
# 125
par(mar = c(5.1, 5.1, 4.1, 2.1))
with(divTemp[divTemp$Size == 125,], plot(1:26, IF_ShannonWiener[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Shannon-Wiener Diversity", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(1.95, 2.6), cex.lab = 1.5, las = 1, cex.axis = 1.1))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
with(divTemp[c(row.nam$s125c, row.nam$d125c),], abline(h = IF_ShannonWiener, lty = 2, col = ((Analysis != "Slide")*3 + 1)))
# 150
with(divTemp[divTemp$Size == 150,], points(1:26, IF_ShannonWiener[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
with(divTemp[c(row.nam$s150c, row.nam$d150c),], abline(h = IF_ShannonWiener, col = ((Analysis != "Slide")*3 + 1)))
legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

png("Figures/Fig7_Evenness.png", 800, 500)
# 125
par(mar = c(5.1, 5.1, 4.1, 2.1))
with(divTemp[divTemp$Size == 125,], plot(1:26, IF_Evenness[match(ord.div, Person)], pch = 1, xaxt = "n", xlab = "Participant", ylab = "Evenness", col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1), ylim = c(0.35, 0.6), cex.lab = 1.5, las = 1, cex.axis = 1.1))
axis(1, at = 1:26, labels = ord.div, cex.axis = 1.1)
with(divTemp[c(row.nam$s125c, row.nam$d125c),], abline(h = IF_Evenness, lty = 2, col = ((Analysis != "Slide")*3 + 1)))
# 150
with(divTemp[divTemp$Size == 150,], points(1:26, IF_Evenness[match(ord.div, Person)], pch = 16, col = ((Analysis[match(ord.div, Person)] != "Slide")*3 + 1)))
with(divTemp[c(row.nam$s150c, row.nam$d150c),], abline(h = IF_Evenness, col = ((Analysis != "Slide")*3 + 1)))
legend("topleft", legend = c("Slide 125", "Slide 150", "Digital 125", "Digital 150"), pch = c(1, 16, 1, 16), col = c(1, 1, 4, 4))
par(mar = c(5.1, 4.1, 4.1, 2.1))
dev.off()

# # 8c. Shannon Wiener comparison -----------------------------------------
png("Figures/Fig8_SWcomp.png")
# if I run this using Nadia's data, I get similar results, however, this should be plotted with my data.
plot(divTemp$IF_ShannonWiener[c(row.nam$s125, row.nam$d125)], divTemp$IF_ShannonWiener[c(row.nam$s150, row.nam$d150)], pch = 16, xlab = "Shannon-Wiener 125", ylab = "Shannon-Wiener 150")
summary(lm(divTemp$IF_ShannonWiener[c(row.nam$s150, row.nam$d150)] ~ divTemp$IF_ShannonWiener[c(row.nam$s125, row.nam$d125)])) # r2 = 0.1558

abline(lm(divTemp$IF_ShannonWiener[c(row.nam$s150, row.nam$d150)] ~ divTemp$IF_ShannonWiener[c(row.nam$s125, row.nam$d125)]))

# Figure 8
points(divTemp$IF_ShannonWiener[divTemp$Size == 125 & divTemp$Person %in% c("10", "13", "15", "B")], divTemp$IF_ShannonWiener[divTemp$Size == 150 & divTemp$Person %in% c("10", "13", "15", "B")], pch = 16, xlab = "Shannon-Wiener 125", ylab = "Shannon-Wiener 150", col = "red")


summary(lm(divTemp$IF_ShannonWiener[divTemp$Size == 150 & !(divTemp$Person %in% c("10", "13", "15", "B"))] ~ divTemp$IF_ShannonWiener[divTemp$Size == 125 & !(divTemp$Person %in% c("10", "13", "15", "B"))]))

abline(lm(divTemp$IF_ShannonWiener[divTemp$Size == 150 & !(divTemp$Person %in% c("10", "13", "15", "B"))] ~ divTemp$IF_ShannonWiener[divTemp$Size == 125 & !(divTemp$Person %in% c("10", "13", "15", "B"))]), lty = 2)
dev.off()

# the other variables show mixed results
# richness r2 = 0.1849, p = 0.0283
# evenness r2 = 0.1891, p = 0.0264
# dominance r2 = 0.0551, p = 0.248

# 9. Outliers -------------------------------------------------------------
# Generate a dataframe of outliers. 
# Table 7

# I decided not to try to do this as it was done here. Instead, I will rank species based on their distance from the mean
outliers <- data.frame(PersonID = accuracyFull$PersonID, Analysis = accuracyFull$Analysis)

# rank for MDS 125
outliers$MDS125 <- NA
tmp.pt <- nmds$IF125f$points["SCID", ]
tmp.tab <- sqrt((nmds$IF125f$points[,1] - tmp.pt[1])^2 + (nmds$IF125f$points[,2] - tmp.pt[2])^2)
outliers$MDS125[outliers$Analysis == "Slide"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))]))
tmp.pt <- nmds$IF125f$points["DCID", ]
tmp.tab <- sqrt((nmds$IF125f$points[,1] - tmp.pt[1])^2 + (nmds$IF125f$points[,2] - tmp.pt[2])^2)
outliers$MDS125[outliers$Analysis == "Digital"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))]))
rm(tmp.pt, tmp.tab)

# and mds 150
outliers$MDS150 <- NA
tmp.pt <- nmds$IF150f$points["SCID", ]
tmp.tab <- sqrt((nmds$IF150f$points[,1] - tmp.pt[1])^2 + (nmds$IF150f$points[,2] - tmp.pt[2])^2)
outliers$MDS150[outliers$Analysis == "Slide"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.tab))]))
tmp.pt <- nmds$IF150f$points["DCID", ]
tmp.tab <- sqrt((nmds$IF150f$points[,1] - tmp.pt[1])^2 + (nmds$IF150f$points[,2] - tmp.pt[2])^2)
outliers$MDS150[outliers$Analysis == "Digital"][order(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))])] <- sort(rank(tmp.tab[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.tab))]))
rm(tmp.pt, tmp.tab)

# for percentage accuracy
outliers$IF_PtAc125 <- NA
outliers$IF_PtAc125[outliers$Analysis == "Slide"][order(100-accuracySlide$IF_PtAc125[match(outliers$PersonID[outliers$Analysis == "Slide"], accuracySlide$PersonID)])] <- sort(rank(100-accuracySlide$IF_PtAc125))
outliers$IF_PtAc125[outliers$Analysis == "Digital"][order(100-accuracyDigital$IF_PtAc125[match(outliers$PersonID[outliers$Analysis == "Digital"], accuracySlide$PersonID)])] <- sort(rank(100-accuracyDigital$IF_PtAc125))

outliers$IF_PtAc150 <- NA
outliers$IF_PtAc150[outliers$Analysis == "Slide"][order(100-accuracySlide$IF_PtAc150[match(outliers$PersonID[outliers$Analysis == "Slide"], accuracyDigital$PersonID)])] <- sort(rank(100-accuracySlide$IF_PtAc150))
outliers$IF_PtAc150[outliers$Analysis == "Digital"][order(100-accuracyDigital$IF_PtAc150[match(outliers$PersonID[outliers$Analysis == "Digital"], accuracyDigital$PersonID)])] <- sort(rank(100-accuracyDigital$IF_PtAc150))

# for mean pairwise agreement
# for percentage accuracy
outliers$mnPA125 <- NA
outliers$mnPA125[outliers$Analysis == "Slide"][order(100-accuracySlide$mnPA125[match(outliers$PersonID[outliers$Analysis == "Slide"], accuracySlide$PersonID)])] <- sort(rank(100-accuracySlide$mnPA125))
outliers$mnPA125[outliers$Analysis == "Digital"][order(100-accuracyDigital$mnPA125[match(outliers$PersonID[outliers$Analysis == "Digital"], accuracySlide$PersonID)])] <- sort(rank(100-accuracyDigital$mnPA125))

outliers$mnPA150 <- NA
outliers$mnPA150[outliers$Analysis == "Slide"][order(100-accuracySlide$mnPA150[match(outliers$PersonID[outliers$Analysis == "Slide"], accuracyDigital$PersonID)])] <- sort(rank(100-accuracySlide$mnPA150))
outliers$mnPA150[outliers$Analysis == "Digital"][order(100-accuracyDigital$mnPA150[match(outliers$PersonID[outliers$Analysis == "Digital"], accuracyDigital$PersonID)])] <- sort(rank(100-accuracyDigital$mnPA150))

# SST
outliers$SST <- NA
tmp.pt <- abs(divTemp$SST10m[divTemp$Analysis == "Slide" & divTemp$Size == 150 & nchar(divTemp$Person) < 5] - 21.76)
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Slide" & divTemp$Size == 150 & nchar(divTemp$Person) < 5]
outliers$SST[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))
tmp.pt <- abs(divTemp$SST10m[divTemp$Analysis == "Digital" & divTemp$Size == 150 & nchar(divTemp$Person) < 5] - 21.76)
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Digital" & divTemp$Size == 150 & nchar(divTemp$Person) < 5]
outliers$SST[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# richness
outliers$Richness125 <- NA
tmp.pt <- abs(divTemp$Richness[divTemp$Analysis == "Slide" & divTemp$Size == 125 & nchar(divTemp$Person) < 5] - divTemp$Richness[divTemp$Analysis == "Slide" & divTemp$Size == 125 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Slide" & divTemp$Size == 125 & nchar(divTemp$Person) < 5]
outliers$Richness125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$Richness[divTemp$Analysis == "Digital" & divTemp$Size == 125 & nchar(divTemp$Person) < 5] - divTemp$Richness[divTemp$Analysis == "Digital" & divTemp$Size == 125 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Digital" & divTemp$Size == 125 & nchar(divTemp$Person) < 5]
outliers$Richness125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$Richness150 <- NA
tmp.pt <- abs(divTemp$Richness[divTemp$Analysis == "Slide" & divTemp$Size == 150 & nchar(divTemp$Person) < 5] - divTemp$Richness[divTemp$Analysis == "Slide" & divTemp$Size == 150 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Slide" & divTemp$Size == 150 & nchar(divTemp$Person) < 5]
outliers$Richness150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$Richness[divTemp$Analysis == "Digital" & divTemp$Size == 150 & nchar(divTemp$Person) < 5] - divTemp$Richness[divTemp$Analysis == "Digital" & divTemp$Size == 150 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Digital" & divTemp$Size == 150 & nchar(divTemp$Person) < 5]
outliers$Richness150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# Dominance
outliers$Dominance125 <- NA
tmp.pt <- abs(divTemp$Dominance[divTemp$Analysis == "Slide" & divTemp$Size == 125 & nchar(divTemp$Person) < 5] - divTemp$Dominance[divTemp$Analysis == "Slide" & divTemp$Size == 125 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Slide" & divTemp$Size == 125 & nchar(divTemp$Person) < 5]
outliers$Dominance125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$Dominance[divTemp$Analysis == "Digital" & divTemp$Size == 125 & nchar(divTemp$Person) < 5] - divTemp$Dominance[divTemp$Analysis == "Digital" & divTemp$Size == 125 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Digital" & divTemp$Size == 125 & nchar(divTemp$Person) < 5]
outliers$Dominance125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$Dominance150 <- NA
tmp.pt <- abs(divTemp$Dominance[divTemp$Analysis == "Slide" & divTemp$Size == 150 & nchar(divTemp$Person) < 5] - divTemp$Dominance[divTemp$Analysis == "Slide" & divTemp$Size == 150 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Slide" & divTemp$Size == 150 & nchar(divTemp$Person) < 5]
outliers$Dominance150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$Dominance[divTemp$Analysis == "Digital" & divTemp$Size == 150 & nchar(divTemp$Person) < 5] - divTemp$Dominance[divTemp$Analysis == "Digital" & divTemp$Size == 150 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Digital" & divTemp$Size == 150 & nchar(divTemp$Person) < 5]
outliers$Dominance150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# Evenness
outliers$Evenness125 <- NA
tmp.pt <- abs(divTemp$Evenness[divTemp$Analysis == "Slide" & divTemp$Size == 125 & nchar(divTemp$Person) < 5] - divTemp$Evenness[divTemp$Analysis == "Slide" & divTemp$Size == 125 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Slide" & divTemp$Size == 125 & nchar(divTemp$Person) < 5]
outliers$Evenness125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$Evenness[divTemp$Analysis == "Digital" & divTemp$Size == 125 & nchar(divTemp$Person) < 5] - divTemp$Evenness[divTemp$Analysis == "Digital" & divTemp$Size == 125 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Digital" & divTemp$Size == 125 & nchar(divTemp$Person) < 5]
outliers$Evenness125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$Evenness150 <- NA
tmp.pt <- abs(divTemp$Evenness[divTemp$Analysis == "Slide" & divTemp$Size == 150 & nchar(divTemp$Person) < 5] - divTemp$Evenness[divTemp$Analysis == "Slide" & divTemp$Size == 150 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Slide" & divTemp$Size == 150 & nchar(divTemp$Person) < 5]
outliers$Evenness150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$Evenness[divTemp$Analysis == "Digital" & divTemp$Size == 150 & nchar(divTemp$Person) < 5] - divTemp$Evenness[divTemp$Analysis == "Digital" & divTemp$Size == 150 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Digital" & divTemp$Size == 150 & nchar(divTemp$Person) < 5]
outliers$Evenness150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# ShannonWiener
outliers$ShannonWiener125 <- NA
tmp.pt <- abs(divTemp$ShannonWiener[divTemp$Analysis == "Slide" & divTemp$Size == 125 & nchar(divTemp$Person) < 5] - divTemp$ShannonWiener[divTemp$Analysis == "Slide" & divTemp$Size == 125 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Slide" & divTemp$Size == 125 & nchar(divTemp$Person) < 5]
outliers$ShannonWiener125[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$ShannonWiener[divTemp$Analysis == "Digital" & divTemp$Size == 125 & nchar(divTemp$Person) < 5] - divTemp$ShannonWiener[divTemp$Analysis == "Digital" & divTemp$Size == 125 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Digital" & divTemp$Size == 125 & nchar(divTemp$Person) < 5]
outliers$ShannonWiener125[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))

outliers$ShannonWiener150 <- NA
tmp.pt <- abs(divTemp$ShannonWiener[divTemp$Analysis == "Slide" & divTemp$Size == 150 & nchar(divTemp$Person) < 5] - divTemp$ShannonWiener[divTemp$Analysis == "Slide" & divTemp$Size == 150 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Slide" & divTemp$Size == 150 & nchar(divTemp$Person) < 5]
outliers$ShannonWiener150[outliers$Analysis == "Slide"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Slide"], names(tmp.pt))])] <- sort(rank(tmp.pt))

tmp.pt <- abs(divTemp$ShannonWiener[divTemp$Analysis == "Digital" & divTemp$Size == 150 & nchar(divTemp$Person) < 5] - divTemp$ShannonWiener[divTemp$Analysis == "Digital" & divTemp$Size == 150 & divTemp$Person == "consensus"])
names(tmp.pt) <- divTemp$Person[divTemp$Analysis == "Digital" & divTemp$Size == 150 & nchar(divTemp$Person) < 5]
outliers$ShannonWiener150[outliers$Analysis == "Digital"][order(tmp.pt[match(outliers$PersonID[outliers$Analysis == "Digital"], names(tmp.pt))])] <- sort(rank(tmp.pt))
rm(tmp.pt)

# weighted sums
outliers$fullSum <- rowSums(outliers[, 3:17])

outliers$wtSum <- rowSums(outliers[, grep("MDS", names(outliers))])/2 + rowSums(outliers[, grep("MPA", names(outliers))])/2 + outliers$SST + rowSums(outliers[, grep("Richness", names(outliers))])/8 + rowSums(outliers[, grep("Dominance", names(outliers))])/8 + rowSums(outliers[, grep("Evenness", names(outliers))])/8 + rowSums(outliers[, grep("ShannonWiener", names(outliers))])/8

outliers$wtSumPt <- outliers$wtSum / c(rep(4*17, 17), rep(4*9, 9)) * 100

# Table 8 
# add percentage of specimens identified to accuracy
accuracySlide$ptID125 <- apply(slide125[,col.nam$s125], 2, function(x) sum(x != "na") / 300 * 100)[accuracySlide$PersonID]
accuracySlide$ptID150 <- apply(slide150[,col.nam$s150], 2, function(x) sum(x != "na") / 300 * 100)[accuracySlide$PersonID]
accuracyDigital$ptID125 <- apply(digital125[,col.nam$d125], 2, function(x) sum(x != "na") / 300 * 100)[accuracyDigital$PersonID]
accuracyDigital$ptID150 <- apply(digital150[,col.nam$d150], 2, function(x) sum(x != "na") / 300 * 100)[accuracyDigital$PersonID]
accuracyFull <- rbind(accuracySlide, accuracyDigital)
accuracyFull$Analysis <- c(rep("Slide", 17), rep("Digital", 9))

# looking for correlations
pairs(outliers[, 3:ncol(outliers)])
pairs(outliers[, grep("125", names(outliers))])
pairs(outliers[, grep("150", names(outliers))])

# Create a .csv file for the accuracy data
# create a subset of the data for this
tmp.sub <- accuracyFull[, c("PersonID", "Analysis", "Experience", "IF_PtAc125", "IF_PtAc150", "l.IF_PtAc125", "l.IF_PtAc150", "ptID125", "ptID150")]
tmp.sub[, grep("125|150", names(tmp.sub))] <- round(tmp.sub[, grep("125|150", names(tmp.sub))], 2)
write.csv(tmp.sub, "Outputs/Accuracy.csv", row.names = FALSE)
rm(tmp.sub)

# 10. Size vs. maximum agreement ------------------------------------------
# Figure 5
png("Figures/Fig5_size_agreement_125.png")
with(size125, plot(slideAgreement, Length, pch = 16, main = "> 125", las = 1, xlab = "Agreement", ylab = expression(paste("Maximum diamter / ", mu, "m"))))
lines(names(tapply(size125$Length, size125$slideAgreement, max)), tapply(size125$Length, size125$slideAgreement, max), pch = 16)
with(size125, points(digitalAgreement, Length, pch = 16, col = "blue"))
lines(names(tapply(size125$Length, size125$digitalAgreement, max)), tapply(size125$Length, size125$digitalAgreement, max), pch = 16, col = 4)
legend("topleft", col = c(1, 4), pch = 16, legend = c("Slide", "Digital"))
dev.off()

png("Figures/Fig5_size_agreement_150.png")
with(size150, plot(slideAgreement, Length, pch = 16, main = "> 150", las = 1, xlab = "Agreement", ylab = expression(paste("Maximum diamter / ", mu, "m"))))
lines(names(tapply(size150$Length, size150$slideAgreement, max)), tapply(size150$Length, size150$slideAgreement, max), pch = 16)
with(size150, points(digitalAgreement, Length, pch = 16, col = "blue"))
lines(names(tapply(size150$Length, size150$digitalAgreement, max)), tapply(size150$Length, size150$digitalAgreement, max), pch = 16, col = 4)
legend("topleft", col = c(1, 4), pch = 16, legend = c("Slide", "Digital"))
dev.off()

# ex. figure 7
# run a  linear model
lm_size_125 <- lm(size125$slideAgreement ~ size125$Length)
par(mfrow = c(2,2))
plot(lm_size_125)
par(mfrow = c(1,1))

# try logging to see if that helps
lm_size_125 <- lm(size125$slideAgreement ~ log(size125$Length))
par(mfrow = c(2,2))
plot(lm_size_125)
par(mfrow = c(1,1))

# slightly better
with(size125, plot(log(Length), slideAgreement, pch = 16))
abline(lm_size_125)
 # however this doesn't really make sense, as it predicts slide agreement above 100%. Given the y-value is bounded between 0 and 100 (or 0 and 1), I think this should be a binomial model. 

glm_size_125 <- glm(slideAgreement/100 ~ log(Length), data = size125, family = "binomial")
par(mfrow = c(2,2))
plot(glm_size_125) # still not a good QQ plot
par(mfrow = c(1,1))

png("Figures/exFig7_sizeAgg125_glm.png")
with(size125, plot(Length, slideAgreement/100, pch = 16))
pred <- predict(glm_size_125, newdata = data.frame(Length = 100:600), type = "response")
points(100:600, pred, type = "l")
dev.off()

glm_size_150 <- glm(slideAgreement/100 ~ log(Length), data = size150, family = "binomial")
par(mfrow = c(2,2))
plot(glm_size_150) # still not a good QQ plot
par(mfrow = c(1,1))

png("Figures/exFig7_sizeAgg150_glm.png")
with(size150, plot(Length, slideAgreement/100, pch = 16))
pred <- predict(glm_size_150, newdata = data.frame(Length = 100:700), type = "response")
points(100:700, pred, type = "l")
dev.off()

# 11. Comparison of different tests ---------------------------------------
# ex. Figure 12
png("Figures/exFig12_Consensus frequency.png", 500, 700)
par(mfrow = c(2, 1))
# for slide
barplot(t(cbind(table(slide125$IFmaxCon), table(slide150$IFmaxCon))), beside = TRUE, xlab = "Maximum consensus", legend.text = c("125", "150"), args.legend = c(x = 7, y = 70), main = "Slide")
abline(v = 18.5, lty = 2)

# for digital
barplot(t(cbind(table(digital125$IFmaxCon), table(digital150$IFmaxCon))), beside = TRUE, xlab = "Maximum consensus", legend.text = c("125", "150"), args.legend = c(x = 4.25, y = 70), main = "Digital")
abline(v = 9.5, lty = 2)
par(mfrow = c(1, 1))
dev.off()

# accuracy as percentage
table(slide125$IFmaxCon)/300 * 100
table(slide150$IFmaxCon)/300 * 100
table(digital125$IFmaxCon)/300 * 100
table(digital150$IFmaxCon)/300 * 100

# fraction with > 50% aggreement
sum(slide125$IFmaxCon > c50_cutoff$slide)/300 * 100 # 72.3%
sum(slide150$IFmaxCon > c50_cutoff$slide)/300 * 100 # 85.3%
sum(digital125$IFmaxCon > c50_cutoff$digital)/300 * 100 # 66.3%
sum(digital150$IFmaxCon > c50_cutoff$digital)/300 * 100 # 76.0%

# sum of unidentified
sum(slide125$IFmaxCon <= c50_cutoff$slide) # 83
sum(slide150$IFmaxCon <= c50_cutoff$slide) # 44
sum(digital125$IFmaxCon <= c50_cutoff$digital) # 101
sum(digital150$IFmaxCon <= c50_cutoff$digital) # 72

# 12. Interpreting the diversity metrics ----------------------------------

# 12a. Comparison with MARGO -----------------------------------------------
# see how the datasets compare with the MARGO values for that latitude
load("../../../../../../../Dropbox/Documents/PhD/Project/MARGO/Outputs/Diversity_measures.RData")

png("Figures/Div_cf_Margo.png", 550, 600)
par(mfrow = c(2,2))
# species richness
plot(ldg.margo.data$Latitude, ldg.margo.data$sp.rich, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Richness")
points(rep(30.2, nrow(divTemp)), divTemp$IF_Richness, col = "blue", pch = 16)
points(rep(30.2, 4), divTemp$IF_Richness[divTemp$Person == "consensus"], col = "red", pch = 16)

# check the metrics are the same
sum(ldg.margo.data$sp.rich != specnumber(ldg.margo.data[,14:44])) # should be 0, so yes

# ShannonWiener
tmp.sw <- diversity(ldg.margo.data[,14:44])
plot(ldg.margo.data$Latitude, tmp.sw, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Shannon Wiener")
points(rep(30.2, nrow(divTemp)), divTemp$IF_ShannonWiener, col = "blue", pch = 16)
points(rep(30.2, 4), divTemp$IF_ShannonWiener[divTemp$Person == "consensus"], col = "red", pch = 16)

# Dominance
tmp.dom <- (1 - diversity(ldg.margo.data[,14:44], index = "simpson"))
plot(ldg.margo.data$Latitude, tmp.dom, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Dominance")
points(rep(30.2, nrow(divTemp)), divTemp$IF_Dominance, col = "blue", pch = 16)
points(rep(30.2, 4), divTemp$IF_Dominance[divTemp$Person == "consensus"], col = "red", pch = 16)


# Evenness
tmp.eve <- (exp(diversity(ldg.margo.data[,14:44])) / specnumber(ldg.margo.data[,14:44]))
plot(ldg.margo.data$Latitude, tmp.eve, pch = 16, col = "grey50", xlab = "Latitude", ylab = "Evenness")
points(rep(30.2, nrow(divTemp)), divTemp$IF_Evenness, col = "blue", pch = 16)
points(rep(30.2, 4), divTemp$IF_Evenness[divTemp$Person == "consensus"], col = "red", pch = 16)
par(mfrow = c(1,1))
dev.off()
rm(tmp.sw, tmp.dom, tmp.eve)


# 12b. Simulations for studying changes -----------------------------------
tmp <- data.frame("a" = c(10, 10, 0, 20, 12, 10, 10, 5, 9, 10, 8, 8, 10, 11), "b" = c(10, 10, 10, 0, 10, 0, 10, 10, 10, 10, 12, 10, 10, 10), "c" = c(2, 2, 2, 2, 0, 12, 4, 2, 2, 1, 2, 4, 3, 1), "d" = c(2, 0, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 1, 2), "e" = c(0, 2, 10, 0, 0, 0, 0, 5, 1, 1, 0, 0, 0, 0))
rownames(tmp) <- c("Orig", "ChangeSpListR", "ChangeSpListC", "LumpCc", "LumpCr", "LumpRc", "LumpRr", "SplitCe", "SplitCu", "SplitR", "BoundShiftCc", "BoundShift2Rc", "BoundShift3Rr", "BoundShift4Cr")
tmp$Rich <- specnumber(tmp[, 1:5]) # richness
tmp$SW <- diversity(tmp[, 1:5]) # Shannon-Wiener
tmp$Dom <- (1 - diversity(tmp[, 1:5], index = "simpson")) # dominance
tmp$Eve <- exp(diversity(tmp[, 1:5])) / specnumber(tmp[, 1:5]) # evenness
tmp
rm(tmp)


# 12c. Direction of change ------------------------------------------------
# add some columns to the diversity dataframe to investigate the magnitude of the change. 
divTemp$Cen_IFE <- divTemp$Cen_IFD <- divTemp$Cen_IFSW <- divTemp$Cen_IFR <- divTemp$Cen_SST <- NA
divTemp[divTemp$Analysis == "Slide" & divTemp$Size == 125 & divTemp$Person != "consensus", grep("Cen", names(divTemp))] <- sweep(data.matrix(divTemp[divTemp$Analysis == "Slide" & divTemp$Size == 125 & divTemp$Person != "consensus", c(8, 10:13)]), 2, as.numeric(divTemp[divTemp$Analysis == "Slide" & divTemp$Size == 125 & divTemp$Person == "consensus", c(8, 10:13)]))

divTemp[divTemp$Analysis == "Slide" & divTemp$Size == 150 & divTemp$Person != "consensus", grep("Cen", names(divTemp))] <- sweep(data.matrix(divTemp[divTemp$Analysis == "Slide" & divTemp$Size == 150 & divTemp$Person != "consensus", c(8, 10:13)]), 2, as.numeric(divTemp[divTemp$Analysis == "Slide" & divTemp$Size == 150 & divTemp$Person == "consensus", c(8, 10:13)]))

divTemp[divTemp$Analysis == "Digital" & divTemp$Size == 125 & divTemp$Person != "consensus", grep("Cen", names(divTemp))] <- sweep(data.matrix(divTemp[divTemp$Analysis == "Digital" & divTemp$Size == 125 & divTemp$Person != "consensus", c(8, 10:13)]), 2, as.numeric(divTemp[divTemp$Analysis == "Digital" & divTemp$Size == 125 & divTemp$Person == "consensus", c(8, 10:13)]))

divTemp[divTemp$Analysis == "Digital" & divTemp$Size == 150 & divTemp$Person != "consensus", grep("Cen", names(divTemp))] <- sweep(data.matrix(divTemp[divTemp$Analysis == "Digital" & divTemp$Size == 150 & divTemp$Person != "consensus", c(8, 10:13)]), 2, as.numeric(divTemp[divTemp$Analysis == "Digital" & divTemp$Size == 150 & divTemp$Person == "consensus", c(8, 10:13)]))

pairs(divTemp[, grep("Cen_", names(divTemp))], col = factor(paste(divTemp$Analysis, divTemp$Size, sep = "_")), pch = 16)

# add some columns to the diversity dataframe to investigate the direction of the change. 
divTemp$Dir_IFE <- divTemp$Dir_IFD <- divTemp$Dir_IFSW <- divTemp$Dir_IFR <- divTemp$Dir_SST <- NA
divTemp[divTemp$Analysis == "Slide" & divTemp$Size == 125 & divTemp$Person != "consensus", grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[divTemp$Analysis == "Slide" & divTemp$Size == 125 & divTemp$Person != "consensus", c(8, 10:13)]), 2, as.numeric(divTemp[divTemp$Analysis == "Slide" & divTemp$Size == 125 & divTemp$Person == "consensus", c(8, 10:13)])) > 0, 1, -1)

divTemp[divTemp$Analysis == "Slide" & divTemp$Size == 150 & divTemp$Person != "consensus", grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[divTemp$Analysis == "Slide" & divTemp$Size == 150 & divTemp$Person != "consensus", c(8, 10:13)]), 2, as.numeric(divTemp[divTemp$Analysis == "Slide" & divTemp$Size == 150 & divTemp$Person == "consensus", c(8, 10:13)])) > 0, 1, -1)

divTemp[divTemp$Analysis == "Digital" & divTemp$Size == 125 & divTemp$Person != "consensus", grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[divTemp$Analysis == "Digital" & divTemp$Size == 125 & divTemp$Person != "consensus", c(8, 10:13)]), 2, as.numeric(divTemp[divTemp$Analysis == "Digital" & divTemp$Size == 125 & divTemp$Person == "consensus", c(8, 10:13)])) > 0, 1, -1)

divTemp[divTemp$Analysis == "Digital" & divTemp$Size == 150 & divTemp$Person != "consensus", grep("Dir", names(divTemp))] <- ifelse(sweep(data.matrix(divTemp[divTemp$Analysis == "Digital" & divTemp$Size == 150 & divTemp$Person != "consensus", c(8, 10:13)]), 2, as.numeric(divTemp[divTemp$Analysis == "Digital" & divTemp$Size == 150 & divTemp$Person == "consensus", c(8, 10:13)])) > 0, 1, -1)

table(divTemp$Dir_SST, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
table(divTemp$Dir_IFR, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
table(divTemp$Dir_IFSW, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
table(divTemp$Dir_IFD, paste(divTemp$Analysis, divTemp$Size, sep = "_"))
table(divTemp$Dir_IFE, paste(divTemp$Analysis, divTemp$Size, sep = "_"))

tmp.div <- divTemp[, grep("Analysis|Size|Person|SST10m|SD|IF_", names(divTemp))]
tmp.div$ID <- paste(tmp.div$Analysis, tmp.div$Person)
head(tmp.div)
tmp.div <- reshape(tmp.div, direction = "wide", v.names = grep("SST10m|SD|IF_", names(tmp.div), value = TRUE), timevar = "Size", idvar = "ID")
tmp.div <- tmp.div[, names(tmp.div) != "ID"]
tmp.div[, grepl("1", names(tmp.div)) & !grepl("Rich", names(tmp.div))] <- round(tmp.div[, grepl("1", names(tmp.div)) & !grepl("Rich", names(tmp.div))], 2)
write.csv(tmp.div, file = "Outputs/Diversity_temperature.csv", row.names = FALSE)
