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
people <- as.data.frame(read_excel("Data/PeopleMetadata.xlsx", na = "NA"))

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
write.csv(slide125[which(slide125$consensus50 != slide125$IFc50), ], "Outputs/Slide125Cmismatch.csv", row.names = FALSE)
slide150[which(slide150$consensus50 != slide150$IFc50), ]
digital125[which(digital125$consensus50 != digital125$IFc50), ] # 2 errors in Nadia's results
digital150[which(digital150$consensus50 != digital150$IFc50), ] # matches for everything

# does removing in NAs help?
tmp <- apply(slide125[, nchar(names(slide125)) < 5], 1, function (x) ifelse(length(which(table(x) > (c50_cutoff_slide - sum(x == "na")))) > 0, names(table(x))[which(table(x) > (c50_cutoff_slide - sum(x == "na")))], "nc")) 
write.csv(slide125[which(slide125$consensus50 != tmp), ], "Outputs/Slide125CmismatchNA.csv", row.names = FALSE)
rm(tmp) # no, it makes it worse

# 2b. Consensus 20 --------------------------------------------------------
# check whether 20% is actually the maximum for the consensus (so what is the highest count per specimen) ignoring na's
slide125$IFmaxCon <- apply(slide125[, nchar(names(slide125)) < 5], 1, function (x)  max(table(x[x != 'na']))) 
slide150$IFmaxCon <- apply(slide150[, nchar(names(slide150)) < 5], 1, function (x) max(table(x[x != 'na']))) 
digital125$IFmaxCon <- apply(digital125[, nchar(names(digital125)) < 5], 1, function (x) max(table(x[x != 'na']))) 
digital150$IFmaxCon <- apply(digital150[, nchar(names(digital150)) < 5], 1, function (x) max(table(x[x != 'na']))) 

# look at the summaries of these
table(slide125$IFmaxCon) # so actually for the slides, it should be consensus-18 (3/17)
table(slide150$IFmaxCon)
table(digital125$IFmaxCon) # and here it would be consensus-22 (2/9)
table(digital150$IFmaxCon)

# calculate consensus-20 
# I'm actually calculating this a the minimum consensus. So take most frequent name - if there are multiple take the first alphabetically.
# where the maximum is 'na', Nadia has ignored that
# I'll see if that is what Nadia has done.
slide125$IFcMin <- apply(slide125[, nchar(names(slide125)) < 5], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][1]) 
slide150$IFcMin <- apply(slide150[, nchar(names(slide150)) < 5], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][1])
digital125$IFcMin <- apply(digital125[, nchar(names(digital125)) < 5], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][1])
digital150$IFcMin <- apply(digital150[, nchar(names(digital150)) < 5], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][1])

# see if this matches Nadia's results
slide125[which(slide125$consensus20 != slide125$IFcMin), ] # I think they are errors
slide150[which(slide150$consensus20 != slide150$IFcMin), ] # I think they are errors
digital125[which(digital125$consensus20 != digital125$IFcMin), ] # pachS vs rub
digital150[which(digital150$consensus20 != digital150$IFcMin), ] # I think these were both errors (assuming ruber pink comes after ruber)

# by alphabetical, was that based on abbreviations or on the full names. It was based on the abbreviations (which is what I am doing)
slide150[order(slide150$IFmaxCon),][1,]
apply(slide150[146, nchar(names(slide150)) < 5], 1, function (x) table(x))
sp.abb[sp.abb$Abbreviation %in% c("cal", "falc"),] # if it were on full names, then consensus name would be falconensis, but it is calida

# was the alphabet only used when there wasn't a single maximum
tmp <- apply(slide125[, nchar(names(slide125)) < 5], 1, function (x) names(table(x))[which(table(x) > 2)][1])
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
  tmp <- table(slide125[,i]) 
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
  size125[, grep("slideCon", names(size125))[i]] <- apply(slide125[, nchar(names(slide125)) < 5], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][i]) 
  size125[, grep("digitalCon", names(size125))[i]] <- apply(digital125[, nchar(names(digital125)) < 5], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][i])
  # and the consensus values for 150
  if (length(grep("slideCon", names(size150))) >= i) 
    size150[, grep("slideCon", names(size150))[i]] <- apply(slide150[, nchar(names(slide150)) < 5], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][i])
  if (length(grep("digitalCon", names(size150))) >= i)
    size150[, grep("digitalCon", names(size150))[i]] <- apply(digital150[, nchar(names(digital150)) < 5], 1, function (x) names(table(x[x != 'na']))[table(x[x != 'na']) == max(table(x[x != 'na']))][i])
}

# add in the specimen agreement values
size125$slideAgreement <- slide125$IFmaxCon/17*100
size150$slideAgreement <- slide150$IFmaxCon/17*100
size125$digitalAgreement <- digital125$IFmaxCon/9*100
size150$digitalAgreement <- digital150$IFmaxCon/9*100

# save it to compare with Nadia's
write.csv(size125, file = "Outputs/Supp3_125.csv", row.names = FALSE)
write.csv(size150, file = "Outputs/Supp3_150.csv", row.names = FALSE)


# 3. Agreement between workers (pairwise comparisons) ---------------------

# 3a. C20_score_participant -----------------------------------------------
accuracySlide <- data.frame(PersonID = names(slide125)[nchar(names(slide125)) < 5][c(1, 3, 2, 4:length(names(slide125)[nchar(names(slide125)) < 5]))], stringsAsFactors = FALSE)
accuracyDigital <- data.frame(PersonID = names(digital125)[nchar(names(digital125)) < 5], stringsAsFactors = FALSE)

# based on Nadia's consensus
accuracySlide$NA_PA125 <- apply(slide125[,nchar(names(slide125)) < 5], 2, function(x) sum(x == slide125$consensus20) / 300 * 100)[accuracySlide$PersonID]
accuracySlide$NA_PA150 <- apply(slide150[,nchar(names(slide150)) < 5], 2, function(x) sum(x == slide150$consensus20) / 300 * 100)[accuracySlide$PersonID]
accuracyDigital$NA_PA125 <- apply(digital125[,nchar(names(digital125)) < 5], 2, function(x) sum(x == digital125$consensus20) / 300 * 100)
accuracyDigital$NA_PA150 <- apply(digital150[,nchar(names(digital150)) < 5], 2, function(x) sum(x == digital150$consensus20) / 300 * 100)
# n.b. these values agree with Nadia's (with the odd rounding error), which is good

# based on my consensus
accuracySlide$IF_PA125 <- apply(slide125[,nchar(names(slide125)) < 5], 2, function(x) sum(x == slide125$IFcMin) / 300 * 100)[accuracySlide$PersonID]
accuracySlide$IF_PA150 <- apply(slide150[,nchar(names(slide150)) < 5], 2, function(x) sum(x == slide150$IFcMin) / 300 * 100)[accuracySlide$PersonID]
accuracyDigital$IF_PA125 <- apply(digital125[,nchar(names(digital125)) < 5], 2, function(x) sum(x == digital125$IFcMin) / 300 * 100)
accuracyDigital$IF_PA150 <- apply(digital150[,nchar(names(digital150)) < 5], 2, function(x) sum(x == digital150$IFcMin) / 300 * 100)

# 3b. C20_score_group -----------------------------------------------------
# the mean value for the group
# Nadia's value
NA_PA125s_group <- mean(accuracySlide$NA_PA125)
NA_PA150s_group <- mean(accuracySlide$NA_PA150)
NA_PA125d_group <- mean(accuracyDigital$NA_PA125)
NA_PA150d_group <- mean(accuracyDigital$NA_PA150)
# My value
IF_PA125s_group <- mean(accuracySlide$IF_PA125)
IF_PA150s_group <- mean(accuracySlide$IF_PA150)
IF_PA125d_group <- mean(accuracyDigital$IF_PA125)
IF_PA150d_group <- mean(accuracyDigital$IF_PA150)

# and sd
NA_PA125s_group_sd <- sd(accuracySlide$NA_PA125)
NA_PA150s_group_sd <- sd(accuracySlide$NA_PA150)
NA_PA125d_group_sd <- sd(accuracyDigital$NA_PA125)
NA_PA150d_group_sd <- sd(accuracyDigital$NA_PA150)
# My value
IF_PA125s_group_sd <- sd(accuracySlide$IF_PA125)
IF_PA150s_group_sd <- sd(accuracySlide$IF_PA150)
IF_PA125d_group_sd <- sd(accuracyDigital$IF_PA125)
IF_PA150d_group_sd <- sd(accuracyDigital$IF_PA150)

# mean consensus value based on doing a full pairwise comparison - these are the same as calculated above, so no need to do them.
sum(slide125$IFcMin == slide125[,nchar(names(slide125)) < 5]) / (300*(sum(nchar(names(slide125)) < 5))) * 100
sum(slide150$IFcMin == slide150[,nchar(names(slide150)) < 5]) / (300*(sum(nchar(names(slide150)) < 5))) * 100
sum(digital125$IFcMin == digital125[,nchar(names(digital125)) < 5]) / (300*(sum(nchar(names(digital125)) < 5))) * 100
sum(digital150$IFcMin == digital150[,nchar(names(digital150)) < 5]) / (300*(sum(nchar(names(digital150)) < 5))) * 100


# 3c. Average pairwise agreement scores -----------------------------------
# mean
accuracySlide$MPA125 <- apply(slide125[,nchar(names(slide125)) < 5], 2, function(x) (sum(x == slide125[,nchar(names(slide125)) < 5]) - 300) / (300*(sum(nchar(names(slide125)) < 5) - 1)) * 100)[accuracySlide$PersonID]
accuracySlide$MPA150 <- apply(slide150[,nchar(names(slide150)) < 5], 2, function(x) (sum(x == slide150[,nchar(names(slide150)) < 5]) - 300) / (300*(sum(nchar(names(slide150)) < 5) - 1)) * 100)[accuracySlide$PersonID]
accuracyDigital$MPA125 <- apply(digital125[,nchar(names(digital125)) < 5], 2, function(x) (sum(x == digital125[,nchar(names(digital125)) < 5]) - 300) / (300*(sum(nchar(names(digital125)) < 5) - 1)) * 100)
accuracyDigital$MPA150 <- apply(digital150[,nchar(names(digital150)) < 5], 2, function(x) (sum(x == digital150[,nchar(names(digital150)) < 5]) - 300) / (300*(sum(nchar(names(digital150)) < 5) - 1)) * 100)
# quite a lot of differences here, though fewer in the digital ones

# SD
accuracySlide$sdPA125 <- NA
accuracySlide$sdPA150 <- NA

for (i in accuracySlide$PersonID) {
  accuracySlide$sdPA125[accuracySlide$PersonID == i] <- sd(apply(slide125[,nchar(names(slide125)) < 5][, i] == slide125[,nchar(names(slide125)) < 5], 2, sum)[which(names(slide125[,nchar(names(slide125)) < 5]) != i)])/300*100
  accuracySlide$sdPA150[accuracySlide$PersonID == i] <- sd(apply(slide150[,nchar(names(slide150)) < 5][, i] == slide150[,nchar(names(slide150)) < 5], 2, sum)[which(names(slide150[,nchar(names(slide150)) < 5]) != i)])/300*100
}

accuracyDigital$sdPA125 <- NA
accuracyDigital$sdPA150 <- NA

for (i in accuracyDigital$PersonID) {
  accuracyDigital$sdPA125[accuracyDigital$PersonID == i] <- sd(apply(digital125[,nchar(names(digital125)) < 5][, i] == digital125[,nchar(names(digital125)) < 5], 2, sum)[which(names(digital125[,nchar(names(digital125)) < 5]) != i)])/300*100
  accuracyDigital$sdPA150[accuracyDigital$PersonID == i] <- sd(apply(digital150[,nchar(names(digital150)) < 5][, i] == digital150[,nchar(names(digital150)) < 5], 2, sum)[which(names(digital150[,nchar(names(digital150)) < 5]) != i)])/300*100
}

# 3d. Look at this pairwise agreement as plots ----------------------------
# compare these as plots
par(mfrow = c(2, 2))
plot(factor(accuracySlide$PersonID), accuracySlide$NA_PA125, main = "Pairwise Agreement 125", ylim = c(45, 90), pch = 16)
points(factor(accuracySlide$PersonID), accuracySlide$IF_PA125, col = "blue", main = "Pairwise Agreement 125", pch = 16)
abline(h = NA_PA125s_group)
abline(h = NA_PA125s_group - NA_PA125s_group_sd, lty = 2)
abline(h = NA_PA125s_group + NA_PA125s_group_sd, lty = 2)
abline(h = IF_PA125s_group, col = "blue")
abline(h = IF_PA125s_group - IF_PA125s_group_sd, lty = 2, col = "blue")
abline(h = IF_PA125s_group + IF_PA125s_group_sd, lty = 2, col = "blue")

plot(factor(accuracySlide$PersonID), accuracySlide$NA_PA150, main = "Pairwise Agreement 150", ylim = c(45, 90))
points(factor(accuracySlide$PersonID), accuracySlide$IF_PA150, col = "blue", main = "Pairwise Agreement 150", pch = 16)
abline(h = NA_PA150s_group)
abline(h = NA_PA150s_group - NA_PA150s_group_sd, lty = 2)
abline(h = NA_PA150s_group + NA_PA150s_group_sd, lty = 2)
abline(h = IF_PA150s_group, col = "blue")
abline(h = IF_PA150s_group - IF_PA150s_group_sd, lty = 2, col = "blue")
abline(h = IF_PA150s_group + IF_PA150s_group_sd, lty = 2, col = "blue")

plot(factor(accuracyDigital$PersonID), accuracyDigital$NA_PA125, main = "Pairwise Agreement 125", ylim = c(45, 90))
points(factor(accuracyDigital$PersonID), accuracyDigital$IF_PA125, col = "blue", main = "Pairwise Agreement 125", pch = 16)
abline(h = NA_PA125d_group)
abline(h = NA_PA125d_group - NA_PA125d_group_sd, lty = 2)
abline(h = NA_PA125d_group + NA_PA125d_group_sd, lty = 2)
abline(h = IF_PA125d_group, col = "blue")
abline(h = IF_PA125d_group - IF_PA125d_group_sd, lty = 2, col = "blue")
abline(h = IF_PA125d_group + IF_PA125d_group_sd, lty = 2, col = "blue")

plot(factor(accuracyDigital$PersonID), accuracyDigital$NA_PA150, main = "Pairwise Agreement 150", ylim = c(45, 90))
points(factor(accuracyDigital$PersonID), accuracyDigital$IF_PA150, col = "blue", main = "Pairwise Agreement 150", pch = 16)
abline(h = NA_PA150d_group)
abline(h = NA_PA150d_group - NA_PA150d_group_sd, lty = 2)
abline(h = NA_PA150d_group + NA_PA150d_group_sd, lty = 2)
abline(h = IF_PA150d_group, col = "blue")
abline(h = IF_PA150d_group - IF_PA150d_group_sd, lty = 2, col = "blue")
abline(h = IF_PA150d_group + IF_PA150d_group_sd, lty = 2, col = "blue")
par(mfrow = c(1,1))

# recreate figure 3
err_bar <- function(mean, sd, xpos, length = 0.05) {
  for(i in 1:length(mean)) {
    arrows(xpos[i], mean[i] - sd[i], xpos[i], mean[i] + sd[i], angle = 90, code = 3, length = length)
  }
}

accuracySlide$Experience <- people$ExperienceSlideA[match(accuracySlide$PersonID, people$SlideID)]
accuracySlide$Experience[accuracySlide$PersonID %in% c("1a", "2a")] <- people$ExperienceSlideA[1:2]
accuracySlide$Experience[accuracySlide$PersonID %in% c("1b", "2b")] <- people$ExperienceSlideB[1:2]
accuracyDigital$Experience <- people$ExperienceDigital[match(accuracyDigital$PersonID, people$DigitalID)]

png("Figures/Fig3_Pairwise_agreement.png", 500, 800)
par(mfrow = c(4, 1), mar = c(2.5, 4.1, .5, 1))
with(accuracySlide, plot(Experience, MPA125, ylim = c(30, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, IF_PA125s_group - IF_PA125s_group_sd, 45, IF_PA125s_group + IF_PA125s_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = IF_PA125s_group, lty = 4)
with(accuracySlide, points(Experience, MPA125, pch = 16))
with(accuracySlide, err_bar(MPA125, sdPA125, Experience))
with(accuracySlide[1, ], text(Experience - 0.5, MPA125, labels = PersonID, cex = 0.7))
with(accuracySlide[2:nrow(accuracySlide), ], text(Experience + 0.5, MPA125, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(MPA125[PersonID == "1a"], MPA125[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(MPA125[PersonID == "2a"], MPA125[PersonID == "2b"])))
text(38, 35, "Slide 125", cex = 1.5)

with(accuracySlide, plot(Experience, MPA150, ylim = c(30, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, IF_PA150s_group - IF_PA150s_group_sd, 45, IF_PA150s_group + IF_PA150s_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = IF_PA150s_group, lty = 4)
with(accuracySlide, points(Experience, MPA150, pch = 16))
with(accuracySlide, err_bar(MPA150, sdPA150, Experience))
with(accuracySlide[1, ], text(Experience - 0.5, MPA150, labels = PersonID, cex = 0.7))
with(accuracySlide[2:nrow(accuracySlide), ], text(Experience + 0.5, MPA150, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(MPA150[PersonID == "1a"], MPA150[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(MPA150[PersonID == "2a"], MPA150[PersonID == "2b"])))
text(38, 35, "Slide 150", cex = 1.5)

with(accuracyDigital, plot(Experience, MPA125, ylim = c(30, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, IF_PA125d_group - IF_PA125d_group_sd, 45, IF_PA125d_group + IF_PA125d_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = IF_PA125d_group, lty = 4)
with(accuracyDigital, points(Experience, MPA125, pch = 16))
with(accuracyDigital, err_bar(MPA125, sdPA125, Experience))
with(accuracyDigital[1, ], text(Experience - 0.5, MPA125, labels = PersonID, cex = 0.7))
with(accuracyDigital[2:nrow(accuracyDigital), ], text(Experience + 0.5, MPA125, labels = PersonID, cex = 0.7))
text(38, 35, "Digital 125", cex = 1.5)

with(accuracyDigital, plot(Experience, MPA150, ylim = c(30, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, IF_PA150d_group - IF_PA150d_group_sd, 45, IF_PA150d_group + IF_PA150d_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = IF_PA150d_group, lty = 4)
with(accuracyDigital, points(Experience, MPA150, pch = 16))
with(accuracyDigital, err_bar(MPA150, sdPA150, Experience))
with(accuracyDigital[1, ], text(Experience - 0.5, MPA150, labels = PersonID, cex = 0.7))
with(accuracyDigital[2:nrow(accuracyDigital), ], text(Experience + 0.5, MPA150, labels = PersonID, cex = 0.7))
text(38, 35, "Digital 150", cex = 1.5)
mfrow = c(1,1)
dev.off()

# is there a more useful way to plot this? What does it actually mean?
# compare consensus PA and MPA
png("Figures/cPA_MPA.png")
plot(1,2, xlim = c(40, 90), ylim = c(38, 74), type = "n", xlab = "consensus PA", ylab = "MPA")
points(accuracySlide$IF_PA125, accuracySlide$MPA125, pch = 4)
abline(lm(accuracySlide$MPA125 ~ accuracySlide$IF_PA125))

points(accuracySlide$NA_PA125, accuracySlide$MPA125, pch = 3)
abline(lm(accuracySlide$MPA125 ~ accuracySlide$NA_PA125), lty = 2)

points(accuracySlide$IF_PA150, accuracySlide$MPA150, pch = 4, col = "red")
abline(lm(accuracySlide$MPA150 ~ accuracySlide$IF_PA150), col = "red")

points(accuracySlide$NA_PA150, accuracySlide$MPA150, pch = 3, col = "red")
abline(lm(accuracySlide$MPA150 ~ accuracySlide$NA_PA150), lty = 2, col = "red")

points(accuracyDigital$IF_PA125, accuracyDigital$MPA125, pch = 4, col = "blue")
abline(lm(accuracyDigital$MPA125 ~ accuracyDigital$IF_PA125), col = "blue")

points(accuracyDigital$NA_PA125, accuracyDigital$MPA125, pch = 3, col = "blue")
abline(lm(accuracyDigital$MPA125 ~ accuracyDigital$NA_PA125), lty = 2, col = "blue")

points(accuracyDigital$IF_PA150, accuracyDigital$MPA150, pch = 4, col = "purple")
abline(lm(accuracyDigital$MPA150 ~ accuracyDigital$IF_PA150), col = "purple")

points(accuracyDigital$NA_PA150, accuracyDigital$MPA150, pch = 3, col = "purple")
abline(lm(accuracyDigital$MPA150 ~ accuracyDigital$NA_PA150), lty = 2, col = "purple")

legend("topleft", legend = c("Slide125", "Slide150", "Digital125", "Digital150", "Isabel (x)", "Nadia (+)"), lty = c(rep(1, 4), 1, 2), col = c("black", "red", "blue", "purple", "black", "black"))
dev.off()

# is it more useful to plot each person relative to the consensus?
png("Figures/Fig3_Consensus_agreement.png", 500, 800)
par(mfrow = c(4, 1), mar = c(2.5, 4.1, .5, 1))
with(accuracySlide, plot(Experience, IF_PA125, ylim = c(45, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, IF_PA125s_group - IF_PA125s_group_sd, 45, IF_PA125s_group + IF_PA125s_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = IF_PA125s_group, lty = 4)
with(accuracySlide, points(Experience, IF_PA125, pch = 16))
with(accuracySlide[1, ], text(Experience - 0.5, IF_PA125, labels = PersonID, cex = 0.7))
with(accuracySlide[2:nrow(accuracySlide), ], text(Experience + 0.5, IF_PA125, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(IF_PA125[PersonID == "1a"], IF_PA125[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(IF_PA125[PersonID == "2a"], IF_PA125[PersonID == "2b"])))
text(38, 35, "Slide 125", cex = 1.5)

with(accuracySlide, plot(Experience, IF_PA150, ylim = c(45, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, IF_PA150s_group - IF_PA150s_group_sd, 45, IF_PA150s_group + IF_PA150s_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = IF_PA150s_group, lty = 4)
with(accuracySlide, points(Experience, IF_PA150, pch = 16))
with(accuracySlide[1, ], text(Experience - 0.5, IF_PA150, labels = PersonID, cex = 0.7))
with(accuracySlide[2:nrow(accuracySlide), ], text(Experience + 0.5, IF_PA150, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(IF_PA150[PersonID == "1a"], IF_PA150[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(IF_PA150[PersonID == "2a"], IF_PA150[PersonID == "2b"])))
text(38, 35, "Slide 150", cex = 1.5)

with(accuracyDigital, plot(Experience, IF_PA125, ylim = c(45, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, IF_PA125d_group - IF_PA125d_group_sd, 45, IF_PA125d_group + IF_PA125d_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = IF_PA125d_group, lty = 4)
with(accuracyDigital, points(Experience, IF_PA125, pch = 16))
with(accuracyDigital[1, ], text(Experience - 0.5, IF_PA125, labels = PersonID, cex = 0.7))
with(accuracyDigital[2:nrow(accuracyDigital), ], text(Experience + 0.5, IF_PA125, labels = PersonID, cex = 0.7))
text(38, 35, "Digital 125", cex = 1.5)

with(accuracyDigital, plot(Experience, IF_PA150, ylim = c(45, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, IF_PA150d_group - IF_PA150d_group_sd, 45, IF_PA150d_group + IF_PA150d_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = IF_PA150d_group, lty = 4)
with(accuracyDigital, points(Experience, IF_PA150, pch = 16))
with(accuracyDigital[1, ], text(Experience - 0.5, IF_PA150, labels = PersonID, cex = 0.7))
with(accuracyDigital[2:nrow(accuracyDigital), ], text(Experience + 0.5, IF_PA150, labels = PersonID, cex = 0.7))
text(38, 35, "Digital 150", cex = 1.5)
mfrow = c(1,1)
dev.off()

# or plotted against Person ID, not experience
png("Figures/Fig3_Consensus_agreement_ID.png", 500, 800)
par(mfrow = c(4, 1), mar = c(2.5, 4.1, .5, 1))
with(accuracySlide, plot(1:17, IF_PA125, ylim = c(45, 90), pch = 16, type = "n", xaxt = "n"))
axis(1, at = 1:17, labels = accuracySlide$PersonID)
rect(0, IF_PA125s_group - IF_PA125s_group_sd, 18, IF_PA125s_group + IF_PA125s_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = IF_PA125s_group, lty = 4)
with(accuracySlide, points(1:17, IF_PA125, pch = 16))
with(accuracySlide, lines(c(1:2), c(IF_PA125[PersonID == "1a"], IF_PA125[PersonID == "1b"])))
with(accuracySlide, lines(c(3:4), c(IF_PA125[PersonID == "2a"], IF_PA125[PersonID == "2b"])))
text(38, 35, "Slide 125", cex = 1.5)

with(accuracySlide, plot(1:17, IF_PA150, ylim = c(45, 90), pch = 16, type = "n", xaxt = "n"))
axis(1, at = 1:17, labels = accuracySlide$PersonID)
rect(-5, IF_PA150s_group - IF_PA150s_group_sd, 45, IF_PA150s_group + IF_PA150s_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = IF_PA150s_group, lty = 4)
with(accuracySlide, points(1:17, IF_PA150, pch = 16))
with(accuracySlide, lines(c(1:2), c(IF_PA150[PersonID == "1a"], IF_PA150[PersonID == "1b"])))
with(accuracySlide, lines(c(3:4), c(IF_PA150[PersonID == "2a"], IF_PA150[PersonID == "2b"])))
text(38, 35, "Slide 150", cex = 1.5)

with(accuracyDigital, plot(1:9, IF_PA125, ylim = c(45, 90), pch = 16, type = "n", xaxt = "n"))
axis(1, at = 1:9, labels = accuracyDigital$PersonID)
rect(-5, IF_PA125d_group - IF_PA125d_group_sd, 45, IF_PA125d_group + IF_PA125d_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = IF_PA125d_group, lty = 4)
with(accuracyDigital, points(1:9, IF_PA125, pch = 16))
text(38, 35, "Digital 125", cex = 1.5)

with(accuracyDigital, plot(1:9, IF_PA150, ylim = c(45, 90), pch = 16, type = "n", xaxt = "n"))
axis(1, at = 1:9, labels = accuracyDigital$PersonID)
rect(-5, IF_PA150d_group - IF_PA150d_group_sd, 45, IF_PA150d_group + IF_PA150d_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = IF_PA150d_group, lty = 4)
with(accuracyDigital, points(1:9, IF_PA150, pch = 16))
mfrow = c(1,1)
dev.off()

# 3e. Lumped pairwise agreement scores ------------------------------------
# Table 4
# what happens to agreements if morphologically similar species are grouped. So:
l.slide125 <- slide125
l.slide150 <- slide150
l.digital125 <- digital125
l.digital150 <- digital150

l.slide125[l.slide125 == "falc"] <- "bull" # bulloides with falconensis
l.slide125[l.slide125 == "cal"] <- "siph"# siphonifera with calida
l.slide125[l.slide125 == "tri"] <- "sac"# sacculifer with trilobus

l.slide150[l.slide150 == "falc"] <- "bull" # bulloides with falconensis
l.slide150[l.slide150 == "cal"] <- "siph"# siphonifera with calida
l.slide150[l.slide150 == "tri"] <- "sac"# sacculifer with trilobus

l.digital125[l.digital125 == "falc"] <- "bull" # bulloides with falconensis
l.digital125[l.digital125 == "cal"] <- "siph"# siphonifera with calida
l.digital125[l.digital125 == "tri"] <- "sac"# sacculifer with trilobus

l.digital150[l.digital150 == "falc"] <- "bull" # bulloides with falconensis
l.digital150[l.digital150 == "cal"] <- "siph"# siphonifera with calida
l.digital150[l.digital150 == "tri"] <- "sac"# sacculifer with trilobus

table(l.slide125$IFcMin)
table(l.slide150$IFcMin)

# now rerun the previous analyses
# based on Nadia's consensus
accuracySlide$l.NA_PA125 <- apply(l.slide125[,nchar(names(l.slide125)) < 5], 2, function(x) sum(x == l.slide125$consensus20) / 300 * 100)[accuracySlide$PersonID]
accuracySlide$l.NA_PA150 <- apply(l.slide150[,nchar(names(l.slide150)) < 5], 2, function(x) sum(x == l.slide150$consensus20) / 300 * 100)[accuracySlide$PersonID]
accuracyDigital$l.NA_PA125 <- apply(l.digital125[,nchar(names(l.digital125)) < 5], 2, function(x) sum(x == l.digital125$consensus20) / 300 * 100)
accuracyDigital$l.NA_PA150 <- apply(l.digital150[,nchar(names(l.digital150)) < 5], 2, function(x) sum(x == l.digital150$consensus20) / 300 * 100)
# n.b. these values agree with Nadia's (with the odd rounding error), which is good

# based on my consensus
accuracySlide$l.IF_PA125 <- apply(l.slide125[,nchar(names(l.slide125)) < 5], 2, function(x) sum(x == l.slide125$IFcMin) / 300 * 100)[accuracySlide$PersonID]
accuracySlide$l.IF_PA150 <- apply(l.slide150[,nchar(names(l.slide150)) < 5], 2, function(x) sum(x == l.slide150$IFcMin) / 300 * 100)[accuracySlide$PersonID]
accuracyDigital$l.IF_PA125 <- apply(l.digital125[,nchar(names(l.digital125)) < 5], 2, function(x) sum(x == l.digital125$IFcMin) / 300 * 100)
accuracyDigital$l.IF_PA150 <- apply(l.digital150[,nchar(names(l.digital150)) < 5], 2, function(x) sum(x == l.digital150$IFcMin) / 300 * 100)

# C20_score_group
# the mean value for the group
# Nadia's value
l.NA_PA125s_group <- mean(accuracySlide$l.NA_PA125)
l.NA_PA150s_group <- mean(accuracySlide$l.NA_PA150)
l.NA_PA125d_group <- mean(accuracyDigital$l.NA_PA125)
l.NA_PA150d_group <- mean(accuracyDigital$l.NA_PA150)
# My value
l.IF_PA125s_group <- mean(accuracySlide$l.IF_PA125)
l.IF_PA150s_group <- mean(accuracySlide$l.IF_PA150)
l.IF_PA125d_group <- mean(accuracyDigital$l.IF_PA125)
l.IF_PA150d_group <- mean(accuracyDigital$l.IF_PA150)

# and sd
l.NA_PA125s_group_sd <- sd(accuracySlide$l.NA_PA125)
l.NA_PA150s_group_sd <- sd(accuracySlide$l.NA_PA150)
l.NA_PA125d_group_sd <- sd(accuracyDigital$l.NA_PA125)
l.NA_PA150d_group_sd <- sd(accuracyDigital$l.NA_PA150)
# My value
l.IF_PA125s_group_sd <- sd(accuracySlide$l.IF_PA125)
l.IF_PA150s_group_sd <- sd(accuracySlide$l.IF_PA150)
l.IF_PA125d_group_sd <- sd(accuracyDigital$l.IF_PA125)
l.IF_PA150d_group_sd <- sd(accuracyDigital$l.IF_PA150)

# mean consensus value based on doing a full pairwise comparison - these are the same as calculated above, so no need to do them.
sum(l.slide125$IFcMin == l.slide125[,nchar(names(l.slide125)) < 5]) / (300*(sum(nchar(names(l.slide125)) < 5))) * 100
sum(l.slide150$IFcMin == l.slide150[,nchar(names(l.slide150)) < 5]) / (300*(sum(nchar(names(l.slide150)) < 5))) * 100
sum(l.digital125$IFcMin == l.digital125[,nchar(names(l.digital125)) < 5]) / (300*(sum(nchar(names(l.digital125)) < 5))) * 100
sum(l.digital150$IFcMin == l.digital150[,nchar(names(l.digital150)) < 5]) / (300*(sum(nchar(names(l.digital150)) < 5))) * 100


# Average pairwise agreement scores
# mean
accuracySlide$l.MPA125 <- apply(l.slide125[,nchar(names(l.slide125)) < 5], 2, function(x) (sum(x == l.slide125[,nchar(names(l.slide125)) < 5]) - 300) / (300*(sum(nchar(names(l.slide125)) < 5) - 1)) * 100)[accuracySlide$PersonID]
accuracySlide$l.MPA150 <- apply(l.slide150[,nchar(names(l.slide150)) < 5], 2, function(x) (sum(x == l.slide150[,nchar(names(l.slide150)) < 5]) - 300) / (300*(sum(nchar(names(l.slide150)) < 5) - 1)) * 100)[accuracySlide$PersonID]
accuracyDigital$l.MPA125 <- apply(l.digital125[,nchar(names(l.digital125)) < 5], 2, function(x) (sum(x == l.digital125[,nchar(names(l.digital125)) < 5]) - 300) / (300*(sum(nchar(names(l.digital125)) < 5) - 1)) * 100)
accuracyDigital$l.MPA150 <- apply(l.digital150[,nchar(names(l.digital150)) < 5], 2, function(x) (sum(x == l.digital150[,nchar(names(l.digital150)) < 5]) - 300) / (300*(sum(nchar(names(l.digital150)) < 5) - 1)) * 100)
# quite a lot of differences here, though fewer in the digital ones

# SD
accuracySlide$l.sdPA125 <- NA
accuracySlide$l.sdPA150 <- NA

for (i in accuracySlide$PersonID) {
  accuracySlide$l.sdPA125[accuracySlide$PersonID == i] <- sd(apply(l.slide125[,nchar(names(l.slide125)) < 5][, i] == l.slide125[,nchar(names(l.slide125)) < 5], 2, sum)[which(names(l.slide125[,nchar(names(l.slide125)) < 5]) != i)])/300*100
  accuracySlide$l.sdPA150[accuracySlide$PersonID == i] <- sd(apply(l.slide150[,nchar(names(l.slide150)) < 5][, i] == l.slide150[,nchar(names(l.slide150)) < 5], 2, sum)[which(names(l.slide150[,nchar(names(l.slide150)) < 5]) != i)])/300*100
}

accuracyDigital$l.sdPA125 <- NA
accuracyDigital$l.sdPA150 <- NA

for (i in accuracyDigital$PersonID) {
  accuracyDigital$l.sdPA125[accuracyDigital$PersonID == i] <- sd(apply(l.digital125[,nchar(names(l.digital125)) < 5][, i] == l.digital125[,nchar(names(l.digital125)) < 5], 2, sum)[which(names(l.digital125[,nchar(names(l.digital125)) < 5]) != i)])/300*100
  accuracyDigital$l.sdPA150[accuracyDigital$PersonID == i] <- sd(apply(l.digital150[,nchar(names(l.digital150)) < 5][, i] == l.digital150[,nchar(names(l.digital150)) < 5], 2, sum)[which(names(l.digital150[,nchar(names(l.digital150)) < 5]) != i)])/300*100
}

# Look at this pairwise agreement as plots
# col.MPAre these as plots
par(mfrow = c(2, 2))
plot(factor(accuracySlide$PersonID), accuracySlide$l.NA_PA125, main = "Pairwise Agreement 125", ylim = c(45, 90), pch = 16)
points(factor(accuracySlide$PersonID), accuracySlide$l.IF_PA125, col = "blue", main = "Pairwise Agreement 125", pch = 16)
abline(h = l.NA_PA125s_group)
abline(h = l.NA_PA125s_group - l.NA_PA125s_group_sd, lty = 2)
abline(h = l.NA_PA125s_group + l.NA_PA125s_group_sd, lty = 2)
abline(h = l.IF_PA125s_group, col = "blue")
abline(h = l.IF_PA125s_group - l.IF_PA125s_group_sd, lty = 2, col = "blue")
abline(h = l.IF_PA125s_group + l.IF_PA125s_group_sd, lty = 2, col = "blue")

plot(factor(accuracySlide$PersonID), accuracySlide$l.NA_PA150, main = "Pairwise Agreement 150", ylim = c(45, 90))
points(factor(accuracySlide$PersonID), accuracySlide$l.IF_PA150, col = "blue", main = "Pairwise Agreement 150", pch = 16)
abline(h = l.NA_PA150s_group)
abline(h = l.NA_PA150s_group - l.NA_PA150s_group_sd, lty = 2)
abline(h = l.NA_PA150s_group + l.NA_PA150s_group_sd, lty = 2)
abline(h = l.IF_PA150s_group, col = "blue")
abline(h = l.IF_PA150s_group - l.IF_PA150s_group_sd, lty = 2, col = "blue")
abline(h = l.IF_PA150s_group + l.IF_PA150s_group_sd, lty = 2, col = "blue")

plot(factor(accuracyDigital$PersonID), accuracyDigital$l.NA_PA125, main = "Pairwise Agreement 125", ylim = c(45, 90))
points(factor(accuracyDigital$PersonID), accuracyDigital$l.IF_PA125, col = "blue", main = "Pairwise Agreement 125", pch = 16)
abline(h = l.NA_PA125d_group)
abline(h = l.NA_PA125d_group - l.NA_PA125d_group_sd, lty = 2)
abline(h = l.NA_PA125d_group + l.NA_PA125d_group_sd, lty = 2)
abline(h = l.IF_PA125d_group, col = "blue")
abline(h = l.IF_PA125d_group - l.IF_PA125d_group_sd, lty = 2, col = "blue")
abline(h = l.IF_PA125d_group + l.IF_PA125d_group_sd, lty = 2, col = "blue")

plot(factor(accuracyDigital$PersonID), accuracyDigital$l.NA_PA150, main = "Pairwise Agreement 150", ylim = c(45, 90))
points(factor(accuracyDigital$PersonID), accuracyDigital$l.IF_PA150, col = "blue", main = "Pairwise Agreement 150", pch = 16)
abline(h = l.NA_PA150d_group)
abline(h = l.NA_PA150d_group - l.NA_PA150d_group_sd, lty = 2)
abline(h = l.NA_PA150d_group + l.NA_PA150d_group_sd, lty = 2)
abline(h = l.IF_PA150d_group, col = "blue")
abline(h = l.IF_PA150d_group - l.IF_PA150d_group_sd, lty = 2, col = "blue")
abline(h = l.IF_PA150d_group + l.IF_PA150d_group_sd, lty = 2, col = "blue")
par(mfrow = c(1,1))

# recreate figure 3
accuracySlide$Experience <- people$ExperienceSlideA[match(accuracySlide$PersonID, people$SlideID)]
accuracySlide$Experience[accuracySlide$PersonID %in% c("1a", "2a")] <- people$ExperienceSlideA[1:2]
accuracySlide$Experience[accuracySlide$PersonID %in% c("1b", "2b")] <- people$ExperienceSlideB[1:2]
accuracyDigital$Experience <- people$ExperienceDigital[match(accuracyDigital$PersonID, people$DigitalID)]

png("Figures/Fig3_Pairwise_agreement_lump.png", 500, 800)
par(mfrow = c(4, 1), mar = c(2.5, 4.1, .5, 1))
with(accuracySlide, plot(Experience, l.MPA125, ylim = c(30, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, l.IF_PA125s_group - l.IF_PA125s_group_sd, 45, l.IF_PA125s_group + l.IF_PA125s_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = l.IF_PA125s_group, lty = 4)
with(accuracySlide, points(Experience, l.MPA125, pch = 16))
with(accuracySlide, err_bar(l.MPA125, l.sdPA125, Experience))
with(accuracySlide[1, ], text(Experience - 0.5, l.MPA125, labels = PersonID, cex = 0.7))
with(accuracySlide[2:nrow(accuracySlide), ], text(Experience + 0.5, l.MPA125, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(l.MPA125[PersonID == "1a"], l.MPA125[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(l.MPA125[PersonID == "2a"], l.MPA125[PersonID == "2b"])))
text(38, 35, "Slide 125", cex = 1.5)

with(accuracySlide, plot(Experience, l.MPA150, ylim = c(30, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, l.IF_PA150s_group - l.IF_PA150s_group_sd, 45, l.IF_PA150s_group + l.IF_PA150s_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = l.IF_PA150s_group, lty = 4)
with(accuracySlide, points(Experience, l.MPA150, pch = 16))
with(accuracySlide, err_bar(l.MPA150, l.sdPA150, Experience))
with(accuracySlide[1, ], text(Experience - 0.5, l.MPA150, labels = PersonID, cex = 0.7))
with(accuracySlide[2:nrow(accuracySlide), ], text(Experience + 0.5, l.MPA150, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(l.MPA150[PersonID == "1a"], l.MPA150[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(l.MPA150[PersonID == "2a"], l.MPA150[PersonID == "2b"])))
text(38, 35, "Slide 150", cex = 1.5)

with(accuracyDigital, plot(Experience, l.MPA125, ylim = c(30, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, l.IF_PA125d_group - l.IF_PA125d_group_sd, 45, l.IF_PA125d_group + l.IF_PA125d_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = l.IF_PA125d_group, lty = 4)
with(accuracyDigital, points(Experience, l.MPA125, pch = 16))
with(accuracyDigital, err_bar(l.MPA125, l.sdPA125, Experience))
with(accuracyDigital[1, ], text(Experience - 0.5, l.MPA125, labels = PersonID, cex = 0.7))
with(accuracyDigital[2:nrow(accuracyDigital), ], text(Experience + 0.5, l.MPA125, labels = PersonID, cex = 0.7))
text(38, 35, "Digital 125", cex = 1.5)

with(accuracyDigital, plot(Experience, l.MPA150, ylim = c(30, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, l.IF_PA150d_group - l.IF_PA150d_group_sd, 45, l.IF_PA150d_group + l.IF_PA150d_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = l.IF_PA150d_group, lty = 4)
with(accuracyDigital, points(Experience, l.MPA150, pch = 16))
with(accuracyDigital, err_bar(l.MPA150, l.sdPA150, Experience))
with(accuracyDigital[1, ], text(Experience - 0.5, l.MPA150, labels = PersonID, cex = 0.7))
with(accuracyDigital[2:nrow(accuracyDigital), ], text(Experience + 0.5, l.MPA150, labels = PersonID, cex = 0.7))
text(38, 35, "Digital 150", cex = 1.5)
mfrow = c(1,1)
dev.off()

# is there a more useful way to plot this? What does it actually mean?
# col.MPAre consensus PA and l.MPA
png("Figures/cPA_l.MPA_lump.png")
plot(1,2, xlim = c(40, 90), ylim = c(38, 74), type = "n", xlab = "consensus PA", ylab = "l.MPA")
points(accuracySlide$l.IF_PA125, accuracySlide$l.MPA125, pch = 4)
abline(lm(accuracySlide$l.MPA125 ~ accuracySlide$l.IF_PA125))

points(accuracySlide$l.NA_PA125, accuracySlide$l.MPA125, pch = 3)
abline(lm(accuracySlide$l.MPA125 ~ accuracySlide$l.NA_PA125), lty = 2)

points(accuracySlide$l.IF_PA150, accuracySlide$l.MPA150, pch = 4, col = "red")
abline(lm(accuracySlide$l.MPA150 ~ accuracySlide$l.IF_PA150), col = "red")

points(accuracySlide$l.NA_PA150, accuracySlide$l.MPA150, pch = 3, col = "red")
abline(lm(accuracySlide$l.MPA150 ~ accuracySlide$l.NA_PA150), lty = 2, col = "red")

points(accuracyDigital$l.IF_PA125, accuracyDigital$l.MPA125, pch = 4, col = "blue")
abline(lm(accuracyDigital$l.MPA125 ~ accuracyDigital$l.IF_PA125), col = "blue")

points(accuracyDigital$l.NA_PA125, accuracyDigital$l.MPA125, pch = 3, col = "blue")
abline(lm(accuracyDigital$l.MPA125 ~ accuracyDigital$l.NA_PA125), lty = 2, col = "blue")

points(accuracyDigital$l.IF_PA150, accuracyDigital$l.MPA150, pch = 4, col = "purple")
abline(lm(accuracyDigital$l.MPA150 ~ accuracyDigital$l.IF_PA150), col = "purple")

points(accuracyDigital$l.NA_PA150, accuracyDigital$l.MPA150, pch = 3, col = "purple")
abline(lm(accuracyDigital$l.MPA150 ~ accuracyDigital$l.NA_PA150), lty = 2, col = "purple")

legend("topleft", legend = c("l.slide125", "l.slide150", "l.digital125", "l.digital150", "Isabel (x)", "Nadia (+)"), lty = c(rep(1, 4), 1, 2), col = c("black", "red", "blue", "purple", "black", "black"))
dev.off()

# is it more useful to plot each person relative to the consensus?
png("Figures/Fig3_Consensus_agreement_lump.png", 500, 800)
par(mfrow = c(4, 1), mar = c(2.5, 4.1, .5, 1))
with(accuracySlide, plot(Experience, l.IF_PA125, ylim = c(45, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, l.IF_PA125s_group - l.IF_PA125s_group_sd, 45, l.IF_PA125s_group + l.IF_PA125s_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = l.IF_PA125s_group, lty = 4)
with(accuracySlide, points(Experience, l.IF_PA125, pch = 16))
with(accuracySlide[1, ], text(Experience - 0.5, l.IF_PA125, labels = PersonID, cex = 0.7))
with(accuracySlide[2:nrow(accuracySlide), ], text(Experience + 0.5, l.IF_PA125, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(l.IF_PA125[PersonID == "1a"], l.IF_PA125[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(l.IF_PA125[PersonID == "2a"], l.IF_PA125[PersonID == "2b"])))
text(38, 35, "Slide 125", cex = 1.5)

with(accuracySlide, plot(Experience, l.IF_PA150, ylim = c(45, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, l.IF_PA150s_group - l.IF_PA150s_group_sd, 45, l.IF_PA150s_group + l.IF_PA150s_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = l.IF_PA150s_group, lty = 4)
with(accuracySlide, points(Experience, l.IF_PA150, pch = 16))
with(accuracySlide[1, ], text(Experience - 0.5, l.IF_PA150, labels = PersonID, cex = 0.7))
with(accuracySlide[2:nrow(accuracySlide), ], text(Experience + 0.5, l.IF_PA150, labels = PersonID, cex = 0.7))
with(accuracySlide, lines(c(Experience[PersonID == "1a"], Experience[PersonID == "1b"]), c(l.IF_PA150[PersonID == "1a"], l.IF_PA150[PersonID == "1b"])))
with(accuracySlide, lines(c(Experience[PersonID == "2a"], Experience[PersonID == "2b"]), c(l.IF_PA150[PersonID == "2a"], l.IF_PA150[PersonID == "2b"])))
text(38, 35, "Slide 150", cex = 1.5)

with(accuracyDigital, plot(Experience, l.IF_PA125, ylim = c(45, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, l.IF_PA125d_group - l.IF_PA125d_group_sd, 45, l.IF_PA125d_group + l.IF_PA125d_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = l.IF_PA125d_group, lty = 4)
with(accuracyDigital, points(Experience, l.IF_PA125, pch = 16))
with(accuracyDigital[1, ], text(Experience - 0.5, l.IF_PA125, labels = PersonID, cex = 0.7))
with(accuracyDigital[2:nrow(accuracyDigital), ], text(Experience + 0.5, l.IF_PA125, labels = PersonID, cex = 0.7))
text(38, 35, "Digital 125", cex = 1.5)

with(accuracyDigital, plot(Experience, l.IF_PA150, ylim = c(45, 90), pch = 16, type = "n", xlim = c(0, 40)))
rect(-5, l.IF_PA150d_group - l.IF_PA150d_group_sd, 45, l.IF_PA150d_group + l.IF_PA150d_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = l.IF_PA150d_group, lty = 4)
with(accuracyDigital, points(Experience, l.IF_PA150, pch = 16))
with(accuracyDigital[1, ], text(Experience - 0.5, l.IF_PA150, labels = PersonID, cex = 0.7))
with(accuracyDigital[2:nrow(accuracyDigital), ], text(Experience + 0.5, l.IF_PA150, labels = PersonID, cex = 0.7))
text(38, 35, "Digital 150", cex = 1.5)
mfrow = c(1,1)
dev.off()

# or plotted against Person ID, not experience
png("Figures/Fig3_Consensus_agreement_ID_lump.png", 500, 800)
par(mfrow = c(4, 1), mar = c(2.5, 4.1, .5, 1))
with(accuracySlide, plot(1:17, l.IF_PA125, ylim = c(45, 90), pch = 16, type = "n", xaxt = "n"))
axis(1, at = 1:17, labels = accuracySlide$PersonID)
rect(0, l.IF_PA125s_group - l.IF_PA125s_group_sd, 18, l.IF_PA125s_group + l.IF_PA125s_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = l.IF_PA125s_group, lty = 4)
with(accuracySlide, points(1:17, l.IF_PA125, pch = 16))
with(accuracySlide, lines(c(1:2), c(l.IF_PA125[PersonID == "1a"], l.IF_PA125[PersonID == "1b"])))
with(accuracySlide, lines(c(3:4), c(l.IF_PA125[PersonID == "2a"], l.IF_PA125[PersonID == "2b"])))
text(38, 35, "Slide 125", cex = 1.5)

with(accuracySlide, plot(1:17, l.IF_PA150, ylim = c(45, 90), pch = 16, type = "n", xaxt = "n"))
axis(1, at = 1:17, labels = accuracySlide$PersonID)
rect(-5, l.IF_PA150s_group - l.IF_PA150s_group_sd, 45, l.IF_PA150s_group + l.IF_PA150s_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = l.IF_PA150s_group, lty = 4)
with(accuracySlide, points(1:17, l.IF_PA150, pch = 16))
with(accuracySlide, lines(c(1:2), c(l.IF_PA150[PersonID == "1a"], l.IF_PA150[PersonID == "1b"])))
with(accuracySlide, lines(c(3:4), c(l.IF_PA150[PersonID == "2a"], l.IF_PA150[PersonID == "2b"])))
text(38, 35, "Slide 150", cex = 1.5)

with(accuracyDigital, plot(1:9, l.IF_PA125, ylim = c(45, 90), pch = 16, type = "n", xaxt = "n"))
axis(1, at = 1:9, labels = accuracyDigital$PersonID)
rect(-5, l.IF_PA125d_group - l.IF_PA125d_group_sd, 45, l.IF_PA125d_group + l.IF_PA125d_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = l.IF_PA125d_group, lty = 4)
with(accuracyDigital, points(1:9, l.IF_PA125, pch = 16))
text(38, 35, "Digital 125", cex = 1.5)

with(accuracyDigital, plot(1:9, l.IF_PA150, ylim = c(45, 90), pch = 16, type = "n", xaxt = "n"))
axis(1, at = 1:9, labels = accuracyDigital$PersonID)
rect(-5, l.IF_PA150d_group - l.IF_PA150d_group_sd, 45, l.IF_PA150d_group + l.IF_PA150d_group_sd, col = rgb(0, .8, .2, alpha = .5))
abline(h = l.IF_PA150d_group, lty = 4)
with(accuracyDigital, points(1:9, l.IF_PA150, pch = 16))
mfrow = c(1,1)
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
    tmp.tab <- rev(sort(table(unlist(head(slide125[which(slide125$IFcMin == tmp.sp),nchar(names(slide125)) < 5])))))
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

# slide150
for(j in sp.abb$Abbreviation) {
  # for each species
  tmp.sp <- j
  # assuming it occurs in the consensus
  tmp.no <- sum(slide150$IFcMin == tmp.sp)
  if (tmp.no > 0) {
    # sorted list of IDs
    tmp.tab <- rev(sort(table(unlist(head(slide150[which(slide150$IFcMin == tmp.sp),nchar(names(slide150)) < 5])))))
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

# digital125
for(j in sp.abb$Abbreviation) {
  # for each species
  tmp.sp <- j
  # assuming it occurs in the consensus
  tmp.no <- sum(digital125$IFcMin == tmp.sp)
  if (tmp.no > 0) {
    # sorted list of IDs
    tmp.tab <- rev(sort(table(unlist(head(digital125[which(digital125$IFcMin == tmp.sp),nchar(names(digital125)) < 5])))))
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

# digital150
for(j in sp.abb$Abbreviation) {
  # for each species
  tmp.sp <- j
  # assuming it occurs in the consensus
  tmp.no <- sum(digital150$IFcMin == tmp.sp)
  if (tmp.no > 0) {
    # sorted list of IDs
    tmp.tab <- rev(sort(table(unlist(head(digital150[which(digital150$IFcMin == tmp.sp),nchar(names(digital150)) < 5])))))
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

# Figure 4

# 3g. Confusion matrix ----------------------------------------------------





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

# 11. Col.MPArison of different tests
# ex. Figure 12
