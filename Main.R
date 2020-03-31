#This is where our main code for our project will go!
DJI <- read.csv("DJI.csv"); head(DJI)

#Create 'Regime' Categories (Regan, Bush, etc.)
DJI["Regime"] <- "None" #Set blank column

#Set the dates for each regime
#Each 'term' is denoted from inaugaration day 
#to the day before the next regime's inaugaration day
dates <- sapply(DJI["Date"], function(x) as.Date(x)) #Convert CSV date character objects 
#to date objects
DJI[dates < as.Date("1989-01-20"), "Regime"] <- "RR"
DJI[dates >= as.Date("1989-01-20") & dates < as.Date("1993-01-20"), "Regime"] <- "GHWB"
DJI[dates >= as.Date("1993-01-20") & dates < as.Date("2001-01-20"), "Regime"] <- "BC"
DJI[dates >= as.Date("2001-01-20") & dates < as.Date("2009-01-20"), "Regime"] <- "GWB"
DJI[dates >= as.Date("2009-01-20") & dates < as.Date("2017-01-20"), "Regime"] <- "BO"
DJI[dates >= as.Date("2017-01-20"), "Regime"] <- "DJT"

#Create a binary category for the political party in control of the White House 
#(Republican = 1, Democrat = 0)
Republican <- (DJI$Regime == "RR") | (DJI$Regime == "GHWB") | (DJI$Regime == "GWB") | (DJI$Regime == "DJT")
DJI <- data.frame(DJI, Republican); head(DJI); sum(Republican); mean(Republican)

#Create a binary category for expansion vs recession year
DJI["Recession"] <- FALSE

#Extract years from date column
DJI.year <- format(as.Date(DJI$Date),"%Y")

#Set the dates
DJI[DJI.year == 1980, "Recession"] <- TRUE
DJI[DJI.year >= 1981 | DJI.year  <= 1981, "Recession"] <- TRUE
DJI[DJI.year >= 1990 | DJI.year  <= 1991, "Recession"] <- TRUE
DJI[DJI.year == 2001, "Recession"] <- TRUE
DJI[DJI.year >= 2007 | DJI.year <= 2009, "Recession"] <- TRUE
DJI[DJI.year >= 2020, "Recession"] <- TRUE

#Set up indices for permutation test by party and president 
index.Republican <- DJI$Republican ; sum(index.Republican) # 4822, matches
index.RR <- (DJI$Regime == "RR") ; sum(index.RR) # 1005
index.GHWB <- (DJI$Regime == "GHWB") ; sum(index.GHWB) # 1011
index.BC <- (DJI$Regime == "BC") ; sum(index.BC) # 2021
index.GWB <- (DJI$Regime == "GWB") ; sum(index.GWB) # 2009
index.BO <- (DJI$Regime == "BO") ; sum(index.BO) # 2015
index.DJT <- (DJI$Regime == "DJT") ; sum(index.DJT) # 797
dim(DJI) ; sum(c(index.RR, index.GHWB, index.BC, index.GWB, index.BO, index.DJT)) # 8858, matches

# Create vector for difference between daily Open values
Open <- DJI$Open
N <- length(Open) ; diffs <- numeric(N - 1)
for (i in 1:N - 1) {
  diffs[i] <- Open[i + 1] - Open[i]
}
head(diffs) ; length(diffs) ; mean(diffs) ; var(diffs) ; sum(diffs) 
# whopping 13760.49 for variance
hist(diffs, breaks = "fd", prob = TRUE)

chg <- numeric(nrow(DJI))
chg[1] <- 0 ; chg[2:nrow(DJI)] <- diffs; chg <- data.frame(chg)
DJI <- data.frame(DJI,chg); head(DJI)

diffs.Republican <- diffs[index.Republican[-1]]
mu.RepDiffs <- mean(diffs.Republican) # -.006959248
diffs.Democrat <- diffs[!(index.Republican[-1])]
mu.DemDiffs <- mean(diffs.Democrat)
mean(diffs.Democrat) # 4.709857
# The means seem different, but given the large variance, this is doubtful.

#Permutation Test
N <- 10^4; result.Republican <- numeric(N); result.Democrat <- numeric(N)
for (i in 1:N) {
  smpl <- sample(index.Republican[-1], replace = FALSE)
  smpl.Republican <- diffs[smpl]
  smpl.Democrat <- diffs[!smpl]
  result.Republican[i] <- mean(smpl.Republican)
  result.Democrat[i] <- mean(smpl.Democrat)
}

#Republican Result
hist(result.Republican, col = "red")
abline(v = mu.RepDiffs, col = "black", lwd = 3)
mu.result.Republican <- mean(result.Republican) ; mu.result.Republican ; mu.RepDiffs
mean(result.Republican <= mu.RepDiffs) 
#2.88% chance of seeing this statistic thus it is statistically significant.

#Democrat Result
hist(result.Democrat, col = "blue")
abline(v = mu.DemDiffs, col = "black", lwd = 3)
mu.result.Democrat <- mean(result.Democrat) ; mu.result.Democrat ; mu.DemDiffs
mean(result.Democrat >= mu.DemDiffs) 
#2.88% chance. Whoa, both means are equally statistically significant.

#Combined Attempt (I'm not sure if I did this right either! So you guys are welcome to check and delete if wrong)
RepAvg <- sum(DJI$chg*(DJI$Republican == TRUE))/sum(DJI$Republican == TRUE) ; RepAvg
DemAvg <- sum(DJI$chg*(DJI$Republican == FALSE))/sum(DJI$Republican == FALSE) ; DemAvg
Obs <-  DemAvg - RepAvg; Obs

N <- 10^4 #number of simulations
diffs <- numeric(N) #this is the vector that will store our simmulated differences
for (i in 1:N) {
  Rep <- sample(DJI$Republican) #This is our permuted party column
  RepMu <- sum(DJI$chg*(Rep == TRUE))/sum(Rep == TRUE) ; RepMu
  DemMu <- sum(DJI$chg*(Rep == FALSE))/sum(Rep == FALSE) ; DemMu
  diffs[i] <- DemMu - RepMu
}
mean(diffs) #inspecting that these are indeed close to zero
hist(diffs, breaks = "FD", probability = TRUE)
abline(v = Obs, col = "red")
pvalue <- (sum(diffs >= Obs) + 1)/(N + 1) ; pvalue 
# 2.71% chance that this extreme of an observed difference would arise by chance .

## Contingency table with chi-square test for political party and recession. 
obs.tbl <- table(DJI$Republican,DJI$Recession) ; obs.tbl
sum(!DJI$Recession)
## Can you replicate the previous contingency table by regime instead?

