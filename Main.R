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
head(diffs) ; length(diffs) ; mean(diffs) ; var(diffs) ; sum(diffs) # whopping 13760.49 for variance

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
hist(result.Republican)
abline(v = mu.RepDiffs, col = "red")
mean(result >= mu.RepDiffs) #97% chance of seeing this statistic thus it is not statistically significant.

#Democrat Result
hist(result.Democrat)
abline(v = mu.DemDiffs, col = "red")
mean(result.Democrat >= mu.DemDiffs) #2% chance. Whoa, the mean is statistically significant....

