#Retrieve any functions we have made for this project
source("prj_functions.R")

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
DJI[DJI.year >= 1981 & DJI.year  <= 1981, "Recession"] <- TRUE
DJI[DJI.year >= 1990 & DJI.year  <= 1991, "Recession"] <- TRUE
DJI[DJI.year == 2001, "Recession"] <- TRUE
DJI[DJI.year >= 2007 & DJI.year <= 2009, "Recession"] <- TRUE
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
diffs <- diff(DJI$Open) #Get first difference of our data.
head(diffs) # there are some rather large values with double digits before the decimal
length(diffs) # 8857 as expected
sum(diffs)  # 18975.43 so over double the number of observations
mu.chg.open <- mean(diffs) # 2.142422 around what we expected based on the aforementioned values
med.chg.open <- median(diffs) # 3.4297 is higher than mean, indicating more extreme lesser values
var.chg.open <- var(diffs) # whopping 13760.49 for variance, we expect this to increase over time
hist(diffs, breaks = "fd", prob = TRUE) 
# resembles normal distribution with narrow concentration around sharp peak at mean 
# and with long tails, with more extreme negative values
abline(v = mu.chg.open, col = "turquoise")

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

#Combined Attempt
RepAvg <- sum(DJI$chg*(DJI$Republican == TRUE))/sum(DJI$Republican == TRUE) ; RepAvg
DemAvg <- sum(DJI$chg*(DJI$Republican == FALSE))/sum(DJI$Republican == FALSE) ; DemAvg
Obs <-  DemAvg - RepAvg; Obs

N <- 10^4 #number of simulations
result.Combined <- numeric(N) #this is the vector that will store our simulated differences
for (i in 1:N) {
  Rep <- sample(DJI$Republican) #This is our permuted party column
  RepMu <- sum(DJI$chg*(Rep == TRUE))/sum(Rep == TRUE) ; RepMu
  DemMu <- sum(DJI$chg*(Rep == FALSE))/sum(Rep == FALSE) ; DemMu
  result.Combined[i] <- DemMu - RepMu
}
mean(result.Combined) #inspecting that these are indeed close to zero
hist(result.Combined, breaks = "FD", probability = TRUE)
abline(v = Obs, col = "red")

pvalue <-  (sum(result.Combined >= Obs)+1)/(N+1) ; pvalue
# 2.71% chance that this extreme of an observed difference would arise by chance .

## Contingency table with chi-square test for political party and recession. 
sum(DJI$Recession)/length(DJI$Recession) # 17.67% of observations are in recession years
obs.tbl <- table(DJI$Republican,DJI$Recession) # Republican has more Recession
exp.tbl <- outer(rowSums(obs.tbl), colSums(obs.tbl))/sum(obs.tbl)
obs.tbl ; exp.tbl
chisq.test(DJI$Republican,DJI$Recession)
# p-value is less than 2.2e-16, far below our .05 threshold, so there would be a very
# small chance that the observed contingency table would arise by chance.
# Thus, the observations provide sufficient evidence to reject the null hypothesis
# that Republican and Democratic regimes are equally likely to be associated with recession
# years from 1985 to early 2020. 

## Can you make this more granular by leveraging recession by quarter instead of year?
## Can we do a contingency table that is larger than 2x2 by using regime intead of party?

## Let's apply simulation methods to the day-over-day (DOD) change in Open values. 
## We can also run this and the previous analyses on the other numerical columns 
## of DJI to see if we arise at similar or different results, 
## once we finish the first round of analysis the DOD change in Open values. 
## If we set up the simulations appropriately, we should expect to see in our results that 
## a greater sample size or a greater number of simulations yields increasingly higher 
## variance, potentially indicated by fat tails when plotting the distribution of the results. 
## Can you figure out a way to leverage a chi-square test or CLT here?

#Bruno: I think one good way to show that the data has infinite variance is to show that it cannot be
#modeled with a Normal distribution. (Then we can fit something like a Cauchy distribution later)
#I'm open to taking this out, but I thought we could try to get feedback on it at least!
#So what we are checking is if the daily price fluxes for the Dow Jones have infinite variance
hist(DJI$chg, breaks = "FD", probability = TRUE) #already doesn't look super promising
mu <- mean(DJI$chg)
sigma <- sd(DJI$chg) 
# do we multiply by n the number of sample observations to estimate the population variance?
# maybe not since this is the population? or is it?
curve(dnorm(x,mu,sigma), from = -2500, to = 1000, add = TRUE, col = "red")
n1 <- qnorm(0.1, mu, sigma); n1    #10% of the normal distribution lies below this value
pnorm(n1, mu, sigma)       
mean(DJI$chg <= n1) #6.2%, not great, we would epxect something closer to 10% if it were Normal
#Now let's create a vector of deciles so that we can split our data and see if it falls as expected
dec <- qnorm(seq(0.0,1,by = 0.1), mu,sigma); dec  #11 bins
Exp <- rep(length(DJI$chg)/10,10); Exp 
binchg <- numeric(10)
for(i in 1:10){
  binchg[i] <- sum((DJI$chg >= dec[i]) & (DJI$chg <= dec[i+1] )) ; binchg
}
#Finally we can test for uniformity using a chi-squared test.
ChiStat <- sum((binchg - Exp)^2/Exp); ChiStat #3581.397
#We estimated two parameters (using the sample mean and standard deviation), which costs two degrees of freedom, 
#and we have set the total days to match our sample which costs another so we have 10 - 3 = 7 degrees of freedom
curve(dchisq(x, df = 7), from = 0, to = ChiStat + 5)
abline(v=ChiStat, col = "red") 
pchisq(ChiStat, df = 7, lower.tail = FALSE) # 0
#Given this extremely low chi-square value, it seems that the normal distribution is not a good model (at all)
#for the daily fluxes in the Dow Jones Industrial Average. So let's now check how a model with infinite variance
#fits the data. 

## Generate a RW model with a drift using arima.sim
# Choose one of the two lines of code below, one for mean and one for median of our DOD value changes.
#Get the standard dviation of our differences
sd.diff <- sqrt(mean(diffs^2) - mean(diffs)); sd_diff

#The sequence we will use to cut our data with.
rw.seq <- seq(from = 0, to = max(DJI$Open), by = max(DJI$Open) / 10); rw.seq

#Observed values
rw.obs <- as.vector(table(cut(DJI$Open, breaks = rw.seq))); rw.obs

#Set up Random Walk Expectation for our model
rw.exp <- rep(mean(rw.obs), 10); rw.expect

#ChiSq test for our Observed and Expected data
rw.cs <- ChiSq(rw.obs, rw.exp); rw.cs

#Set up Random Walk Simulation
N <- 10^4; results.RW <- numeric(N)
for(i in 1:N){
rw.drift <- arima.sim(model = list(order = c(0,1,0)), 
                      length(diffs), mean = mu.chg.open,
                      sd = sd_diff)
rw.sim <- as.vector(table(cut(rw.drift, breaks = rw.seq)))
results.RW[i] <- ChiSq(rw.sim, rw.exp)
}

hist(results.RW)
abline(v = rw.cs, col = "red", lwd = 3)

rw.pvalue <- mean(rw.cs >= results.RW); rw.pvalue #.3815
#Fail to reject the null hypothesis. There is about a 40% chance that DJI$Open came
#from a Random Walk.


# Plot rw.drift
plot(DJI$Open, type = "l", xlab = "Time", 
     ylab = "Random Daily Opens", main = "Random Walk Model",
     ylim = c(-10000,50000))

#Graphical volatility test (This is a phrase I made up)
for (i in 1:100){
  rw.drift <- arima.sim(model = list(order = c(0,1,0)), 
                        length(diffs), mean = mu.chg.open,
                        sd = sd_diff)
  lines(rw.drift, col = rgb(runif(1,0,1),runif(1,0,1),runif(1,0,1)))
}
#Extremely volatile -- Infinite variance as time goes forward!!


# Calculate the first difference series rw.drift_diff
N <- length(rw.drift); rw.drift_diff <- numeric(N)
for (i in 2:length(rw.drift)) {
  rw.drift_diff[i] <- rw.drift[i] - rw.drift[i - 1] }
rw.drift_diff # no negative values, does not model DJI DOD changes very well
#Mean of the First Difference Series
mean(rw.drift_diff) # 2.222279 when using mean, 3.263219 when using median
#Plot a histogram of the values and draw the mean/median using abline
hist(rw.drift_diff); abline(v = mean(rw.drift_diff), col = "red", lwd = 3)
#Plot the values of the time series and draw the mean/median using abline - and we get a constant!
plot(rw.drift_diff, type = "l", xlab = "Time", ylab = "Differences in Random Daily Opens",
     main = "First Difference Series"); abline(h = mean(rw.drift_diff), col = "red", lwd = 3)
hist(rw.drift_diff, prob = TRUE)

# Chi Square test for rw.drift_diff, bins by deciles
# Apply central limit theorem, if data is symmetric, should produce standard normal distribution

# Partial variance
N <- length(Open) ; 
variances.normal <- numeric(N - 1)
variances.cauchy <- numeric(N - 1)
variances.Open <- numeric(N - 1)
variances.flux.Open <- numeric(N - 1)
sample.normal <- rnorm(N) ; sample.cauchy <- rcauchy(N)
Open <- DJI$Open ; flux.Open <- DJI$chg
index <- 1:(N - 1)
for (i in 2:N) {
 variances.normal[i - 1] <- var(sample.normal[1:i])
 variances.cauchy[i - 1] <- var(sample.cauchy[1:i])
 variances.Open[i - 1] <- var(Open[1:i])
 variances.flux.Open[i - 1] <- var(flux.Open[1:i])
}
variances.flux.Open <- variances.flux.Open[-1]
par(mfrow = c(2,2))
plot(index,variances.normal, type = "l", col = "steelblue", log = "x")
plot(index,variances.cauchy, type = "l", col = "firebrick", log = "xy")
plot(index,variances.Open, type = "l", col = "yellowgreen", log = "xy")
plot(head(index,-1),variances.flux.Open, type = "l", col = "slategray", log = "xy")
par(mfrow = c(1,1))
summary(variances.normal)
summary(variances.cauchy)
summary(variances.Open)
summary(variances.flux.Open)

# Normal seems to fit random walk, but not first differences (which model fits?) 
# What is a narrow highly concentrated around the mean but with really large variance?
# random walk doesn't seem to have infinite variance, but how does it theoretically?

## Should we perform this analysis for High, Low or Close? Does it change anything?
## Does using a logarithmic or exponential model make the data behave better?
## Given that our observed DOD data produces low frequency, extreme values, particularly negative,
## how do we do a better job of fitting a probability distribution model to our observed data? 
## That is, can we slice our data differently, change the parameters of the random walk model,
## or use a different distribution altogether to model our data? We could use Cauchy, Levy, Pareto, 
## T or other distributions with infinite variance. Are there any that look like our histogram?
## Is this necessary in order to test for infinite variance? 
## Does the current random walk model invalidate our test for infinite variance on the observed data?
## Or can we test for infinite variance without worrying about any model at all? 
## Perhaps partial variance can help here?
## Can we build out the partial variance model and visualization per the shared stackexchange post?

