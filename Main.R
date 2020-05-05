#Math 23C Term Project by Rakeen Tanvir, Bruno Kömel, and Derl Clausen

#Retrieve any functions we have made for this project
source("prj_Functions.R")
#Load Dow Jones Industrial dataset and prepare it for analysis
DJI <- read.csv("DJI.csv"); head(DJI)
source("prj_DataPreparation.R")
head(DJI)

## Preliminary Data Exploration
#
#Create vector for difference between daily Open values and add as new column to data frame
Open <- DJI$Open
diffs <- diff(DJI$Open) #Get first difference of our data.
head(diffs) # there are some rather large values with double digits before the decimal
diffs.length <- length(diffs) # 8857 as expected
sum(diffs)  # 18975.43 so over double the number of observations
mu.chg.open <- mean(diffs) # 2.142422 around what we expected based on the aforementioned values
med.chg.open <- median(diffs) # 3.4297 is higher than mean, indicating extreme lower values
var.chg.open <- var(diffs) # whopping 13760.49 for variance, we expect this to increase over time
hist(diffs, breaks = "fd", prob = TRUE) 
# resembles normal distribution with narrow concentration around sharp peak at mean 
# and with long tails, with more extreme negative values than positive values
abline(v = mu.chg.open, col = "turquoise")
chg <- numeric(nrow(DJI))
chg[1] <- 0 ; chg[2:nrow(DJI)] <- diffs; chg <- data.frame(chg)
DJI <- data.frame(DJI,chg); head(DJI)

## Linear and logistic regression models
#
# Compare open price first difference and trade volume variables
x <- log(DJI$Volume)
y <- log(abs(DJI$chg))
plot(x~y) # interesting shape to the graph
linmod <- lm(y~x)
summary(linmod) # R-squared value is only .1854, low explanation for residual errors


## Magnitude of First Differences
#
AbsDiffs <- abs(diffs) # the absolute value of the first differences
mu.AbsDiffs <- mean(AbsDiffs) # 69 is the empirical mean
sum(AbsDiffs/diffs.length) # 69 is also the total value of the contributions to the mean
sum(AbsDiffs[AbsDiffs <= mu.AbsDiffs]/diffs.length)/mu.AbsDiffs 
# only 23.9% of the contributions to the mean come from values at or below the mean
sum(AbsDiffs <= mu.AbsDiffs)/diffs.length
# 67.7% of the values are at or below the mean value
# yet 76.1% of the contributions to the mean value come from values above the mean
max.AbsDiffs <- max(AbsDiffs) # 2419.92 is the max value
max(AbsDiffs)/diffs.length / mu.AbsDiffs # single maximum value contributed .4% to the mean
# there are long tails of extreme values with large contributions to the mean
plot(AbsDiffs) # messy
hist(AbsDiffs, breaks = "fd", prob = TRUE) 
# could be modeled by a non-negative valued, long-tailed distribution

## Empirical Cumulative Distributions
#
AbsDiffsCDF <- ecdf(AbsDiffs)
plot(AbsDiffsCDF) # could be modeled by non-negative, long-tailed distribution
#
LogAbsDiffs <- log(AbsDiffs)
plot(LogAbsDiffs) # could be modeled by random walk with positive drift value
#
CDF.diffs <- ecdf(diffs)
plot(CDF.diffs) # logistic regression model could fit
plot(Open) # exponential model could fit
plot(log(Open)) # linear regression or polynomial model could fit


## Hypothesis Testing With Permutation Test
#
#Set up indices for permutation test by party and president 
index.Republican <- DJI$Republican ; sum(index.Republican) # 4822, matches
index.RR <- (DJI$Regime == "RR") ; sum(index.RR) # 1005
index.GHWB <- (DJI$Regime == "GHWB") ; sum(index.GHWB) # 1011
index.BC <- (DJI$Regime == "BC") ; sum(index.BC) # 2021
index.GWB <- (DJI$Regime == "GWB") ; sum(index.GWB) # 2009
index.BO <- (DJI$Regime == "BO") ; sum(index.BO) # 2015
index.DJT <- (DJI$Regime == "DJT") ; sum(index.DJT) # 797
dim(DJI) ; sum(c(index.RR, index.GHWB, index.BC, index.GWB, index.BO, index.DJT)) # 8858, matches
diffs.Republican <- diffs[index.Republican[-1]]
mu.RepDiffs <- mean(diffs.Republican) # -.006959248
diffs.Democrat <- diffs[!(index.Republican[-1])]
mu.DemDiffs <- mean(diffs.Democrat)
mean(diffs.Democrat) # 4.709857
# The means seem different, but given the large variance, this is doubtful.
#
#Permutation Test
N <- 10^4; result.Republican <- numeric(N); result.Democrat <- numeric(N)
for (i in 1:N) {
  smpl <- sample(index.Republican[-1], replace = FALSE)
  smpl.Republican <- diffs[smpl]
  smpl.Democrat <- diffs[!smpl]
  result.Republican[i] <- mean(smpl.Republican)
  result.Democrat[i] <- mean(smpl.Democrat)
}
#
#Republican Result
hist(result.Republican, col = "red")
abline(v = mu.RepDiffs, col = "black", lwd = 3)
mu.result.Republican <- mean(result.Republican) ; mu.result.Republican ; mu.RepDiffs
mean(result.Republican <= mu.RepDiffs) 
#2.83% chance of seeing this statistic thus it is statistically significant.
#
#Democrat Result
hist(result.Democrat, col = "blue")
abline(v = mu.DemDiffs, col = "black", lwd = 3)
mu.result.Democrat <- mean(result.Democrat) ; mu.result.Democrat ; mu.DemDiffs
mean(result.Democrat >= mu.DemDiffs) 
#2.83% chance. Whoa, both means are equally statistically significant.
#
#Rerun as a Combined Permutation Test
RepAvg <- sum(DJI$chg*(DJI$Republican == TRUE))/sum(DJI$Republican == TRUE) ; RepAvg
DemAvg <- sum(DJI$chg*(DJI$Republican == FALSE))/sum(DJI$Republican == FALSE) ; DemAvg
Obs <-  DemAvg - RepAvg; Obs
#
N <- 10^4 #number of simulations
result.Combined <- numeric(N) #this is the vector that will store our simulated differences
for (i in 1:N) {
  Rep <- sample(DJI$Republican) #This is our permuted party column
  RepMu <- sum(DJI$chg*(Rep == TRUE))/sum(Rep == TRUE) ; RepMu
  DemMu <- sum(DJI$chg*(Rep == FALSE))/sum(Rep == FALSE) ; DemMu
  result.Combined[i] <- DemMu - RepMu
}
mean(result.Combined) #inspecting that these are indeed close to zero
hist(result.Combined, breaks = "FD", probability = TRUE, col = "steelblue")
abline(v = Obs, col = "red")
pvalue <-  (sum(result.Combined >= Obs) + 1)/(N + 1) ; pvalue # +1 to counts because of our Observed value
# 3.08% chance that this extreme of an observed difference would arise by chance .

## Hypothesis Testing: Contingency table with chi-square test for political party and recession. 
## 
p <- sum(DJI$Recession)/length(DJI$Recession) # 17.67% of observations are in recession years
obs.tbl <- table(DJI$Republican,DJI$Recession) # Republican has more Recession
exp.tbl <- outer(rowSums(obs.tbl), colSums(obs.tbl))/sum(obs.tbl)
obs.tbl ; exp.tbl
chisq.test(DJI$Republican,DJI$Recession)
# p-value is less than 2.2e-16, far below our .05 threshold, so there would be a very
# small chance that the observed contingency table would arise by chance.
# Thus, the observations provide sufficient evidence to reject the null hypothesis
# that Republican and Democratic regimes are equally likely to be associated with recession
# years from 1985 to early 2020. 
#
#Running this as chi-square test of contingency table including all regimes:
obs.tbl <- table(DJI$Recession, DJI$Regime); obs.tbl #GWB had the most recession days
exp.tbl <- outer(rowSums(obs.tbl), colSums(obs.tbl))/sum(obs.tbl); exp.tbl
chisqvalue <- sum((obs.tbl - exp.tbl)^2/exp.tbl)
pchisq(chisqvalue, df = (2 - 1) * (6 - 1), lower.tail = FALSE) # p-value is 0
#We reject the null hypothesis that recession years arise by chance across regimes. 
#
#Running this as chi-square test specific to each regime with p the observed probabilty of recession:
q <- 1 - p # 0.8233235 probability of not being in a recession
prob <- (DJI$Recession*p + (!DJI$Recession)*q) / sum(DJI$Recession*p + (!DJI$Recession)*q) 
min(prob) ; max(prob) ; sum(prob) 
for (i in unique(DJI$Regime)) {
  print(chisq.test(DJI$Recession, DJI$Regime == i, p = prob))
}
# Null hypothesis is that each regime has the observed probability p of recession across regimes. 
# Note: There could exist carryover/lingering effects of recession or otherwise from one regime to the next..
#Each p-value is less than 2.2e-16, far below the .05 threshold. This indicates
#that no individual regime is equally likely to be associated with recessions
#from the years 1985 to early 2020. (We should look into if our statistical methods
#are correct here.)

# Exploratory Data Analysis of Partial Variance to test for convergence of variance
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
# variance does not converge neither for our Open values nor for our first differences

#Let us now try to fit our DJI data with some other distribution whose properties are familiar:
#Given the lack of convergence in the variance shown above, it appears that we should model the data
#with some distribution with infinite variance. Thus, we will see how well our data is modeled by a 
#Cauchy distribution. 
#
#install.packages("fitdistrplus")
library("fitdistrplus")

fitdist(DJI$chg, distr = "cauchy", method = c("mle"))

hist(DJI$chg, breaks = "FD", probability = TRUE) 

locat <- fitdist(DJI$chg, distr = "cauchy", method = c("mle"))$estimate[1]
shap <- fitdist(DJI$chg, distr = "cauchy", method = c("mle"))$estimate[2]

curve(dcauchy(x, location = locat, scale = shap, log = FALSE), col = "blue", add = TRUE)
n1 <- qcauchy(.1, location = locat, scale = shap, lower.tail = TRUE); n1
pcauchy(n1, location = locat, scale = shap)
mean(DJI$chg <= n1) #about 10.1% of the data is in this bin
abline(v = n1, col = "red")

#Now let's create a vector of deciles so that we can split our data and see if it falls as expected
dec <- qcauchy(seq(0.0,1,by = 0.1), location = locat, scale = shap); dec  #10 binned intervals
Exp <- rep(length(DJI$chg)/10,10); Exp 
binchg <- numeric(10)
for(i in 1:10) {
  binchg[i] <- sum((DJI$chg >= dec[i]) & (DJI$chg <= dec[i + 1] )) ; binchg
}
#Finally we can test for uniformity using a chi-squared test.
ChiStat <- sum((binchg - Exp)^2/Exp); ChiStat #219,8663
#We estimated two parameters (using the sample mean and standard deviation), which costs two degrees of freedom, 
#and we have set the total days to match our sample which costs another so we have 10 - 3 = 7 degrees of freedom
curve(dchisq(x, df = 7), from = 0, to = ChiStat + 5)
abline(v = ChiStat, col = "red") 
pchisq(ChiStat, df = 7, lower.tail = FALSE) # 0
#Given this extremely low chi-square value, it seems that the cauchy distribution is not a good model (at all)
#for the daily fluxes in the Dow Jones Industrial Average. So let's now check how a model with infinite variance
#fits the data. 

#Pareto
library(MASS)
#install.packages("actuar")
library(actuar)

hist(AbsDiffs , breaks = "FD", probability = TRUE) 
shape.pareto <- fitdist(AbsDiffs, "pareto", start=list(shape = 1, scale = 500))$estimate[1]
scale.pareto <- fitdist(AbsDiffs, "pareto", start=list(shape = 1, scale = 500))$estimate[2]
curve(dpareto(x,shape = shape.pareto , scale = scale.pareto ), col = "red", add = TRUE)

n1 <- qpareto(.1, shape = shape.pareto, scale = scale.pareto); n1
ppareto(n1,shape = shape.pareto, scale = scale.pareto)
mean(AbsDiffs <= n1) #about 11% of the data is in this bin
abline(v = n1, col = "blue") #very close to 0, as expected from a pareto distribution

#Now let's create a vector of deciles so that we can split our data and see if it falls as expected
dec <- qpareto(seq(0.0,1,by = 0.1), shape = shape.pareto, scale = scale.pareto); dec  #10 bins
Exp <- rep(length(AbsDiffs)/10,10); Exp 
binabschg <- numeric(10)
for(i in 1:10){
  binabschg[i] <- sum((AbsDiffs >= dec[i]) & (AbsDiffs <= dec[i+1] )) ; binabschg
}
#Finally we can test for uniformity using a chi-squared test.
ChiStatAbs <- sum((binabschg - Exp)^2/Exp); ChiStatAbs #460.8169
#We estimated two parameters (using the sample mean and standard deviation), which costs two degrees of freedom, 
#and we have set the total days to match our sample which costs another so we have 10 - 3 = 7 degrees of freedom
curve(dchisq(x, df = 7), from = 0, to = ChiStatAbs + 5)
abline(v=ChiStatAbs, col = "red") 
pchisq(ChiStatAbs, df = 7, lower.tail = FALSE) # 0
#Given this extremely low chi-square value, it seems that the pareto distribution is not a good model (at all)
#for the absolute daily fluxes in the Dow Jones Industrial Average. 


#Pareto
library(MASS)
#install.packages("actuar")
library(actuar)

hist(AbsDiffs , breaks = "FD", probability = TRUE) 

curve(dpareto(x,shape = shape.pareto , scale = scale.pareto ), col = "red", add = TRUE)

shape.pareto <- fitdist(AbsDiffs, "pareto", start=list(shape = 1, scale = 500))$estimate[1]
scale.pareto <- fitdist(AbsDiffs, "pareto", start=list(shape = 1, scale = 500))$estimate[2]

n1 <- qpareto(.1, shape = shape.pareto, scale = scale.pareto); n1
ppareto(n1,shape = shape.pareto, scale = scale.pareto)
mean(AbsDiffs <= n1) #about 11% of the data is in this bin
abline(v = n1, col = "blue") #very close to 0, as expected from a pareto distribution

#Now let's create a vector of deciles so that we can split our data and see if it falls as expected
dec <- qpareto(seq(0.0,1,by = 0.1), shape = shape.pareto, scale = scale.pareto); dec  #11 bins
Exp <- rep(length(AbsDiffs)/10,10); Exp 
binabschg <- numeric(10)
for(i in 1:10){
  binabschg[i] <- sum((AbsDiffs >= dec[i]) & (AbsDiffs <= dec[i+1] )) ; binabschg
}
#Finally we can test for uniformity using a chi-squared test.
ChiStatAbs <- sum((binabschg - Exp)^2/Exp); ChiStatAbs #460.8169
#We estimated two parameters (using the sample mean and standard deviation), which costs two degrees of freedom, 
#and we have set the total days to match our sample which costs another so we have 10 - 3 = 7 degrees of freedom
curve(dchisq(x, df = 7), from = 0, to = ChiStatAbs + 5)
abline(v=ChiStatAbs, col = "red") 
pchisq(ChiStatAbs, df = 7, lower.tail = FALSE) # 0
#Given this extremely low chi-square value, it seems that the pareto distribution is not a good model (at all)
#for the absolute daily fluxes in the Dow Jones Industrial Average. 

#Taking a look at rare events:
#Let's consider each of the days in our data and examine how far from the mean flux in value each of
#them was. 
mu <- mean(diffs); mu
sigma <- sd(diffs); sigma

N <- length(diffs)
SDs <- numeric(N)
for (i in 1:N){
  SDs[i] <- (diffs[i]-mu)/sigma
}
head(SDs)
head(diffs)

SDs.data <- c(0,SDs); head(SDs.data) 
DJI <- data.frame(DJI,SDs.data)
idx <- which(abs(SDs.data) > 5); head(idx)
unusual <- DJI[idx,]; View(unusual) 
#As we can see, there are 32 days in which the price flux for the Dow Jones was larger than 5 standard 
#deviations away from the mean. To show just how bad of a fit the normal distribution is for our data
#consider the p-values for each of these events. 
N <- nrow(unusual)
pvals <- numeric(N)
for (i in 1:(N)) {
 pvals[i] <- pnorm(abs(unusual$SDs.data[i]*sigma), mean = mu, sd = sigma, lower.tail = FALSE)
}
head(pvals)
rare <- max(pvals) #2.306794e-07, which is pretty much 0. 

#In other words, if the DJI proces followed a normal distribution, we would expect to see the least rare
#of these rare events once in 4,335,021 years.
1/rare

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

## Let's apply simulation methods to the day-over-day (DOD) change in Open values. 
## We can also run this and the previous analyses on the other numerical columns 
## of DJI to see if we arise at similar or different results, 
## once we finish the first round of analysis the DOD change in Open values. 
## If we set up the simulations appropriately, we should expect to see in our results that 
## a greater sample size or a greater number of simulations yields increasingly higher 
## variance, potentially indicated by fat tails when plotting the distribution of the results. 
## Can you figure out a way to leverage a chi-square test or CLT here? 