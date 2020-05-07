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
first.diffs <- diff(DJI$Open) #Get first difference of our data.
diffs <- numeric(nrow(DJI))
diffs[1] <- 0 ; diffs[2:nrow(DJI)] <- first.diffs; head(diffs)
DJI <- data.frame(DJI,diffs) ; head(DJI)
head(diffs) # there are some rather large values with double digits before the decimal
diffs.length <- length(diffs) # 8857 as expected
sum(diffs)  # 18975.43 so over double the number of observations
mu.chg.open <- mean(diffs); mu.chg.open # 2.142422 around what we expected based on the aforementioned values
med.chg.open <- median(diffs); med.chg.open # 3.4297 is higher than mean, indicating extreme lower values
var.chg.open <- var(diffs); var.chg.open # whopping 13760.49 for variance, we expect this to increase over time
hist(diffs, breaks = "fd", prob = TRUE, main = "Histogram of First Differences", xlab = "First Differences") 
# resembles normal distribution with narrow concentration around sharp peak at mean 
# and with long tails, with more extreme negative values than positive values
abline(v = mu.chg.open, col = "turquoise")

## Linear and logistic regression models
#
# Compare open price first difference and trade volume variables
x <- log(DJI$Volume)
y <- log(abs(DJI$diffs)); zidx <- which(y == -Inf); y[zidx] <- 0
plot(y~x, main = "Plot of First Differences by Volume", xlab = "Volume", ylab = "First Differences") # interesting shape to the graph
linmod <- lm(y~x)
summary(linmod) # R-squared value is only .2802, low explanation for residual errors
abline(linmod$coefficients[1], linmod$coefficients[2], col = "red")

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
max.AbsDiffs <- max(AbsDiffs); max.AbsDiffs # 2419.92 is the max value
max(AbsDiffs)/diffs.length / mu.AbsDiffs # single maximum value contributed .4% to the mean
# there are long tails of extreme values with large contributions to the mean
hist(AbsDiffs, breaks = "fd", prob = TRUE) 
# could be modeled by a non-negative valued, long-tailed distribution (see below for Pareto analysis)

## Empirical Cumulative Distributions
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


## Central Limit Theorem, Bootstrap and Empirical Cumulative Distribution
##
## Sources: https://stats.stackexchange.com/questions/2504/test-for-finite-variance 
## and Chihara and Hesterberg's Mathematical Statistics with Resampling
## 
# If our sample is an independent and identically distributed realization 
# from a light-tailed distribution, the central limit theorem should hold.
# Use bootstrap resampling to demonstrate this.
## Bootstrap For A Single Population
N <- 10 ;  n <- 10
meanY <- varY <- Z <- numeric(N) ; X <- diffs ; meanX <- mean(X)
plot(function(x) pnorm(x), xlim = c(-3,3), main = "CDF for a Normal Distribution")
# Perform N boostrap resamples of size n from sample X
for (i in 1:N) {
for (i in 1:N) {
  Y <- sample(X,n,replace = TRUE) # Resample
  meanY[i] <- mean(Y)
  varY[i] <- var(Y)
  Z[i] <- (mean(Y) - meanX) / (sd(Y)/sqrt(n))# Compute Z test statistic
}
  plot(ecdf(Z), add = TRUE, col = "red", lwd = 1, cex = .1)
}
hist(Z, breaks = "fd", prob = TRUE, main = "Histogram of Standardized Random Variable")# approximates standard normal distribution by CLT
hist(varY, breaks = "fd", prob = TRUE, main = "Histogram of First Diffs. Sample Variances") # long tailed distribution with right skew
hist(meanY, breaks = "fd", prob = TRUE, main = "Histogram of First Diffs. Sample Mean") # approximately normal with center around meanX

# Boostrap resampling distribution of sample means has smaller spread than sample
# This is expected since we expect Var(Ybar) = Var(diffs.mu)/n. And as N approaches infinity, 
# the bootstrap sampling distribution of the sample mean of Y should 
# be increasingly normally distributed according to the central limit theorem. 
# It should also approximate the shape and spread of the sampling distribution 
# of sample means from the underlying population,
# though its center will retain the bias for sample mean from the original sample. 
# Next: "perform a large number of bootstraps and compare the empirical distribution 
# function of the observed Z's with the e.d.f. of a N(0,1). A natural way to make this 
# comparison is the Kolmogorov–Smirnov test."


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
mu.RepDiffs <- mean(diffs.Republican) #  0.007355863
diffs.Democrat <- diffs[!(index.Republican[-1])]
mu.DemDiffs <- mean(diffs.Democrat)
mean(diffs.Democrat) # 4.692755
# The means seem different, but given the large variance, this is doubtful.
#Another way to look at the data is to consider the total gains in the DJI during republican and democrat regimes. 
rep.idx <- which(DJI$Republican == TRUE) 
rep.gains <- sum(DJI$diffs[rep.idx]); rep.gains
dem.gains <- sum(DJI$diffs[-rep.idx]); dem.gains
sum(DJI$diffs); sum(rep.gains,dem.gains)
rep.dem.gains <- cbind(rep.gains, dem.gains); rep.dem.gains
library(RColorBrewer)
coul <- brewer.pal(5, "Set2") 
name <- c("Republican Gains", "Democrat Gains" )
barplot(rep.dem.gains, col = coul, main = "Barplot of Cummulative Gains by Party", names = name) #This seems to indicate that maybe the democrat regimes saw more economic prosperity. 

#Let us check with a series of permutation tests:
#First, let us consider if the observed means for republicans and democrats, respectively, are 
#statistically significant:
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
RepAvg <- sum(DJI$diffs*(DJI$Republican == TRUE))/sum(DJI$Republican == TRUE) ; RepAvg
DemAvg <- sum(DJI$diffs*(DJI$Republican == FALSE))/sum(DJI$Republican == FALSE) ; DemAvg
Obs <-  DemAvg - RepAvg; Obs
#
N <- 10^4 #number of simulations
result.Combined <- numeric(N) #this is the vector that will store our simulated differences
for (i in 1:N) {
  Rep <- sample(DJI$Republican) #This is our permuted party column
  RepMu <- sum(DJI$diffs*(Rep == TRUE))/sum(Rep == TRUE) ; RepMu
  DemMu <- sum(DJI$diffs*(Rep == FALSE))/sum(Rep == FALSE) ; DemMu
  result.Combined[i] <- DemMu - RepMu
}
mean(result.Combined) #inspecting that these are indeed close to zero
hist(result.Combined, breaks = "FD", probability = TRUE, col = "steelblue")
abline(v = Obs, col = "red")
pvalue <-  (sum(result.Combined >= Obs) + 1)/(N + 1) ; pvalue # +1 to counts because of our Observed value
# 2.87% chance that this extreme of an observed difference would arise by chance, so it appears that the DJI performed
#better during democratic regimes, a result that is statistically signifficant.

## Hypothesis Testing: Contingency table with chi-square test for political party and recession. 
## 
p <- sum(DJI$Recession)/length(DJI$Recession) # 17.67% of observations are in recession years
obs.tbl <- table(DJI$Republican,DJI$Recession)# Republican has more Recession
colnames(obs.tbl) <- c("Expansion", "Recession")
rownames(obs.tbl) <- c("Democrat", "Republican")
exp.tbl <- outer(rowSums(obs.tbl), colSums(obs.tbl))/sum(obs.tbl)
colnames(exp.tbl) <- c("Expansion", "Recession")
rownames(exp.tbl) <- c("Democrat", "Republican")
obs.tbl ; exp.tbl
chisq.test(DJI$Republican,DJI$Recession)
# p-value is less than 2.2e-16, far below our .05 threshold, so there would be a very
# small chance that the observed contingency table would arise by chance.
# Thus, the observations provide sufficient evidence to reject the null hypothesis
# that Republican and Democratic regimes are equally likely to be associated with recession
# years from 1985 to early 2020. 
#
#Running this as chi-square test of contingency table including all regimes:
obs.tbl <- table(DJI$Recession, DJI$Regime); rownames(obs.tbl) <- c("Expansion", "Recession"); obs.tbl #GWB had the most recession days
exp.tbl <- outer(rowSums(obs.tbl), colSums(obs.tbl))/sum(obs.tbl); rownames(exp.tbl) <- c("Expansion", "Recession"); exp.tbl
chisqvalue <- sum((obs.tbl - exp.tbl)^2/exp.tbl)
pchisq(chisqvalue, df = (2 - 1) * (6 - 1), lower.tail = FALSE) # p-value is 0
#We reject the null hypothesis that recession years arise by chance across regimes. 
#
#Running this as chi-square test specific to each regime with p the observed probabilty of recession:
q <- 1 - p; q # 0.8233235 probability of not being in a recession
prob <- (DJI$Recession*p + (!DJI$Recession)*q) / sum(DJI$Recession*p + (!DJI$Recession)*q) 
min(prob) ; max(prob) ; sum(prob) 
for (i in unique(DJI$Regime)) {
  print(chisq.test(DJI$Recession, DJI$Regime == i, p = prob))
}
# Null hypothesis is that each regime has the observed probability p of recession across regimes. 
# Note: There could exist carryover/lingering effects of recession or otherwise from one regime to the next
#Each p-value is less than 2.2e-16, far below the .05 threshold. This indicates
#that no individual regime is equally likely to be associated with recessions
#from the years 1985 to early 2020. (We should look into if our statistical methods
#are correct here.)

# Exploratory Data Analysis of Partial Variance to test for convergence of variance
N <- length(Open) ; 
variances.normal <- numeric(N - 1)
variances.cauchy <- numeric(N - 1)
variances.Open <- numeric(N - 1)
variances.diffs <- numeric(N - 1)
sample.normal <- rnorm(N) ; sample.cauchy <- rcauchy(N)
Open <- DJI$Open ; diffs.Open <- DJI$diffs
index <- 1:(N - 1)
for (i in 2:N) {
 variances.normal[i - 1] <- var(sample.normal[1:i])
 variances.cauchy[i - 1] <- var(sample.cauchy[1:i])
 variances.Open[i - 1] <- var(Open[1:i])
 variances.diffs[i - 1] <- var(diffs.Open[1:i])
}
variances.diffs <- variances.diffs[-1]
par(mfrow = c(2,2)) # create 2x2 plot matrix
plot(index,variances.normal, type = "l", col = "steelblue", log = "x", ylab = "Normal Variance", xlab = "Sample Size") # converges
plot(index,variances.cauchy, type = "l", col = "firebrick", log = "xy",ylab = "Cauchy Variance", xlab = "Sample Size") # diverges jagged
plot(index,variances.Open, type = "l", col = "yellowgreen", log = "xy",ylab = "Open Variance", xlab = "Sample Size") # diverges
plot(head(index,-1),variances.diffs, type = "l", col = "slategray", log = "xy", ylab = "Firs Diff. Variance", xlab = "Sample Size") # diverges
par(mfrow = c(1,1)) # revert to 1x1 plot matrix
summary(variances.normal) # data is centered closely around mean and median
summary(variances.cauchy) # seems to be large spread
summary(variances.Open) # extremely large spread
summary(variances.diffs) # spread is larger than it is for Cauchy but less than Open prices
# variance seems to diverge for both Open prices and for first differences

#Let us now try to fit our DJI data with some other distribution whose properties are familiar:
#Given the lack of convergence in the variance shown above, it appears that we should model the data
#with some distribution with infinite variance. Thus, we will see how well our data is modeled by a 
#Cauchy distribution. 
#
#install.packages("fitdistrplus")
library("fitdistrplus")
hist(diffs, prob = TRUE, breaks = "FD")
#Paramaters for Cauchy thanks to the paper by M. Mahdizadeh, and Ehsan Zamanzade.
#https://www.sciencedirect.com/science/article/pii/S1018364718313193?via%3Dihub
#Median:
diffs.median <- median(diffs); diffs.median #
#Half Interquartile Range:
diffs.hiq <-  (quantile(diffs)[[4]] - quantile(diffs)[[2]]) /2; diffs.hiq # 36.41016


#Checking our paramaters against the fitdist paramaters (nearly equal). But, because fit.diffs
#uses better paramater estimation, we will use our fit.diffs values.
fit.diffs <- as.vector(fitdist(diffs, "cauchy")$estimate); fit.diffs
curve(dcauchy(x, location = fit.diffs[1], scale = fit.diffs[2]), add = TRUE, lwd = 3, col = "red")
curve(dcauchy(x, location = diffs.median, scale = diffs.hiq), add = TRUE, lwd = 1, col = "cyan")
#And this does appear to be a good estimator

#Chi-Square Test for our Cauchy Distribution. We begin by setting up the breaks for the bins:
cauchy.breaks <- qcauchy((0:4)*.25, location = diffs.median, scale = diffs.hiq)
#Get observed data
cauchy.obs <- table(cut(diffs, breaks = cauchy.breaks)); cauchy.obs
#Get expected data:
cauchy.exp <- rep(length(diffs)/4, 4); cauchy.exp 

#Get initial Chi Square Statistic:
cauchy.cs <- ChiSq(cauchy.obs, cauchy.exp); cauchy.cs

#Run simulation: 
N <- 10^4; results.Cauchy <- numeric(N)
for (i in 1:N) {
  obs.sim.cauchy <- rcauchy(1:length(diffs), fit.diffs[1], fit.diffs[2])
  obs.sim.cauchy <- table(cut(obs.sim.cauchy, breaks = cauchy.breaks))
  results.Cauchy[i] <- ChiSq(obs.sim.cauchy, cauchy.exp)
}
hist(results.Cauchy, breaks = "FD")
abline(v = cauchy.cs, col = "red", lwd = 3)
cauchy.pvalue <- mean(results.Cauchy >= cauchy.cs); cauchy.pvalue # P value of approximately 24%.
#We fail to reject the null hypothesis. There is significant evidence to suggest that our data pulls
#from (at the very least) a family-member of heavy tail stable distributions. Although it is definitely
#NOT a gaussian distribution. It may or may not be Cauchy, but since all Stable distributions besides Gaussian
#have infinite variance, it appears our data's distribution also has infinite variance.

#################
#using Octiles
hist(diffs, prob = TRUE, breaks = "FD")
#Paramaters for Cauchy thanks to the paper by M. Mahdizadeh, and Ehsan Zamanzade.
#https://www.sciencedirect.com/science/article/pii/S1018364718313193?via%3Dihub
#Median:
diffs.median <- median(diffs); diffs.median #
#Half Interquartile Range:
diffs.hiq <- (quantile(diffs, c(seq(0.0,1,by = 0.125)))[[7]] - quantile(diffs, c(seq(0.0,1,by = 0.1)))[[3]])/2; diffs.hiq #44.32993

#Checking our paramaters against the fitdist paramaters (nearly equal). But, because fit.diffs
#uses better paramater estimation, we will use our fit.diffs values.
fit.diffs <- as.vector(fitdist(diffs, "cauchy")$estimate); fit.diffs
curve(dcauchy(x, location = fit.diffs[1], scale = fit.diffs[2]), add = TRUE, lwd = 3, col = "red")

#Chi-Square Test for our Cauchy Distribution. We begin by setting up the breaks for the bins:
cauchy.breaks <- qcauchy((0:8)*.125, location = diffs.median, scale = diffs.hiq)
#Get observed data
cauchy.obs <- table(cut(diffs, breaks = cauchy.breaks)); cauchy.obs
#Get expected data:
cauchy.exp <- rep(length(diffs)/8, 8); cauchy.exp 

#Get initial Chi Square Statistic:
cauchy.cs <- ChiSq(cauchy.obs, cauchy.exp); cauchy.cs

#Run simulation: 
N <- 10^4; results.Cauchy <- numeric(N)
for (i in 1:N) {
  obs.sim.cauchy <- rcauchy(1:length(diffs), fit.diffs[1], fit.diffs[2])
  obs.sim.cauchy <- table(cut(obs.sim.cauchy, breaks = cauchy.breaks))
  results.Cauchy[i] <- ChiSq(obs.sim.cauchy, cauchy.exp)
}
hist(results.Cauchy, breaks = "FD")
abline(v = cauchy.cs, col = "red", lwd = 3)
cauchy.pvalue <- mean(results.Cauchy >= cauchy.cs); cauchy.pvalue # P value of approximately 12%.
#We fail to reject the null hypothesis. There is significant evidence to suggest that our data pulls
#from (at the very least) a family-member of heavy tail stable distributions. Although it is definitely
#NOT a gaussian distribution. It may or may not be Cauchy, but since all Stable distributions besides Gaussian
#have infinite variance, it appears our data's distribution also has infinite variance.
#Pareto
library(MASS)
#install.packages("actuar")
library(actuar)

#The distribution of the absolute value of the changes appears to be well modeled by a Pareto distribution. Let us 
#investigate if that is the case:
hist(AbsDiffs , breaks = "FD", probability = TRUE) 
shape.pareto <- fitdist(AbsDiffs, "pareto", start=list(shape = 1, scale = 500))$estimate[1]
scale.pareto <- fitdist(AbsDiffs, "pareto", start=list(shape = 1, scale = 500))$estimate[2]
curve(dpareto(x,shape = shape.pareto , scale = scale.pareto ), col = "red", add = TRUE) #This looks pretty good

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
ChiStatAbs <- sum((binabschg - Exp)^2/Exp); ChiStatAbs #73.3071
#We estimated two parameters (using the sample mean and standard deviation), which costs two degrees of freedom, 
#and we have set the total days to match our sample which costs another so we have 10 - 3 = 7 degrees of freedom
curve(dchisq(x, df = 7), from = 0, to = ChiStatAbs + 5)
abline(v=ChiStatAbs, col = "red") 
pchisq(ChiStatAbs, df = 7, lower.tail = FALSE) # 0
#Given this extremely low p-value, it seems that the pareto distribution is not a good model (at all)
#for the absolute daily fluxes in the Dow Jones Industrial Average. 

#Let's shift now to taking a look at rare events, and how they are impacting some of our results:
#Let's consider each of the days in our data and examine how far from the mean flux in value each of
#them was. 
mu <- mean(diffs); mu
sigma <- sd(diffs); sigma

N <- length(diffs)
SDs <- numeric(N)
for (i in 1:N){
  SDs[i] <- (diffs[i]-mu)/sigma
}
head(SDs);head(diffs)
length(SDs)
SDs.data <- c(0,SDs[2:length(SDs)]); head(SDs.data) 
DJI <- data.frame(DJI,SDs.data)
idx <- which(abs(SDs.data) > 5); head(idx)
unusual <- DJI[idx,]; head(unusual) 
#As we can see, there are 32 days in which the price flux for the Dow Jones was larger than 5 standard 
#deviations away from the mean. To show just how bad of a fit the normal distribution is for our data
#consider the p-values for each of these events. 
N <- nrow(unusual)
pvals <- numeric(N)
for (i in 1:(N)) {
 pvals[i] <- pnorm(abs(unusual$SDs.data[i]*sigma), mean = mu, sd = sigma, lower.tail = FALSE)
}
head(pvals)
rare <- max(pvals);rare #2.303366e-07, which is pretty much 0. 

(1/rare)/365
#If we interpret the p-value as the probability of an event taking place, and our events are measured in days,
#then a given p-vakue tells us the probability of seeing that extreme of an event on any given days (i.e. a p-value of 
#0.05 would correspond to an event that we'd "expect" to see once every 20 days, since it has a 1/20 chance of arising).
#In other words, if the DJI first differences followed a normal distribution,we would expect to see the least rare 
#of these rare events once in 11894.45 years.

#Now that we have seen that the data follows a model with infinite variance, which somewhat resembles 
#that of a random walk, let us consider the hypothesis that political regimes have an impact on 
#the market, here clearly represented by the Dow Jones industrial Average. We can do this by considering
#the impact of political party on the performance of the market. 

regime <- lm(DJI$diffs~DJI$Republican + DJI$Recession)
#At first glance, it looks like republican regimes tend to be negatively correlated with growth in the Dow.
#But let us look a little closer:
summary(regime)
#Given the large standard error for the Republican coefficient, we cannot safely conclude that republican
#administrations correlate with losses in the Dow. (i.e. the p-value is .3, so we fail to reject the null hypothesis
#that this control variable has no impact on the response variable.)

#Let's see if we can get a statistically signifficant result for any of the individual presidents:
GHWB <- DJI$Regime == "GHWB"
BC <- DJI$Regime == "BC"
GWB <- DJI$Regime == "GWB"
BO <- DJI$Regime == "BO"
DJT <- DJI$Regime == "DJT"

pres.binary <- data.frame(GHWB, BC, GWB, BO, DJT)

regress.data <- data.frame(DJI,pres.binary)

#since we omitted the Ronald Reagan variable, each of the coefficients for the various presidents 
#represents the incremental surplus or deficit in the DJIA that occurred during each of the other 
#presidents' terms
ind.pres <- lm(regress.data$diffs ~ regress.data$Recession + regress.data$GHWB + regress.data$BC  + regress.data$GWB + regress.data$BO + regress.data$DJT)
ind.pres
#These results are interesting becuase they seem out of line with the recent (2018-2019) rhetoric of the Donald 
#Trump administration's success in boosting the DJIA to new heights.
summary(ind.pres)
#here again, we find that no presidents' presense in the White House had a significant impact on the Dow.

#So, let's exclude the days when 5 sigma + events took place, that way we'll keep only the data for "normal/typical"
#days.
regress.data.ne <- regress.data[-idx,]
ind.pres.2 <- lm(regress.data.ne$diffs ~ + regress.data.ne$Recession + regress.data.ne$GHWB + regress.data.ne$BC + regress.data.ne$GWB + regress.data.ne$BO + regress.data.ne$DJT)
summary(ind.pres.2)
#These results seem to indicate that, excluding days with extreme events, and controlling for recessions 
#(and nothing else), when we compare the performance of the DJIA during the previous 6 US presidents, 
#the results indeed appear to be most favorable to Donald Trump, and least favorable to George W Bush.
#It appears that the negative coefficient related to GWB comes from the effects of the 2008 housing crisis,
#which took place in the final months of his second term. Additionally, we should note that the only statistically 
#signifficant coefficient was that of Donald Trump, which had a p-value well below 0.01.

#Logistic Regression
#Given the earlier results, it seems that recessions play a large role on what happens wiht the Dow Jones Industrial Average.
#Let us inspect how economic recessions correlate with the performance of the DJIA. 
plot(DJI$diffs,DJI$Recession, xlim = c(-5000,5000))
MLL <- function(alpha, beta) {
  -sum( log( exp(alpha+beta*DJI$diffs)/(1+exp(alpha+beta*DJI$diffs)) )*DJI$Recession
        + log(1/(1+exp(alpha+beta*DJI$diffs)))*(1-DJI$Recession) )
}
#R has a function that will maximize this function of alpha and beta
#install.packages("stats4")   #needs to be run at most once
library(stats4)
results <- mle(MLL, start = list(alpha = 0, beta = 0)) #an initial guess is required
results@coef
curve( exp(results@coef[1] + results@coef[2]*x) / (1 + exp(results@coef[1] + results@coef[2]*x)),col = "blue", add=TRUE)
#This is a fairly interesting result because its graph looks different than the normal logistic curve.
#Of course, this becomes obvious when one considers the nature of the regression, namely that recessions
#are events that are expected to correlate with negative values of first differences (i.e. price drops).
#In any case, this provides some evidence to support the hypothesis that negative fluxes in the Dow Jones
#correlate with economic recessions.




