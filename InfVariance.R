#Import libraries
library('pracma')
library('fitdistrplus')
library('MASS')
library('ggplot2')
library('gganimate')
library('gifski')
library('fractaldim')
library('fBasics')
library('stabledist')
library('car')

#Import Cleaned Data
source("prj_DataPreparation.R")
source("prj_Functions.R")

#This is a graphical representation of the Dow Jones Index from 1985 to 2020
plot(DJI$Open, type = "l", xlab = "Time", 
     ylab = "Random Daily Opens", main = "Random Walk Model",
     ylim = c(-10000,50000))

#Get the first difference of the time series data
diffs <- diff(DJI$Open)

## Fourier Analysis
 
#Can we capture DJI$Open by Fourier Analysis?
RunFourier(length(DJI$Open)/2, DJI$Open, FALSE) #Perfect Reconstruction 
RunFourier(10, DJI$Open, FALSE) #Using approx 1/400th of the basis vectors
#We capture a general trend of the market using very few basis vectors in our analysis.
#However, the information that has been captureed is relatively useless as there is no
#cyclical pattern to follow. 

#Can the first differences be captured by fourier analysis?
RunFourier(length(diffs) /2, diffs, TRUE) #Perfect Reconstruction
RunFourier(1000, diffs, TRUE) #Pretty good reconstruction using only a fraction of the basis vectors
RunFourier(10, diffs, TRUE) #Useless
#The first difference of the time series does even worse under the Fourier analysis.
#Using a small number of basis vectors yields almost no information. Although we can
#perfectly reconstruct the data, the more basis vectors added to our analysis mostly just
#pick up what I would consider "white noise". There is no noticeable trend or cycle.
#This is our first indication that the market is acting as a random walk.

## Random Walk

#The first step in testing to see if the Dow Jones Index has infinite variance
#is to see if the the index follows a Random Walk.

## We begin with a Chi Square test to test Random Walk model for daily flux in Open price values
#Get the mean of data
mu.chg.open <- mean(diffs); mu.chg.open # 2.142422 around what we expected based on the aforementioned values

#Get the standard deviation of our differences
sd.diff <- sqrt(mean(diffs^2) - mean(diffs)); sd.diff
#The sequence we will use to cut our data with.
rw.seq <- seq(from = 0, to = max(DJI$Open), by = max(DJI$Open) / 10); rw.seq
#Observed values
rw.obs <- as.vector(table(cut(DJI$Open, breaks = rw.seq))); rw.obs
#Set up Random Walk Expectation for our model
rw.exp <- rep(mean(rw.obs), 10); rw.exp
#ChiSq test for our Observed and Expected data
rw.cs <- ChiSq(rw.obs, rw.exp); rw.cs
#Set up Random Walk Simulation
N <- 10^4; results.RW <- numeric(N)
for (i in 1:N) {
  rw.drift <- arima.sim(model = list(order = c(0,1,0)), 
                        length(DJI$Open), mean = mu.chg.open,
                        sd = sd.diff)
  rw.sim <- as.vector(table(cut(rw.drift, breaks = rw.seq)))
  results.RW[i] <- ChiSq(rw.sim, rw.exp)
}
hist(results.RW)
abline(v = rw.cs, col = "red", lwd = 3)
rw.pvalue <- mean(rw.cs >= results.RW); rw.pvalue # .3841
#There is a 38.41% probability that would we encounter this extreme of a test statistic by random chance
#Our null hypothesis was that the daily flux in Open price values is well modeled by a random walk model. 
#We fail to reject this null hypothesis due to the p-value being greater than our .05 threshold for rejection.

# To demonstrate the idea of a random walk having infinite variance,
# we simulate random walk model graphically
plot(DJI$Open, type = "l", xlab = "Time", 
     ylab = "Random Daily Opens", main = "Random Walk Model",
     ylim = c(-10000,50000))

#Graphical volatility test (This is a phrase I made up)
for (i in 1:100) {
  rw.drift <- arima.sim(model = list(order = c(0,1,0)), 
                        length(DJI$Open), mean = mu.chg.open,
                        sd = sd.diff)
  lines(rw.drift, col = rgb(runif(1,0,1),runif(1,0,1),runif(1,0,1)))
}
#Extremely volatile -- the variance between simultation points at the same time is increasing as time
#increases.


## Fractal Analysis

#Using the hurst exponent, we can analyze whether our data is a random walk.

#We will be using the R/S Method to find the Hurst Exponent:
#Set the maximum number of divisons that will occur as we divide our data set in 1/2^n divisions
N <- floor(log(length(diffs), base = 2)) - 1

#Create a table to store the size of each divison (2^0, 2^1, 2^2,..., 2^n)
n <- numeric(N)

#Create table for our R/S result from each 
result.hexp <- numeric(N)

#Loop through each divison
for(i in 1:N){
  #Set division size
  n[i] <- floor(length(diffs)/2^(i-1))
  
  #"Chunk" our data into 2^n chunks
  ch <- split(diffs, cut_number(1:length(diffs), 2^(i - 1)), drop = TRUE)
  
  #Create table for R/S analysis average from each chunk for each division
  rs_avgs <- numeric(length(ch))
  
  #Loop through each chunk
  for (k in 1:length(ch)){
    
    #Set X = to the chunk
    X <- ch[[k]]
    
    #Get the mean of the chunk
    m <- mean(X)
    
    #Mean Adjusted Series:
    Y <- X - m; Y
    
    #Table for Cumulative Deviate Series:
    Z <- numeric(length(Y))
    
    #Calculate Cumulative Deviate Series
    for (j in 1:length(Z)){
      Z[j] <- sum(Y[1:j])
    }
    
    #Calculate range
    r <- max(Z) - min(Z); r
    
    #Calculate standard deviation of the the mean adjusted series:
    s <- sqrt(mean(Y^2))
    
    #Store R/S Analysis average
    rs_avgs[k] <- (r/s)
  }
  
  #Store the mean of the R/S analysis Average for each chunk
  result.hexp[i] <- mean(rs_avgs)
  
}

#Plot the log of the R/S analysis across the log of the size of the division
plot(log(result.hexp) ~ log(n))

#Plot Regression line
rl <- lm(log(result.hexp) ~ log(n)); rl
abline(rl$coefficients[1], rl$coefficients[2])

#Hurst Exponent:
hexp <- rl$coefficients[2]; hexp

#Sanity Test using the built in R Hurst Exponent:
hurstexp(diffs)

#Approximately .5 across all of the analysis. This indicates that our data is,
#in fact, a random walk.

#Since H = 2 - D where H is our Hurst Exponent and D is fractal dimension, our fractal dimension
#of our stock prices should be approximately 1.5
#Using the fractal dimensions package, we can see that our fractal dimension is approximately
#1.5 across 4 separate methods of fractal analysis:
#(For this analysis, we use the stock prices not the stationary first difference)

#install.packages("wavelets")
fracd.estimates <- fd.estimate(DJI$Open, methods = c("variogram", "madogram",
                                                     "hallwood", "wavelet")); fracd.estimates$fd
#One final check for our Hurst Exponent:
2 - fracd.estimates$fd

#Correct again. We get .5 across all 4 methods. The Dow Jones Index seems to be following a
#random walk.

#Standard deviation is undefined for models with infinite variance.
#Therefore, we will have to test our stock market drop against a gaussian random walk.
# Calculate the first difference series rw.drift_diff
rw.drift_diff <- diff(rw.drift)
rw.drift_diff
#Mean of the First Difference Series
rw.sigma <- sd(rw.drift_diff); rw.sigma #117.0062
rw.mu <- mean(rw.drift_diff); rw.mu #1.835033

#Plot a histogram of the values and draw the mean/median using abline
hist(rw.drift_diff); abline(v = mean(rw.drift_diff), col = "red", lwd = 3)
#Plot the values of the time series and draw the mean/median using abline - and we get a constant!
plot(rw.drift_diff, type = "l", xlab = "Time", ylab = "Differences in Random Daily Opens",
     main = "First Difference Series"); abline(h = mean(rw.drift_diff), col = "red", lwd = 3)
hist(rw.drift_diff, prob = TRUE)
curve(dnorm(x,rw.mu,rw.sigma), from = -500, to = 500, add = TRUE, col = "red") 

#Zooming in our plot:
ggplot(DJI[,c("Date", "Open")], aes(x=Date, y=Open, group=1)) + geom_line() +scale_x_discrete(breaks = 0)
ggplot(DJI[8000:8857,c("Date", "Open")], aes(x=Date, y=Open, group=1)) + geom_line() +scale_x_discrete(breaks = 0)
ggplot(DJI[8850:8857,c("Date", "Open")], aes(x=Date, y=Open, group=1)) + geom_line() +scale_x_discrete(breaks = 0)

#Approximately 5 Sigma Event for the average drop of the last 7 days in our data set:
abs(mean(diffs[8851:8857]) / sd.diff) 

#Approximately 20 sigma event for the biggest drop
abs(min(diffs) / sd.diff) 


## Cauchy Distribution Model
hist(diffs, prob = TRUE, breaks = "FD", main = "Histogram of First Differences
     Cauchy Model", xlab = "First Differences")
#Paramaters for Cauchy thanks to the paper by M. Mahdizadeh, and Ehsan Zamanzade.
#https://www.sciencedirect.com/science/article/pii/S1018364718313193?via%3Dihub
#Median:
diffs.median <- median(diffs); diffs.median #
#Half Interquartile Range:
diffs.hiq <-  (quantile(diffs)[[4]] - quantile(diffs)[[2]]) /2; diffs.hiq # 36.41016
# Result: Doesn't look normal, possibly Cauchy. 

#Checking our paramaters against the fitdist paramaters (nearly equal). But, because fit.diffs
#uses better paramater estimation, we will use our fit.diffs values.
fit.diffs <- as.vector(fitdist(diffs, "cauchy")$estimate); fit.diffs
curve(dcauchy(x, location = fit.diffs[1], scale = fit.diffs[2]), add = TRUE, lwd = 3, col = "red")

#Chi-Square Test for our Cauchy Distribution. We begin by setting up the breaks for the bins:
cauchy.breaks <- qcauchy((0:4)*.25, location = diffs.median, scale = diffs.hiq)
#Get observed data
cauchy.obs <- table(cut(diffs, breaks = cauchy.breaks)); cauchy.obs
#Get expected data:
cauchy.exp <- rep(length(diffs)/4, 4); cauchy.exp 
#Get initial Chi Square Statistic:
cauchy.cs <- ChiSq(cauchy.obs, cauchy.exp); cauchy.cs

#Run simulation: 
N <- 10^4; results.Cauchy <- pvalues.Chisq <- numeric(N)
for (i in 1:N) {
  obs.sim.cauchy <- rcauchy(1:length(diffs), fit.diffs[1], fit.diffs[2])
  obs.sim.cauchy <- table(cut(obs.sim.cauchy, breaks = cauchy.breaks))
  results.Cauchy[i] <- ChiSq(obs.sim.cauchy, cauchy.exp)
  pvalues.Chisq[i] <- pchisq(results.Cauchy[i],4-3,lower.tail=F)
}
hist(results.Cauchy, breaks = "FD", probability = TRUE, main = "Chi-Square Values Cauchy Simulation", xlab = "Chi-Square Stat")
abline(v = cauchy.cs, col = "red", lwd = 3)
cauchy.pvalue <- mean(results.Cauchy >= cauchy.cs); cauchy.pvalue# P value of approximately 24%.
#This would lead us fail to reject the null hypothesis. But the built in chi-square test gives us a different result:
pchisq(cauchy.cs, df = 4 - 3, lower.tail = FALSE)
#So we will try to use other goodness of fit tests that give more conclusive results.
#But in any way, there seems to be significant evidence to suggest that our data pulls from (at the very least) 
#a family-member of heavy tail stable distributions. Although it is definitely NOT a gaussian distribution.
#It may or may not be Cauchy, but since all Stable distributions besides Gaussian have infinite variance,
#it appears our data's distribution also has infinite variance.

#using Octiles
hist(diffs, prob = TRUE, breaks = "FD", main = "Histogram of First Differences
     Cauchy Model", xlab = "First Differences")
#Paramaters for Cauchy thanks to the paper by M. Mahdizadeh, and Ehsan Zamanzade.
#https://www.sciencedirect.com/science/article/pii/S1018364718313193?via%3Dihub
#Median:
diffs.median <- median(diffs); diffs.median #
#Half Interquartile Range:
diffs.hiq <- (quantile(diffs, c(seq(0.0,1,by = 0.125)))[[7]] - quantile(diffs, c(seq(0.0,1,by = 0.1)))[[3]])/2; diffs.hiq #44.32993
#Chi-Square Test for our Cauchy Distribution. We begin by setting up the breaks for the bins using estimated parameters:
cauchy.breaks <- qcauchy((0:8)*.125, location = diffs.median, scale = diffs.hiq)
cauchy.obs <- table(cut(diffs, breaks = cauchy.breaks)); cauchy.obs #Get observed data
cauchy.exp <- rep(length(diffs)/8, 8); cauchy.exp #Get expected data
cauchy.cs <- ChiSq(cauchy.obs, cauchy.exp); cauchy.cs #Get initial Chi Square Statistic
#Run simulation: 
N <- 10^4; results.Cauchy <- numeric(N)
for (i in 1:N) {
  obs.sim.cauchy <- rcauchy(1:length(diffs), fit.diffs[1], fit.diffs[2])
  obs.sim.cauchy <- table(cut(obs.sim.cauchy, breaks = cauchy.breaks))
  results.Cauchy[i] <- ChiSq(obs.sim.cauchy, cauchy.exp)
}
hist(results.Cauchy, breaks = "FD", main = "Chi-Square Values Cauchy Simulation", xlab = "Chi-Square Stat")
abline(v = cauchy.cs, col = "red", lwd = 3)
cauchy.pvalue <- mean(results.Cauchy >= cauchy.cs); cauchy.pvalue # P value of approximately 12%.
#This again would lead u to failt to reject the null hypothesis. But the built in chi-square test gives us a different result:
pchisq(cauchy.cs, df = 8 - 3, lower.tail = FALSE)
#So we will try to use other goodness of fit tests that give more conclusive results.
#But in any way, there seems to be significant evidence to suggest that our data pulls from (at the very least) 
#a family-member of heavy tail stable distributions. Although it is definitely NOT a gaussian distribution.
#It may or may not be Cauchy, but since all Stable distributions besides Gaussian have infinite variance,
#it appears our data's distribution also has infinite variance.

## Central Limit Theorem, Bootstrap and Empirical Cumulative Distribution
#
# Sources: 
# https://stats.stackexchange.com/questions/2504/test-for-finite-variance 
# Chihara and Hesterberg's Mathematical Statistics with Resampling
#
# If our sample is an independent and identically distributed realization 
# from a light-tailed distribution, the central limit theorem should hold
# even at smaller sample sizes. If we have a heavy-tailed distribution, 
# larger sample sizes will be needed to approximate the standard normal distribution.
# Use bootstrap resampling to demonstrate this.

## Bootstrap For A Single Population
par(mfrow = c(2,2)) # create 2x2 plot matrix
sampsize <- length(diffs) # starting sample size to draw from
N <- 100 # number of boostrap samples to run
n <-  100 # bootstrap sample size to draw
xlima <- -3 ; xlimb <- 3

# Perform N boostrap resamples of size n from sample X
meanY <- varY <- Z <- numeric(N)
plot(function(x) pnorm(x), xlim = c(xlima,xlimb), lwd = 5, main = "eCDF of Z from First Differences")
for (i in 1:N) {
  X <- diffs
  for (i in 1:N) {
    Y <- sample(X,n,replace = TRUE) # Resample
    meanY[i] <- mean(Y)
    varY[i] <- var(Y)
    Z[i] <- (mean(Y) - mean(X)) / (sd(Y)/sqrt(n))# Compute Z test statistic
  }
  lines(ecdf(Z), col = rgb(runif(1,0,1),runif(1,0,1),runif(1,0,1)), cex = .1)
}
meanY.diffs <- meanY
varY.diffs <- varY
Z.diffs <- Z
plot(function(x) pnorm(x), lwd = 3, add = TRUE)

# Perform N boostrap resamples of size n from sample Cauchy with interquartile parameters
meanY <- varY <- Z <- numeric(N)
plot(function(x) pnorm(x), xlim = c(xlima,xlimb), lwd = 5, main = "eCDF of Z from Cauchy (Interquartile)")
for (i in 1:N) {
  X <- rcauchy(sampsize, location = diffs.median, scale = diffs.hiq)
  for (i in 1:N) {
    Y <- sample(X,n,replace = TRUE) # Resample
    meanY[i] <- mean(Y)
    varY[i] <- var(Y)
    Z[i] <- (mean(Y) - mean(X)) / (sd(Y)/sqrt(n))# Compute Z test statistic
  }
  lines(ecdf(Z), col = rgb(runif(1,0,1),runif(1,0,1),runif(1,0,1)), cex = .1)
}
meanY.iqrCauchy <- meanY
varY.iqrCauchy <- varY
Z.iqrCauchy <- Z
plot(function(x) pnorm(x), lwd = 3, add = TRUE)

# Perform N boostrap resamples of size n from Cauchy with fitdist parameters
meanY <- varY <- Z <- numeric(N)
plot(function(x) pnorm(x), xlim = c(xlima,xlimb), lwd = 5, main = "eCDF of Z from Cauchy (FitDist)")
for (i in 1:N) {
  X <- rcauchy(sampsize, location = fit.diffs[1], scale = fit.diffs[2])
  for (i in 1:N) {
    Y <- sample(X,n,replace = TRUE) # Resample
    meanY[i] <- mean(Y)
    varY[i] <- var(Y)
    Z[i] <- (mean(Y) - mean(X)) / (sd(Y)/sqrt(n))# Compute Z test statistic
  }
  lines(ecdf(Z), col = rgb(runif(1,0,1),runif(1,0,1),runif(1,0,1)), cex = .1)
}
meanY.fdCauchy <- meanY
varY.fdCauchy <- varY
Z.fdCauchy <- Z
plot(function(x) pnorm(x), lwd = 3, add = TRUE)

# Perform N boostrap resamples of size n from normal distribution
meanY <- varY <- Z <- numeric(N)
plot(function(x) pnorm(x), xlim = c(xlima,xlimb), lwd = 5, main = "eCDF of Z from Normal")
for (i in 1:N) {
  X <- rnorm(sampsize, mean(diffs), sd(diffs)) 
  for (i in 1:N) {
    Y <- sample(X,n,replace = TRUE) # Resample
    meanY[i] <- mean(Y)
    varY[i] <- var(Y)
    Z[i] <- (mean(Y) - mean(X)) / (sd(Y)/sqrt(n))# Compute Z test statistic
  }
  lines(ecdf(Z), col = rgb(runif(1,0,1),runif(1,0,1),runif(1,0,1)), cex = .1)
}
meanY.norm <- meanY
varY.norm <- varY
Z.norm <- Z
plot(function(x) pnorm(x), lwd = 3, add = TRUE)

#########################################################
## Result: Empirical cumulative distribution for standardized random variable from
## bootstrap sampling distribution seems to indicate similarity between first differences
## and standard normal as well as dissimilarity between first differences and Cauchy. 
## Though there are indicates of wider dispersion in the eCDFs for the first differences
## than there are in the eCDFs for the normal, indicating the possibility of 
## significantly greater Kolmogorov-Smirnov test statistics. This requires testing. 


# Bootstrap sampling distribution of standardized sample observations
hist(Z.diffs, breaks = "fd", prob = TRUE) # looks normal
hist(Z.iqrCauchy, breaks = "fd", prob = TRUE) # somewhat normal, inconsistent center on iteration
hist(Z.fdCauchy, breaks = "fd", prob = TRUE) # somewhat normal, inconsistent center on iteration
hist(Z.norm, breaks = "fd", prob = TRUE) # looks normal
# Result: Central Limit Theorem explains approximation of standard normal for larger sample sizes. 
# However, due to infinite variance of Cauchy distribution, extreme values from heavy tails
# require very large sample size to approximate standard normal distribution. 
# Cauchy has inconsistent center, shape and spread on different iterations. 

# Bootstrap sampling distribution of sample variances
hist(varY.diffs, breaks = "fd", prob = TRUE) # somewhat normal with somewhat heavy right tail
hist(varY.iqrCauchy, breaks = "fd", prob = TRUE) # heavy right tail
hist(varY.fdCauchy, breaks = "fd", prob = TRUE) # heavy right tail
hist(varY.norm, breaks = "fd", prob = TRUE) # looks normal
# Result: Variance for first difference somewhat resembles shape and spread of Cauchy variance, not normal.

# Bootstrap sampling distribution of sample means
hist(meanY.diffs, breaks = "fd", prob = TRUE) # looks normal
hist(meanY.iqrCauchy, breaks = "fd", prob = TRUE) # looks like stable distribution with large variance
hist(meanY.fdCauchy, breaks = "fd", prob = TRUE) # looks like stable distribution with large variance
hist(meanY.norm, breaks = "fd", prob = TRUE) # looks normal
# Result: Mean for first difference resembles shape and spread of normal sample mean, not Cauhcy.

## Overall Result: Our first differences data may lie somewhere between Gaussian and Cauchy 
## on the parameter scale for stable distributions. However, the data is a sample from 
## an underlying population, so the limited accessibilility to a sample size of only 8857 
## limits the capture of extreme values from possibly heavy right tails. 

par(mfrow = c(1,1)) # reset to 1x1 plot matrix

# Partial Variance to test for convergence of variance
N <- length(Open) - 1; 
variances.normal <- variances.cauchy <- variances.Open <- variances.diffs <- numeric(N)
sample.normal <- rnorm(N + 1) ; sample.cauchy <- rcauchy(N + 1)
Open <- DJI$Open ; diffs.Open <- DJI$diffs
index <- 1:N
for (i in 2:(N + 1)) {
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
### Result: As index increases, partial variance converges for normal distribution,
### but it diverges in jagged jumps for Cauchy distribution and
### in smoother curves for both Open values and first differences. 
### This indicates that our data may have undefined or infinite variance. 




## Stable Distribution Model

#Paramaterized using:"Simple Consistent Estimators of Stable Distribution Parameters," Communications in Statistics - Simulation and Computation 15 (1986): 1109-36.
#https://www.asc.ohio-state.edu/mcculloch.2/papers/stabparm.pdf by J. Huston McCulloch
stable.Xs <- quantile(diffs, c(.05, .25, .5, .75, .95))

#Calculate V's
stable.V_a <- (stable.Xs[[5]] - stable.Xs[[1]]) / (stable.Xs[[4]] - stable.Xs[[2]]); stable.V_a
stable.V_b <- (stable.Xs[[5]] + stable.Xs[[3]] - (2*stable.Xs[[3]])) / (stable.Xs[[5]] - stable.Xs[[1]]); stable.V_b

#Using Table we calculate alpha and beta
stable.a <- 1.13
stable.b <- 0

#Calculate Phi_3 
stable.phi_3 <- 2.312

#Use phi_3 to calculate scale and then location is found from the table 
stable.c <- (stable.Xs[[4]] - stable.Xs[[2]]) / stable.phi_3; stable.c
stable.location <- median(diffs)

#Plot:
hist(diffs, breaks = "FD", freq = FALSE)
curve(dstable(x, (stable.a),(stable.b),stable.c, stable.location), add = TRUE, col = "cyan")

#Breaks:
stable.breaks <- c(-Inf, qstable((1:9) * .1, stable.a,stable.b,stable.c,stable.location), Inf)
chisq.test(stable.obs, table(cut(rstable(1:length(diffs), (stable.a),(stable.b - .5),stable.c, stable.location), stable.breaks)))

#Both Values
stable.obs <- table(cut(diffs, stable.breaks)); stable.obs
stable.exp <- rep(mean(stable.obs), length(stable.obs)); stable.exp

#Stable
stable.cs <- ChiSq(stable.obs, stable.exp)

#Simulation:
N <- 10^4; results.Stable <- numeric(N)
for (i in 1:N) {
  obs.sim.stable <- rstable(1:length(diffs), stable.a,stable.b,stable.c,stable.location)
  obs.sim.stable <- table(cut(obs.sim.stable, breaks = stable.breaks))
  results.Stable[i] <- ChiSq(obs.sim.stable, stable.exp)
}

mean(results.Stable >= stable.cs) #We reject the null hypothesis, our data fails to pull from the 
#Stable distribution. Our data is therefore, between a Cauchy Distribution and a Gaussian Distribution
#leaning toward a Cauchy distribution.

## Monte Carlo Simulation to Estimate Parameters for Stable Distribution
N <- 10^3 ; cs <- Inf ; idx <- 1 # set upper limit for chisquare and starting index
tempdata <- diffs[(diffs>-750)&(diffs<750)] # select data for chisquare test
alpha <- runif(N,0,2) # alpha parameter has support (0,2]
beta <- runif(N,-1,1) # beta parameter has support [-1,1]
gamma <- rexp(N,1) # gamma parameter has support (0,Inf)
delta <- median(tempdata) # delta parameter has support , estimated from data
for (i in 1:N) {
  stable.breaks <- c(-Inf, qstable( ((1:9) * .1), alpha[i],beta[i],gamma[i],delta), Inf)
  stable.obs <- table(cut(tempdata,stable.breaks)) ; 
  stable.exp <- rep(mean(stable.obs), length(stable.obs)); stable.exp   
  csCurrent <- ChiSq(stable.obs,stable.exp)
  if (csCurrent < cs) { 
    cs <- csCurrent # set lower chisquare value
    idx <- i # record the index
  }
}
cs ; pchisq(cs,df=5, lower.tail=FALSE) # high chisquare value, 0 p-value, poor fit from simulation
hist(tempdata, breaks = "fd", prob=TRUE) # histogram of data to model
curve(dstable(x,alpha[idx],beta[idx],gamma[idx],delta), add=TRUE, col = "red") # overlay density curve
alpha[idx] ; beta[idx] ; gamma[idx] ; delta # stable distribution parameters


#Kolmogorov-Smirnov Test
nn <- 500 #sample size. 
rand <- rcauchy(nn, location = diffs.median, scale = diffs.hiq); head(rand)
hist(rand, probability = TRUE, breaks = "FD")
curve(dcauchy(x, location = diffs.median, scale = diffs.hiq), add = TRUE, lwd = 3, col = "blue")

#using interquartile parameters
N <- 10^4
ks.stats <- numeric(N)
for (i in 1:N) {
  rand <- rcauchy(nn, location = diffs.median, scale = diffs.hiq); head(rand)
  diff.samp <- sample(diffs,nn, replace = TRUE)
  ks.stats[i] <- ks.test(diff.samp,rand)$p.value
}
mean(ks.stats)#mean p-values = 0.1976958
sum(ks.stats > 0.05)/length(ks.stats)#About 72% of the time our random samples generate a p-value greater than 0.05
{hist(ks.stats, main = "Histogram of K-S p-values", xlab = "K-S p-values")
  abline(v = 0.05, col = "red")}

#Using fitdist parameters
rand2 <- rcauchy(nn, location = fit.diffs[1], scale = fit.diffs[2]); head(rand)
hist(rand2, probability = TRUE, breaks = "FD")
curve(dcauchy(x, location = fit.diffs[1], scale = fit.diffs[2]), add = TRUE, lwd = 3, col = "blue") 

N <- 10^3
ks.stats2 <- numeric(N)
for (i in 1:N) {
  rand2 <- rcauchy(nn, location = fit.diffs[1], scale = fit.diffs[2]); head(rand)
  diff.samp2 <- sample(diffs,nn, replace = TRUE)
  ks.stats2[i] <- ks.test(diff.samp2,rand2)$p.value
}
mean(ks.stats2) #mean pvalues =  0.4248467
sum(ks.stats2 > 0.05)/length(ks.stats2)#About 92% of the time our random samples generate a p-value greater than 0.05
{hist(ks.stats2, breaks = "fd", main = "Histogram of K-S p-values", xlab = "K-S p-values")
  abline(v = 0.05, col = "red")}

#Taking N samples from diffs and using fitdist to fit each of those. Then comparing each of those samples to an 
#rCauchy dist with those same parameters, and then taking chi-squares and looking at the histogram
N <- 10^3
samp.rand <- numeric(N)
nn <- 500
ks.samp.rand <- numeric(N)
for(i in 1:N){
  diff.samp <- sample(diffs, nn, replace = TRUE)
  fit.samp <- fitdist(diff.samp, "cauchy", "mle")
  cauchy.samp <- rcauchy(nn, location = fit.samp$estimate[1], scale = fit.samp$estimate[2])
  ks.samp.rand[i] <- ks.test(diff.samp, cauchy.samp)$p.value
}
mean(ks.samp.rand)
sum(ks.samp.rand > 0.05)/length(ks.samp.rand)#About 99% of the time our random samples generate a p-value greater than 0.05
{hist(ks.samp.rand, breaks = "fd", main = "Histogram of K-S P-values", xlab = "K-S p-values")
  abline(v = 0.05, col = "red")}

#doing the same thing as above, but instead of Kolmogorov-Smirnov, doing chi-square test
samp.rand.cs <- numeric(N); pvals.chi <- numeric(N)
for(i in 1:N){
  diff.samp <- sample(diffs, nn, replace = TRUE)
  fit.samp <- fitdist(diff.samp, "cauchy", "mle")
  cauchy.samp <- rcauchy(nn, location = fit.samp$estimate[1], scale = fit.samp$estimate[2])
  cauchy.breaks <- qcauchy((0:10)*.1, location = fit.samp$estimate[1], scale = fit.samp$estimate[2])
  cauchy.obs <- table(cut(diff.samp, breaks = cauchy.breaks)); cauchy.obs
  cauchy.exp <- rep(nn/10, 10); cauchy.exp 
  samp.rand.cs[i] <- ChiSq(cauchy.obs, cauchy.exp); cauchy.cs
  pvals.chi[i] <-pchisq(samp.rand.cs[i], df = 7, lower.tail = FALSE)
}
hist(samp.rand.cs, breaks = "FD", main = "Histogram of Chi-Squared Stats", xlab = "K-S Stats")
{hist(pvals.chi, breaks = "FD",  main = "Histogram of Chi-Squared p-values", xlab = "p-values")
  abline(v = 0.05, col = "red")}
sum(pvals.chi > 0.05)/length(pvals.chi)# p-values are greater than 0.05 about 62% of the time

#Checking how the K-S statistic is affected by sample size
N <- 6000
ks.stat <- numeric(N)
for (i in 1:N) {
  diff.samp <- sample(diffs, i+1, replace = TRUE)
  fit.samp <- fitdist(diff.samp, "cauchy", "mle")
  cauchy.samp <- rcauchy(i + 1, location = fit.diffs[1], scale = fit.diffs[2])
  ks.stat[i] <- ks.test(diff.samp, cauchy.samp)$p.value
}
plot(ks.stat, pch = ".", main = "Plot of K-S p-values as Sample Size Increases", xlab = "Sample Size")
hist(ks.stat, breaks = "FD", probability = TRUE, main = "Histogram of K-S p-values", xlab = "K-S p-values" )

#Combining a few ideas:
plot(DJI$Open[seq(from = 1, to = length(diffs), by = 3)], typ = "l")
plot(DJI$Open[seq(from = 1, to = length(diffs), by = 7)], typ = "l")
plot(DJI$Open[seq(from = 1, to = length(diffs), by = 365)], typ = "l")
plot(DJI$Open[seq(from = 1, to = length(diffs), by = 365 * 3)], typ = "l")

#Hurstexponent
result.hexptrunc <- numeric(length(DJI$Open) / 5)
for (i in 1:(length(DJI$Open) / 5)){
  print(i)
  hexp <- hurstexp(diff(DJI$Open[seq(from = 1, to = length(diffs), by = i)]))
  result.hexptrunc[i] <- hexp$Hs
}
plot(result.hexptrunc, type = "l")

lmhurst <- lm(result.hexptrunc~seq(from = 1, to = length(result.hexptrunc), by = 1))
lmhurst <- as.vector(lmhurst$coefficients)
curve(x*lmhurst[2] +lmhurst[1], add = TRUE, col = "red", lwd = 3)

#KS.Test - Gaussian vs Cauchy at various "details"
nn <- (length(DJI$Open) / 2) - 1
result.kstest.gauss <- numeric(nn)
result.kstest.cauchy <- numeric(nn)
for (i in 1:(nn)) {
  samp <- diff(DJI$Open[seq(from = 1, to = length(diffs), by = i)]); (samp <- samp - mean(samp) / sd(samp))
  samp.median <- median(samp)
  samp.hiq <- (quantile(samp)[[4]] - quantile(samp)[[2]]) / 2
  samp.mean <- mean(samp)
  samp.sd <- sd(samp)
  kstest.cauchy <- ks.test(samp, "pcauchy", samp.median, samp.hiq)$p.value
  kstest.gauss <- ks.test(samp, "pnorm", samp.mean, samp.sd)$p.value
  result.kstest.cauchy[length(result.kstest.cauchy) - i] <- kstest.cauchy
  result.kstest.gauss[length(result.kstest.gauss) - i] <- kstest.gauss
}
plot(result.kstest.gauss, type = "l"); abline(h = 0.05, col = "red")
mean(result.kstest.gauss >= .05)
plot(result.kstest.cauchy, type = "l"); abline(h = 0.05, col = "red")
mean(result.kstest.cauchy >= .05)

#QQ Plot Test:
library(car)
par(mfrow = c(1,3))
qqPlot(diffs, "cauchy", id = FALSE);mtext("QQ-Plot of Cauchy, Normal, and Stable Distributions", side = 3, line = -2.5, outer = TRUE)
qqPlot(diffs, "norm", id = FALSE)
qqPlot(diffs, "stable", alpha = stable.a, beta = stable.b)
par(mfrow = c(1,1)) # reset plot partitions
#As we can see more quantiles from our data fit the Cauchy and Stables line and
#fall within the confidence interval when compared to a normal distribution.


bins <- 6

## Logarithm of Absolute First Differences of Open Values
logDiffs <- log(abs(diffs[diffs != 0]))
# Fit normal, stable, and Cauchy distributions to data and run chi-square goodness of fit test.

#Normal
hist(logDiffs, breaks = "FD", freq = FALSE)
curve(dnorm(x, mean(logDiffs), sd(logDiffs)), col = "red", lwd = 2, add = TRUE)
# Chi-square test
norm.breaks <- qnorm((0:bins) * 1/bins, mean(logDiffs), sd(logDiffs))
norm.obs <- table(cut(logDiffs, breaks = norm.breaks)); norm.obs
norm.exp <- rep(length(logDiffs) / bins, length(norm.obs)); norm.exp
norm.cs <- ChiSq(norm.obs, norm.exp); norm.cs
pchisq(norm.cs, df = bins - 3, lower.tail = F) # 0, Reject null

#Stable
stable.Xs <- quantile(logDiffs, c(.05, .25, .5, .75, .95))
#Calculate V's
stable.V_a <- (stable.Xs[[5]] - stable.Xs[[1]]) / (stable.Xs[[4]] - stable.Xs[[2]]); stable.V_a
stable.V_b <- (stable.Xs[[5]] + stable.Xs[[1]] - (2*stable.Xs[[3]])) / (stable.Xs[[5]] - stable.Xs[[1]]); stable.V_b
#Using Table we calculate alpha and beta
stable.a <- 2
stable.b <- -1
#Calculate Phi_3 
stable.phi_3 <- 1.908
#Use phi_3 to calculate scale and then location is found from the table 
stable.c <- (stable.Xs[[4]] - stable.Xs[[2]]) / stable.phi_3; stable.c
stable.location <- median(logDiffs)
curve(dstable(x, stable.a, stable.b, stable.c, stable.location), add = TRUE, lwd = 2, col = "blue")
# Chi-square test
stable.breaks <- c(-Inf, qstable((1:(bins - 1)) * 1/bins, alpha = stable.a, beta = stable.b, gamma = stable.c, delta = stable.location), Inf)
stable.obs <- table(cut(logDiffs, breaks = stable.breaks)); stable.obs
stable.exp <- rep(length(logDiffs) / bins, length(stable.obs)); stable.exp
stable.cs <- ChiSq(stable.obs, stable.exp); stable.cs # 9.198371
pchisq(stable.cs, df = bins - 3, lower.tail = F) # 0.002422305, reject null

#Cauchy
cauchy.median <- median(logDiffs)
cauchy.hiq <- (quantile(logDiffs)[[4]] - quantile(logDiffs)[[2]]) / 2
curve(dcauchy(x, cauchy.median, cauchy.hiq), add = TRUE, col = "green", lwd = 3)
# Chi-square test
cauchy.breaks <- qcauchy((0:bins) * 1/bins, cauchy.median, cauchy.hiq)
cauchy.obs <- table(cut(logDiffs, breaks = cauchy.breaks)); cauchy.obs
cauchy.exp <- rep(length(logDiffs) / bins, length(cauchy.obs))
cauchy.cs <- ChiSq(cauchy.obs, cauchy.exp); cauchy.cs # 9.198371
pchisq(cauchy.cs, df = bins - 3, lower.tail = F) # 0.002422305, reject null

#QQ Plot
par(mfrow = c(1,3))
qqPlot(logDiffs, "norm"); qqPlot(logDiffs, "cauchy"); qqPlot(logDiffs, "stable", alpha = stable.a, beta = stable.b, gamma = stable.c, delta = stable.location)
hurstexp(logDiffs) # Around 0.9, except for theoretical around 0.53
## Result: Normal, stable, and Cauchy distributions fail to pass the chi-square 
## goodness of fit test for the logarithm of absolute first differences of non-zero 
## open values, though Cauchy seems to perform best and to look the best on the qqPlot. 


## First Differences of Monthly Average of Open Values
monthly <- DJI; monthly$Date <- as.Date(monthly$Date)
monthly$Date <- format(as.Date(monthly$Date, format = "%d/%m/%Y"),"%Y-%m") 
monthly <-  aggregate(monthly[,2:4], list(monthly$Date), mean, drop = TRUE) 
colnames(monthly)[1] <- "Date"; head(monthly)
par(mfrow = c(1,2))
plot(monthly$Open, type = "l")
monthlyDiffs <- diff(monthly$Open); hist(diff(monthly$Open), freq = FALSE, breaks = "FD")
# Fit normal, stable, and Cauchy distributions to data and run chi-square goodness of fit test.

#Normal
curve(dnorm(x, mean(monthlyDiffs), sd(monthlyDiffs)), col = "red", lwd = 2, add = TRUE)
# Chi-square test
norm.breaks <- qnorm((0:bins) * 1/bins, mean(monthlyDiffs), sd(monthlyDiffs))
norm.obs <- table(cut(monthlyDiffs, breaks = norm.breaks)); norm.obs
norm.exp <- rep(length(monthlyDiffs) / bins, length(norm.obs)); norm.exp
norm.cs <- ChiSq(norm.obs, norm.exp); norm.cs # 79.25118
pchisq(norm.cs, df = bins - 3, lower.tail = F) # 0, Reject null

#Stable
stable.Xs <- quantile(monthlyDiffs, c(.05, .25, .5, .75, .95))
#Calculate V's
stable.V_a <- (stable.Xs[[5]] - stable.Xs[[1]]) / (stable.Xs[[4]] - stable.Xs[[2]]); stable.V_a
stable.V_b <- (stable.Xs[[5]] + stable.Xs[[1]] - (2*stable.Xs[[3]])) / (stable.Xs[[5]] - stable.Xs[[1]]); stable.V_b
#Using Table we calculate alpha and beta
stable.a <- 1.279
stable.b <- 0
#Calculate Phi_3 
stable.phi_3 <- 1.955
#Use phi_3 to calculate scale and then location is found from the table 
stable.c <- (stable.Xs[[4]] - stable.Xs[[2]]) / stable.phi_3; stable.c
stable.location <- median(monthlyDiffs); stable.location
curve(dstable(x, stable.a, stable.b, stable.c, stable.location), add = TRUE, lwd = 2, col = "blue")
# Chi-square test
stable.breaks <- c(-Inf, qstable((1:(bins-1)) * 1/bins, alpha = stable.a, beta = stable.b, gamma = stable.c, delta = stable.location), Inf)
stable.obs <- table(cut(monthlyDiffs, breaks = stable.breaks)); stable.obs
stable.exp <- rep(length(monthlyDiffs) / bins, length(stable.obs)); stable.exp
stable.cs <- ChiSq(stable.obs, stable.exp); stable.cs # 13.92417
pchisq(stable.cs, df = bins - 3, lower.tail = F) # 0.0009471195, reject null

#Cauchy
cauchy.median <- median(monthlyDiffs)
cauchy.hiq <- (quantile(monthlyDiffs)[[4]] - quantile(logDiffs)[[2]]) / 2
curve(dcauchy(x, cauchy.median, cauchy.hiq), add = TRUE, col = "green", lwd = 3)
# Chi-square test
cauchy.breaks <- qcauchy((0:bins) * 1/bins, cauchy.median, cauchy.hiq)
cauchy.obs <- table(cut(monthlyDiffs, breaks = cauchy.breaks)); cauchy.obs
cauchy.exp <- rep(length(monthlyDiffs) / bins, length(cauchy.obs))
cauchy.cs <- ChiSq(cauchy.obs, cauchy.exp); cauchy.cs # 10.85308
pchisq(cauchy.cs, df = bins - 3, lower.tail = F) # 0.0009863154, reject null

#QQ Plot
par(mfrow = c(1,3))
qqPlot(monthlyDiffs, "norm"); qqPlot(monthlyDiffs, "cauchy"); qqPlot(monthlyDiffs, "stable", alpha = stable.a, beta = stable.b, gamma = stable.c, delta = stable.location)
hurstexp(monthlyDiffs) # Around 0.6
## Result: Normal, stable, and Cauchy distributions fail to pass the chi-square 
## goodness of fit test for the first differences of monthly average of open values, 
## though the normal distribution seems to perform especially worse than the others.


## First Differences of Yearly Average of Open Values
yearly <- DJI; yearly$Date <- as.Date(yearly$Date)
yearly$Date <- format(as.Date(yearly$Date, format="%d/%m/%Y"),"%Y") 
yearly <-  aggregate(yearly[,2:4], list(yearly$Date), mean, drop = TRUE) 
colnames(yearly)[1] <- "Date"; head(yearly)
par(mfrow = c(1,2))
plot(yearly$Open, type = "l")
yearlyDiffs <- diff(yearly$Open); hist(diff(yearly$Open), freq = FALSE, breaks = "FD")

#Normal
curve(dnorm(x, mean(yearlyDiffs), sd(yearlyDiffs)), col = "red", lwd = 2, add = TRUE)
norm.breaks <- qnorm((0:bins) * 1/bins, mean(yearlyDiffs), sd(yearlyDiffs))
norm.obs <- table(cut(yearlyDiffs, breaks = norm.breaks)); norm.obs
norm.exp <- rep(length(yearlyDiffs) / bins, length(norm.obs)); norm.exp
norm.cs <- ChiSq(norm.obs, norm.exp); norm.cs # 
pchisq(norm.cs, df = bins - 3, lower.tail = F) # fail to reject null

#Stable
stable.Xs <- quantile(yearlyDiffs, c(.05, .25, .5, .75, .95))
#Calculate V's
stable.V_a <- (stable.Xs[[5]] - stable.Xs[[1]]) / (stable.Xs[[4]] - stable.Xs[[2]]); stable.V_a
stable.V_b <- (stable.Xs[[5]] + stable.Xs[[1]] - (2*stable.Xs[[3]])) / (stable.Xs[[5]] - stable.Xs[[1]]); stable.V_b
#Using Table we calculate alpha and beta
stable.a <- 1.388
stable.b <- -0.165
#Calculate Phi_3 
stable.phi_3 <- 1.795
#Use phi_3 to calculate scale and then location is found from the table 
stable.c <- (stable.Xs[[4]] - stable.Xs[[2]]) / stable.phi_3; stable.c
stable.location <- median(yearlyDiffs)
curve(dstable(x, stable.a, stable.b, stable.c, stable.location), add = TRUE, lwd = 2, col = "blue")
# Chi-square test
stable.breaks <- c(-Inf, qstable((1:(bins-1)) * 1/bins, alpha = stable.a, beta = stable.b, gamma = stable.c, delta = stable.location), Inf)
stable.obs <- table(cut(yearlyDiffs, breaks = stable.breaks)); stable.obs
stable.exp <- rep(length(yearlyDiffs) / bins, length(stable.obs)); stable.exp
stable.cs <- ChiSq(stable.obs, stable.exp); stable.cs # 4.885714
pchisq(stable.cs, df = bins - 3, lower.tail = F) # fail to reject null

#Cauchy
cauchy.median <- median(yearlyDiffs)
cauchy.hiq <- (quantile(yearlyDiffs)[[4]] - quantile(yearlyDiffs)[[2]]) / 2
curve(dcauchy(x, cauchy.median, cauchy.hiq), add = TRUE, col = "green", lwd = 3)
# Chi-square test
cauchy.breaks <- qcauchy((0:bins) * 1/bins, cauchy.median, cauchy.hiq)
cauchy.obs <- table(cut(yearlyDiffs, breaks = cauchy.breaks)); cauchy.obs
cauchy.exp <- rep(length(yearlyDiffs) / bins, length(cauchy.obs))
cauchy.cs <- ChiSq(cauchy.obs, cauchy.exp); cauchy.cs # 4.885714
pchisq(cauchy.cs, df = bins - 3, lower.tail = F) # reject null

#QQ Plot
par(mfrow = c(1,3))
qqPlot(yearlyDiffs, "norm"); qqPlot(yearlyDiffs, "cauchy"); qqPlot(yearlyDiffs, "stable", alpha = stable.a, beta = stable.b, gamma = stable.c, delta = stable.location)
hurstexp(yearlyDiffs) # Around 0.65
## Result: Normal and stable distributions pass while Cauchy fails to pass the chi-square
## goodness of fit test for the first differences of yearly average of open values, 
## though the normal distribution seems to perform better than the stable. 

## Comparing normal QQ Plots for yearly, monthly and daily averages of first difference
par(mfrow = c(1,2))
qqPlot(yearlyDiffs, "norm"); qqPlot(rnorm(length(yearlyDiffs), mean(yearlyDiffs), sd(yearlyDiffs)), "norm")
# Good fit on yearly
par(mfrow = c(1,2))
qqPlot(monthlyDiffs, "norm"); qqPlot(rnorm(length(monthlyDiffs), mean(monthlyDiffs), sd(monthlyDiffs)), "norm")
# OK fit on monthly
par(mfrow = c(1,2))
qqPlot(logDiffs, "norm"); qqPlot(rnorm(length(logDiffs), mean(logDiffs), sd(logDiffs)), "norm")
# Poor fit on daily
## Result: Data looks increasingly Gaussian as we average over longer time horizons. 


## First Differences of Logarithm of Open values with Chi-square tests
logDiffs <- diff(log(DJI$Open))

## Normal
hist(logDiffs, breaks = "FD", freq = FALSE, prob = TRUE)
curve(dnorm(x, mean(logDiffs), sd(logDiffs)), col = "red", lwd = 2, add = TRUE)
# Chi-square test
norm.breaks <- qnorm((0:bins) * 1/bins, mean(logDiffs), sd(logDiffs))
norm.obs <- table(cut(logDiffs, breaks = norm.breaks)); norm.obs
norm.exp <- rep(length(logDiffs) / bins, length(norm.obs)); norm.exp
norm.cs <- ChiSq(norm.obs, norm.exp); norm.cs # 779.7705
pchisq(norm.cs, df = bins - 3, lower.tail = F) # 0, Reject null

## Stable
stable.Xs <- quantile(logDiffs, c(.05, .25, .5, .75, .95))
#Calculate V's
stable.V_a <- (stable.Xs[[5]] - stable.Xs[[1]]) / (stable.Xs[[4]] - stable.Xs[[2]]); stable.V_a
stable.V_b <- (stable.Xs[[5]] + stable.Xs[[1]] - (2*stable.Xs[[3]])) / (stable.Xs[[5]] - stable.Xs[[1]]); stable.V_b
#Using Table we calculate alpha and beta
stable.a <- 1.484
stable.b <- 0
#Calculate Phi_3 
stable.phi_3 <- 1.939
#Use phi_3 to calculate scale and then location is found from the table 
stable.c <- (stable.Xs[[4]] - stable.Xs[[2]]) / stable.phi_3; stable.c
stable.location <- median(logDiffs)
curve(dstable(x, stable.a, stable.b, stable.c, stable.location), add = TRUE, lwd = 2, col = "blue")
# Chi-square test
stable.breaks <- c(-Inf, qstable((1:(bins-1)) * 1/bins, alpha = stable.a, beta = stable.b, gamma = stable.c, delta = stable.location), Inf)
stable.obs <- table(cut(logDiffs, breaks = stable.breaks)); stable.obs
stable.exp <- rep(length(logDiffs) / bins, length(stable.obs)); stable.exp
stable.cs <- ChiSq(stable.obs, stable.exp); stable.cs # 443.075
pchisq(stable.cs, df = bins - 3, lower.tail = F) # 0, Reject null

#Cauchy
cauchy.median <- median(logDiffs)
cauchy.hiq <- (quantile(logDiffs)[[4]] - quantile(logDiffs)[[2]]) / 2
curve(dcauchy(x, cauchy.median, cauchy.hiq), add = TRUE, col = "green", lwd = 3)
# Chi-square test
cauchy.breaks <- qcauchy((0:bins) * 1/bins, cauchy.median, cauchy.hiq)
cauchy.obs <- table(cut(logDiffs, breaks = cauchy.breaks)); cauchy.obs
cauchy.exp <- rep(length(logDiffs) / bins, length(cauchy.obs))
cauchy.cs <- ChiSq(cauchy.obs, cauchy.exp); cauchy.cs # 0.6046065
pchisq(cauchy.cs, df = bins - 3, lower.tail = F) # Reject null

#QQ Plot on daily first differences on logarithm of Open values
par(mfrow = c(1,3))
qqPlot(logDiffs, "norm"); qqPlot(logDiffs, "cauchy"); qqPlot(logDiffs, "stable", alpha = stable.a, beta = stable.b, gamma = stable.c, delta = stable.location)
# Cauchy looks like best fit. Most data is within the confidence interval. 
hurstexp(logDiffs) # Around 0.5
## Result: Normal and stable distributions fail to pass while Cauchy passes the chi-square
## goodness of fit test for the first differences of of the logarithm of Open values. 



## Divergent Integration For Calculating Variance

#Now let us try to show that, if our data is modeled by a Cauchy distribution with the location and 
#scale parameters above, that it does indeed have infinite variance. 
#To do that, we'll need to define functions that will describe the integrands that will allow us to
#use the tail-integral theorem. First, we'll define the integrand that will give us E(X).
#Since the underlying distribution can be modeled by a Cauchy distribution, we will use 
#dcauchy with the right parameters as our mu_x:
integrand <- function(x) dcauchy(x, location = fit.diffs[1], scale = fit.diffs[2])*x
#and now we have E(X):
exp.x <- integrate(f = integrand, lower = -Inf, upper = Inf)$value; exp.x 

#In the same manner, we can try to calculate E(X^2) so that we can get Var = E(X^2) - E(X)^2
integrand2 <- function(x) dcauchy(x, location = fit.diffs[1], scale = fit.diffs[2])*x^2
#And E(X^2)
exp.x2 <- integrate(f = integrand2, lower = -Inf, upper = Inf)$value; exp.x2 
#And it appears that the integral is divergent! This means that 
#Var = E(X^2) - E(X)^2 also diverges, and thus Var = Inf!