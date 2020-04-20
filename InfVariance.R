#Import library
library('pracma')
library('MASS')
library('ggplot2')
library('gganimate')
library('gifski')
library('fractaldim')

#Import Cleaned Data
source("prj_DataPreparation.R")
source("prj_Functions.R")

#This is a graphical representation of the Dow Jones Index from 1985 to 2020
plot(DJI$Open, type = "l", xlab = "Time", 
     ylab = "Random Daily Opens", main = "Random Walk Model",
     ylim = c(-10000,50000))

#Get the first difference of the time series data
diffs <- diff(DJI$Open)

#Bruno: Demonstrate infinite variance by showing that the model is not well modeled by a Normal distribution.
#Then try an alternate model like Cauchy later. Check if daily price fluxes have infinite variance.
hist(diffs, breaks = "FD", probability = TRUE) #already doesn't look super promising
mu <- mean(diffs)
sigma <- sd(diffs) 
# do we multiply by n the number of sample observations to estimate the population variance?
# maybe not since this is the population? or is it?
curve(dnorm(x,mu,sigma), from = -2500, to = 1000, add = TRUE, col = "red")
n1 <- qnorm(0.1, mu, sigma); n1    #10% of the normal distribution lies below this value
pnorm(n1, mu, sigma)       
mean(diffs <= n1) #6.2%, not great, we would epxect something closer to 10% if it were Normal
#Now let's create a vector of deciles so that we can split our data and see if it falls as expected
dec <- qnorm(seq(0.0,1,by = 0.1), mu,sigma); dec  #11 bins
Exp <- rep(length(diffs)/10,10); Exp 
binchg <- numeric(10)
for (i in 1:10) {
  binchg[i] <- sum((diffs >= dec[i]) & (diffs <= dec[i + 1] )) ; binchg
}
#Finally we can test for uniformity using a chi-squared test.
ChiStat <- sum((binchg - Exp)^2/Exp); ChiStat #3581.397
#We estimated two parameters (using the sample mean and standard deviation), which costs two degrees of freedom, 
#and we have set the total days to match our sample which costs another so we have 10 - 3 = 7 degrees of freedom
curve(dchisq(x, df = 7), from = 0, to = ChiStat + 5)
abline(v = ChiStat, col = "red") 
pchisq(ChiStat, df = 7, lower.tail = FALSE) # 0
#Given this extremely low chi-square value, it seems that the normal distribution is not a good model (at all)
#for the daily fluxes in the Dow Jones Industrial Average. So let's now check how a model with infinite variance
#fits the data. 

#The first step in testing to see if the Dow Jones Index has infinite variance
#is to see if the the index follows a Random Walk.

## We begin with a Chi Square test to test Random Walk model for daily flux in Open price values
#Get the mean of data
mu.chg.open <- mean(diffs); mu.chg.open # 2.142422 around what we expected based on the aforementioned values

#Get the standard dviation of our differences
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

# To demonstrat the idea of a random walk having infinit variance,
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

#Fractal Analysis
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
fracd.estimates <- fd.estimate(DJI$Open, methods = c("variogram", "madogram",
                                                     "hallwood", "wavelet")); fracd.estimates$fd

#One final check for our Hurst Exponent:
2 - fracd.estimates$fd

#Correct again. We get .5 across all 4 methods. The Down Jones Index seems to be following a
#random walk.



#Checking a 5 sigma Event:
hist(diffs, prob = TRUE, breaks = "FD")
fit.diffs <- as.vector(fitdist(diffs, "cauchy")$estimate); fit.diffs
curve(dcauchy(x, location = fit.diffs[1], scale = fit.diffs[2]), add = TRUE, lwd = 3, col = "red")
abline(v = min(diffs), col = "green")

#Standard deviation and expectation are undefined for our model because it is a Cauchy distribution.
#Therefore, we will have to test
#our stock market drop against a gaussian random walk.

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