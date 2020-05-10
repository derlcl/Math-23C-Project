#Math 23C Term Project by Rakeen Tanvir, Bruno Kömel, and Derl Clausen

#Comments
##Titles
###Results/Conclusions

#Packages to install:
#install.packages("MASS")
#install.packages("pracma")
#install.packages("ggplot2")
#install.packages("gganimate")
#install.packages("fractaldim")
#install.packages("fitdistrplus")
#install.packages("gifski")
#install.packages("fBasics")
#install.packages("stabledist")
#install.packages("car")

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

##Retrieve any functions we have made for this project
#Define the Chi-Squared function
ChiSq <- function(Obs,Exp){
  sum((Obs - Exp)^2/Exp)
}

RunFourier <- function(ncoeff, data){
  myCos <- function(m) cos((1:length(data))*m*2*pi/length(data))
  mySin <- function(m) sin((1:length(data))*m*2*pi/length(data))
  
  #This function gives the Fourier coefficient a_m
  coeffA <- function(m){
    sum(data*2*myCos(m)/length(data))
  }
  #This function gives the Fourier coefficient b_m
  coeffB <- function(m){
    sum(data*2*mySin(m)/length(data))
  }
  
  #Evaluate ncoeff coefficients
  FourierA <- sapply(1:ncoeff,coeffA)
  FourierB <- sapply(1:ncoeff,coeffB)
  #Find the amplitude of the sum of the cosine and sine terms
  Fourier <- sqrt(FourierA^2+FourierB^2)
  Fourier.max <- which.max(Fourier)   #the largest value is for m
  
  #Reconstruction
  recon <- mean(data)   #this is a_0
  for (m in 1:ncoeff) {
    recon <- recon + FourierA[m]*myCos(m)+FourierB[m]*mySin(m)
  }
  
  plot(data,type = "l")
  points(1:length(data),recon, type = "l", col = "red",lwd = 2) 
}


##Load Dow Jones Industrial dataset and prepare it for analysis
DJI <- read.csv("DJI.csv")

## Exploratory Data Analysis
summary(DJI)
# Minimum, 1st quartile, median, mean, 3rd quartile, and maximum values look similar
# across Open, High, Low, Close and Adjusted Close. We will work with Open values.
# Add first differences between daily Open values to data frame
Open <- DJI$Open # retrieve our Open price values
diffs <- diff(DJI$Open) #Get first difference of our Open values. 
diffs.length <- length(diffs) # 8857, concatenate with 0 as starting value to fit data frame
DJI <- data.frame(DJI,"diffs" = c(0,diffs)) ; head(DJI) # Add to data frame

## Center and Shape
summary(diffs)  ; boxplot(diffs) # Five-number summary and boxplot visualization
# Min.        1st Qu.    Median      Mean    3rd Qu.  Max. 
# -2419.920   -30.400     3.430     2.142    42.420   1171.961 
# Median is greater than mean, indicates left skew
# Minimum is absolutely greater than maximum, indicates left skew
mean(diffs, trim = 0.25) # 4.351473 trimmed mean, greater than mean, indicates left skew
max(diffs) - min(diffs) # 3591.881 is the range, suggesting wide spread
IQR(diffs) # 72.82031 is the interquartile range, suggesting concentrated center
var(diffs) # whopping 13760.49 for variance, we expect this to increase over time
sd(diffs) # 117.3051 standard deviation

# Visualize
hist(diffs, breaks = "fd", prob = TRUE, main = "Histogram of First Differences", xlab = "First Differences") 
abline(v = mean(diffs), col = "turquoise")
# resembles normal distribution with narrow concentration around sharp peak at mean 
# and heavy tails, with more extreme negative values than positive values
qqnorm(diffs) # plot data with quantiles of the standard normal on the x-axis.
qqline(diffs) # add straight line through the first and third quartiles of our data
# data does not seem to follow normal distribution at tails

# Magnitude of First Differences
AbsDiffs <- abs(diffs) # the absolute value of the first differences
mu.AbsDiffs <- mean(AbsDiffs) # 69 is the mean, greater than trimmed mean
mean(AbsDiffs, trim = .25) # 41.49568 is the trimmed mean, indicating right skew
sum(AbsDiffs/diffs.length) # 69 is also the total value of the contributions to the mean
sum(AbsDiffs[AbsDiffs <= mu.AbsDiffs]/diffs.length)/mu.AbsDiffs 
# only 23.9% of the contributions to the mean come from values at or below the mean
sum(AbsDiffs <= mu.AbsDiffs)/diffs.length
# 67.7% of the values are at or below the mean value
# yet 76.1% of the contributions to the mean value come from values above the mean
max.AbsDiffs <- max(AbsDiffs); max.AbsDiffs # 2419.92 is the max value
max(AbsDiffs)/diffs.length / mu.AbsDiffs # single maximum value contributed .4% to the mean value
# there are long tails of extreme values with large contributions to the mean
hist(AbsDiffs, breaks = "fd", prob = TRUE) ; abline(v = mu.AbsDiffs)
# could be modeled by a non-negative valued, long-tailed distribution 


## Logarithm of Absolute First Differences
logAbsDiffs <- log(AbsDiffs) ; summary(logAbsDiffs) ; 
logAbsDiffs <- logAbsDiffs[logAbsDiffs > -Inf] ; summary(logAbsDiffs) ; length(logAbsDiffs)
hist(logAbsDiffs, prob = TRUE, breaks = "fd") # could be modeled by a stable distribution with left skew

## Empirical Cumulative Distributions
plot.ecdf(diffs)
AbsDiffsCDF <- ecdf(AbsDiffs)
plot(AbsDiffsCDF) # could be modeled by non-negative, long-tailed distribution
#
LogAbsDiffs <- log(AbsDiffs)
plot(LogAbsDiffs) # could be modeled by random walk with positive drift value
plot.ecdf(LogAbsDiffs) # produces error because of presence of -Inf values
#
CDF.diffs <- ecdf(diffs)
plot(CDF.diffs) # logistic regression model could fit
plot(Open) # exponential model could fit
plot(log(Open)) # linear regression or polynomial model could fit

## Regression
# An interesting correlation that we would not have expected to be significant, 
# but turns out to be so:
# Finding interesting linear relationships
Volume <- DJI$Volume
plot(DJI$diffs~Volume)
plot(DJI$diffs~log(Volume))
plot(log(abs(DJI$diffs))~log(Volume))
plot(DJI$diffs~log(Open))
plot(Open,Volume)
plot(log(Open)~log(Volume)) # looks most promising
# Compare logarithms of open price and trade volume variables
x <- log(Volume)
y <- log(Open)
#Let us first manually run the regression, and then we'll check our results against the built-in
#function. First, let's find the slope of the line
b <- sum( (x-mean(x))*(y-mean(y)) / sum((x-mean(x))^2)) # 0.6014286
#Alternative - works because division by n-1 cancels out
cov(x,y)/var(x) # 0.6014286
#Here is the formula for the intercept
a <- mean(y) - b*mean(x) ; a # -2.041207

#It is quicker to use the built-in R function
linmod <- lm(y~x)
linmod; a; b #And we get the same result
plot(y~x, main = "Plot of First Differences by Volume", xlab = "log(Volume)", ylab = "log(Absolute First Differences)")
abline(a, b, col = "cyan")
abline(linmod$coefficients[1], linmod$coefficients[2], col = "red")
summary(linmod) 
###R-squared is 0.7298, so the linear model explains 73% of the data, and it appears that volume and open
###prices are positively correlated.

##Let's shift now to taking a look at presumably rare events, and how they are impacting our results:
#Let's consider each of the days in our data and examine how far from the mean flux in value 
#each of them was. 
mu <- mean(diffs); mu
sigma <- sd(diffs); sigma

N <- length(diffs)
SDs <- numeric(N)
for (i in 1:N){
  SDs[i] <- (diffs[i]-mu)/sigma
}
head(SDs);head(diffs)
length(SDs)
SDs.data <- c(0,SDs[1:length(SDs)]); head(SDs.data) 
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
###If we interpret the p-value as the probability of an event taking place, and our events are measured in days,
###then a given p-value tells us the probability of seeing that extreme of an event on any given days (i.e. a p-value of 
###0.05 would correspond to an event that we'd "expect" to see once every 20 days, since it has a 1/20 chance of arising).
###In other words, if the DJI first differences followed a normal distribution,we would expect to see the least rare 
###of these rare events once in 11,894.45 years, but from the data it is clear that these events are far more
###common than that.


##Prepare data for political analysis
#
#Create 'Regime' Categories (Regan, Bush, etc.)
DJI["Regime"] <- "None" #Set blank column
#Set the dates for each regime. Each 'term' is denoted from inaugaration day to the day 
#before the next regime's inaugaration day.
dates <- sapply(DJI["Date"], function(x) as.Date(x)) 
#Convert CSV date character objects to date objects.
DJI[dates < as.Date("1989-01-20"), "Regime"] <- "RR"
DJI[dates >= as.Date("1989-01-20") & dates < as.Date("1993-01-20"), "Regime"] <- "GHWB"
DJI[dates >= as.Date("1993-01-20") & dates < as.Date("2001-01-20"), "Regime"] <- "BC"
DJI[dates >= as.Date("2001-01-20") & dates < as.Date("2009-01-20"), "Regime"] <- "GWB"
DJI[dates >= as.Date("2009-01-20") & dates < as.Date("2017-01-20"), "Regime"] <- "BO"
DJI[dates >= as.Date("2017-01-20"), "Regime"] <- "DJT"
#
#Create a binary category for the political party in control of the White House 
#(Republican = 1, Democrat = 0)
Republican <- (DJI$Regime == "RR") | (DJI$Regime == "GHWB") | (DJI$Regime == "GWB") | (DJI$Regime == "DJT")
DJI <- data.frame(DJI, Republican); head(DJI); sum(Republican); mean(Republican)
#
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

#Note that we got this dataset from the FRED website, and this is real GDP growth (which had larger fluxes)
#Note 2: I did this based on quarter over quarter fluxes.
GDPb <- read.csv("GDPC1.csv"); head(GDPb) 
N <- length(GDPb$GDPC1)
GDPg <- numeric(N) #now creating a vector to store the flux in GDP for each quarter
for(i in 2:N){ 
  #since we're not interested the data until 1989, we can start from the second row
  GDPg[i] <- GDPb$GDPC1[i] - GDPb$GDPC1[i-1]
}
Recess <- numeric(N)
for(i in 1:N){
  Recess[i] <- GDPg[i] < 0
}
GDPb <- data.frame(GDPb,GDPg, Recess) ; GDPb$DATE <- as.Date(GDPb$DATE); head(GDPb)
#selecting only the relevant years for our analysis
GDP <- subset(GDPb, DATE >= "1989-01-01" & DATE <= "2019-01-01"); head(GDP) 
sum(GDP$Recess) #So there have been 12 quarters with recessions
GDP[GDP$Recess == 1,]

##Surveying the impact of the White House on the Dow Jones
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
###These results are interesting becuase they seem out of line with the recent (2018-2019) rhetoric of the Donald 
###Trump administration's success in boosting the DJIA to new heights.
summary(ind.pres)
###here again, we find that no presidents' presence in the White House had a significant impact on the Dow.

#Let us now exclude the days when 5 sigma + events took place, that way we'll keep only the data for 
#"normal/typical" days.
regress.data.ne <- regress.data[-idx,]
ind.pres.2 <- lm(regress.data.ne$diffs ~ + regress.data.ne$Recession + regress.data.ne$GHWB + regress.data.ne$BC + regress.data.ne$GWB + regress.data.ne$BO + regress.data.ne$DJT)
summary(ind.pres.2)
###These results seem to indicate that, excluding days with "extremely" events, and controlling for recessions 
###(and nothing else), when we compare the performance of the DJIA during the previous 6 US presidents, 
###the results indeed appear to be most favorable to Donald Trump, and least favorable to George W Bush.
###It appears that the negative coefficient related to GWB comes from the effects of the 2008 housing crisis,
###which took place in the final months of his second term. Additionally, we should note that the only statistically 
###signifficant coefficient was that of Donald Trump, which had a p-value well below 0.01.


##Logistic Regression: Recessions and Results
#Given the earlier results, it seems that recessions have a large impact on 
#first differences in daily Open values for the Dow Jones Industrial Average.
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
###This is a fairly interesting result because its graph looks different than the normal logistic curve.
###Of course, this becomes obvious when one considers the nature of the regression, namely that recessions
###are events that are expected to correlate with negative values of first differences (i.e. price drops).
###In any case, this provides some evidence to support the hypothesis that negative fluxes in the Dow Jones
###correlate with economic recessions.


## Hypothesis Testing With Permutation Test

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

##Combined Permutation Test
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
### 2.87% chance that this extreme of an observed difference would arise by chance, so it appears that the 
###DJI performed better during democratic regimes, a result that is statistically signifficant.

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
#As we can see from this contingency table, Republicans had more days in office during recessions, but 
#they also had more days in office during expansions.
chisq.test(DJI$Republican,DJI$Recession)
### p-value is less than 2.2e-16, far below our .05 threshold, so there would be a very
### small chance that the observed contingency table would arise by chance.
### Thus, the observations provide sufficient evidence to reject the null hypothesis
### that Republican and Democratic regimes are equally likely to be associated with recession
### years from 1985 to early 2020. 
#
#Running this as chi-square test of contingency table including all regimes:
obs.tbl <- table(DJI$Recession, DJI$Regime); rownames(obs.tbl) <- c("Expansion", "Recession"); obs.tbl #GWB had the most recession days
exp.tbl <- outer(rowSums(obs.tbl), colSums(obs.tbl))/sum(obs.tbl); rownames(exp.tbl) <- c("Expansion", "Recession"); exp.tbl
#This table allows us to see a breakdown of how long each president was in office in terms of recessions and 
#expansions
chisqvalue <- sum((obs.tbl - exp.tbl)^2/exp.tbl)
pchisq(chisqvalue, df = (2 - 1) * (6 - 1), lower.tail = FALSE) # p-value is 0
#We reject the null hypothesis that recession years are equaly likely to arise across regimes. 
#
#Running this as chi-square test specific to each regime with p the observed probabilty of recession:
q <- 1 - p; q # 0.8233235 probability of not being in a recession
prob <- (DJI$Recession*p + (!DJI$Recession)*q) / sum(DJI$Recession*p + (!DJI$Recession)*q) 
min(prob) ; max(prob) ; sum(prob) 
for (i in unique(DJI$Regime)) {
  print(chisq.test(DJI$Recession, DJI$Regime == i, p = prob))
}
###Null hypothesis is that each regime has the observed probability p of recession across regimes. 
###Note: There could exist carryover/lingering effects of recession or otherwise from one regime to the next
###Each p-value is less than 2.2e-16, far below the .05 threshold. This indicates
###that no individual regime is equally likely to be associated with recessions
###from the years 1985 to early 2020. (We should look into if our statistical methods
###are correct here.)


# Fractal Tribute to Mandelbrot
# "This is just a tribute!
# You gotta believe it!
# And I wish you were there!
# Just a matter of opinion."
# - Tenacious D

# Source: https://rosettacode.org/wiki/Mandelbrot_set#R 
iterate.until.escape <- function(z, c, trans, cond, max=50, response=dwell) {
  #we iterate all active points in the same array operation,
  #and keeping track of which points are still iterating.
  active <- seq_along(z)
  dwell <- z
  dwell[] <- 0
  for (i in 1:max) {
    z[active] <- trans(z[active], c[active]);
    survived <- cond(z[active])
    dwell[active[!survived]] <- i
    active <- active[survived]
    if (length(active) == 0) break
  }
  eval(substitute(response))
}

re = seq(-2, 1, len=500)
im = seq(-1.5, 1.5, len=500)
c <- outer(re, im, function(x,y) complex(real=x, imaginary=y))
x <- iterate.until.escape(array(0, dim(c)), c,
                          function(z,c)z^2+c, function(z)abs(z) <= 2,
                          max=100)
image(x)
