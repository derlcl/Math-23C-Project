#Math 23C Term Project by Rakeen Tanvir, Bruno Kömel, and Derl Clausen

#For the grader's sake, we have indexed our problem by number of #'s in each comment in the following way:

#Comments 
##Titles
###Results/Conclusions
####{Rubric Cross Reference}

#Additionally, for your convenience, we pasted the grading rubric here, and cross referenced the sections in our script
#that we believe meet the respective requirements:

## Required dataset standards (DONE)
#1. A dataframe. (DJI)
#2. At least two categorical or logical columns. (DJI)
#3. At least two numeric columns. (DJI)
#4. At least 20 rows, preferably more, but real-world data may be limited. (DJI - 8858 rows)

## Required graphical displays (all graphs must be colored and nicely labeled) (TD)
#1. A barplot. (line #422)
#2. A histogram. (line #158 is the first)
#3. A probability density graph overlaid on a histogram. (lines #733 and #734 is one example)
#4. A contingency table. (line #480 is one example)

## Required analysis (TD)
#1. A permutation test. (lines #396-470 we run a series of permutations)
#2. A p-value or other statistic based on a distribution function. (line #268 is one example)
#3. Analysis of a contingency table. (line #482 is one example)
#4. Comparison of analysis by classical methods (chi-square, CLT) and simulation methods. (lines #1035-#1045)

## Required submission uploads (DONE)
#1. A .csv file with the dataset (DJI)
#2. A long, well-commented script that loads the dataset, explores it, and does all the analysis. (Main.R file)
#3. A shorter .Rmd with compiled .pdf or .html file that presents highlights in ten minutes. (Report.Rmd file)
#4. A one-page handout that explains the dataset and summarizes the analysis. (https://www.overleaf.com/project/5e7665484b0d3600011e16d0)


## Additional points for creativity or complexity (ONLY KEPT THE ONES THAT ARE RELEVANT)
# 1. [DJI Accumulated Columns over Main.R] A data set with lots of columns, allowing comparison of many different variables. (by volume, decade, year, month) (Recession by Party and Regime)
#2.[DJI] A data set that is so large that it can be used as a population from which samples are taken. (35 x ~365 rows) (Bootstrap sampling distribution)
#4.[Overleaf] A one-page document that discusses ethical issues raised by conclusions reached from analysis of the data. (What would infinite variance in economic price indices mean for the long-term future of humanity?)
#5.[ECDF] A graphical display that is different from those in the textbook or in the class scripts.
#6.[Line #764] Appropriate use of R functions for a probability distribution other than binomial, normal, or chi-square.
#7.[Line #1415] Appropriate use of integration to calculate a significant result.
#8.[Line #468] A convincing demonstration of a relationship that might not have been statistically significant but that turns out to be so.
#9.[Line #792]A convincing demonstration of a relationship that might have been statistically significant but that turns out not to be so.
#10.[Yes] Professional-looking software engineering (e.g defining and using your own functions).
#11.[Line #738] Nicely labeled graphics using ggplot, with good use of color, line styles, etc., that tell a convincing story.
#12.[Lines #396-470] An example where permutation tests or other computational techniques clearly work better than classical methods. 
#13.[Hurst exponent, stable dist.] Appropriate use of novel statistics (e.g. trimmed mean, maximum or minimum, skewness, ratios). 
#14.[Line #329] Use of linear regression. (keep aberrant data but disregard events?, boolean variable, e.g. coronavirus, black monday)
#15.[Line #385] Calculation and display of a logistic regression curve.
#16.[Line #219] Appropriate use of covariance or correlation.
#18.[Bootstrapping] Use of theoretical knowledge of sampling distributions.
#19.[ECDF] A graphical display that is different from those in the class scripts.
#21.[Line #768] Appropriate use of quantiles to compare distributions. 


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

RunFourier <- function(ncoeff, data, change){
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
  
  if (change == FALSE){
    plot(data,type = "l",xlab = "Trading Day", ylab = "Stock",
         main = "Stock Price over Trading Days")
    points(1:length(data),recon, type = "l", col = "red",
           lwd = 2) 
  }
  else{
  plot(data,type = "l",xlab = "Trading Day", ylab = "Stock Price Change",
       main = "Stock Price Change over Trading Days")
  points(1:length(data),recon, type = "l", col = "red",
         lwd = 2) 
  }
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
####{Histogram Requirement}
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
####{Add'l Point #16}
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
####{P-value requirement}


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
####{Add'l Points #14}
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
####{Add'l Points #15}
###This is a fairly interesting result because its graph looks different than the normal logistic curve.
###Of course, this becomes obvious when one considers the nature of the regression, namely that recessions
###are events that are expected to correlate with negative values of first differences (i.e. price drops).
###In any case, this provides some evidence to support the hypothesis that negative fluxes in the Dow Jones
###correlate with economic recessions.


## Hypothesis Testing With Permutation Test 
####{Permutation Test Requirement}
####{Add'l Points #12}

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
####{Barplot requirement}

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
####{Add'l Points # 8}

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
####{Contingency Table Requirement}
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
####{Overlaying Curve Requirement}

#Zooming in our plot:
ggplot(DJI[,c("Date", "Open")], aes(x=Date, y=Open, group=1)) + geom_line() +scale_x_discrete(breaks = 0)
ggplot(DJI[8000:8857,c("Date", "Open")], aes(x=Date, y=Open, group=1)) + geom_line() +scale_x_discrete(breaks = 0)
ggplot(DJI[8850:8857,c("Date", "Open")], aes(x=Date, y=Open, group=1)) + geom_line() +scale_x_discrete(breaks = 0)
####{Add'l Points #11}

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
####{Add'l Points #6}

#Chi-Square Test for our Cauchy Distribution. We begin by setting up the breaks for the bins:
cauchy.breaks <- qcauchy((0:4)*.25, location = diffs.median, scale = diffs.hiq)
####{Add'l Points #21}
#Get observed data
cauchy.obs <- table(cut(diffs, breaks = cauchy.breaks)); cauchy.obs
#Get expected data:
cauchy.exp <- rep(length(diffs)/4, 4); cauchy.exp 
#Get initial Chi Square Statistic:
cauchy.cs <- ChiSq(cauchy.obs, cauchy.exp); cauchy.cs
pchisq(cauchy.cs, df = 4 - 3, lower.tail = FALSE) # 0, reject null

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
pchisq(cauchy.cs, df = 8-3, lower.tail=FALSE) # 0 , reject null
####{Add'l Points #9}

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
####{Non-Classical Requirement}

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
pchisq(norm.cs, df = bins - 3, lower.tail = F) #reject null

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
####{Add'l Points #7}
