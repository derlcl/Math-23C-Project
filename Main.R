#Math 23C Term Project by Rakeen Tanveer, Bruno Kömel, and Derl Clausen

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
#
# Compare open price first difference and trade volume variables
x <- log(DJI$Volume)
y <- log(abs(DJI$chg))
plot(x~y)
linmod <- lm(y~x)
summary(linmod) # R-squared value is only .1854, low explanation for residual errors
#





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
#
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