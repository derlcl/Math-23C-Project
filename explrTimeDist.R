DJI <- read.csv('DJI.csv')
GSPC <- read.csv('GSPC.csv')
NASDAQ <- read.csv('IXIC.csv')

master <- data.frame(as.Date(DJI$Date), DJI$Open, GSPC$Open, NASDAQ$Open); colnames(master) <- c("Date", "DJI", "SAP500", "NASDAQ"); head(master)

{master_yearly <- master; master_yearly$Date <- format(as.Date(master$Date, format="%d/%m/%Y"),"%Y") 
  master_yearly <-  aggregate(master_yearly[,2:4], list(master_yearly$Date), mean, drop = TRUE) 
  colnames(master_yearly)[1] <- "Date"; head(master_yearly)
}

{master_monthly <- master; master_monthly$Date <- format(as.Date(master$Date, format="%d/%m/%Y"),"%Y-%m") 
  master_monthly <-  aggregate(master_monthly[,2:4], list(master_monthly$Date), mean, drop = TRUE) 
  colnames(master_monthly)[1] <- "Date"; head(master_monthly)
}


plot(master$DJI, type = "l", xlab = "Time", 
     ylab = "Random Daily Opens", main = "Random Walk Model")
plot(master_monthly$DJI, type = "l", xlab = "Time", 
     ylab = "Random Daily Opens", main = "Random Walk Model")
plot(master_yearly$DJI, type = "l", xlab = "Time", 
     ylab = "Random Daily Opens", main = "Random Walk Model")

plot(master$SAP500, type = "l", xlab = "Time", 
     ylab = "Random Daily Opens", main = "Random Walk Model")
plot(master_monthly$SAP500, type = "l", xlab = "Time", 
     ylab = "Random Daily Opens", main = "Random Walk Model")
plot(master_yearly$SAP500, type = "l", xlab = "Time", 
     ylab = "Random Daily Opens", main = "Random Walk Model")

plot(master$NASDAQ, type = "l", xlab = "Time", 
     ylab = "Random Daily Opens", main = "Random Walk Model")
plot(master_monthly$NASDAQ, type = "l", xlab = "Time", 
     ylab = "Random Daily Opens", main = "Random Walk Model")
plot(master_yearly$NASDAQ, type = "l", xlab = "Time", 
     ylab = "Random Daily Opens", main = "Random Walk Model")

diffs <- master$DJI
diffs.median <- median(diffs); diffs.median #
diffs.hiq <- (quantile(diffs, c(seq(0.0,1,by = 0.125)))[[7]] - quantile(diffs, c(seq(0.0,1,by = 0.1)))[[3]])/2; diffs.hiq
hurstexp(diff(master$DJI))
hurstexp(diff(master_monthly$DJI))
hurstexp(diff(master_yearly$DJI))

ks.test(diff(master$DJI), rnorm(length(diff(master$DJI)), mean(diff(master$DJI)), sd(diff(master$DJI))))
ks.test(diff(master_monthly$DJI), rnorm(length(diff(master_monthly$DJI)), mean(diff(master_monthly$DJI)), sd(diff(master_monthly$DJI))))
ks.test(diff(master_yearly$DJI), rnorm(length(diff(master_yearly$DJI)), mean(diff(master_yearly$DJI)), sd(diff(master_yearly$DJI))))

ks.test(diff(master$DJI), rcauchy(length(diff(master$DJI)), diffs.median, diffs.hiq))
ks.test(diff(master_monthly$DJI), rcauchy(length(diff(master_monthly$DJI)), diffs.median, diffs.hiq))
ks.test(diff(master_yearly$DJI), rcauchy(length(diff(master_yearly$DJI)), diffs.median, diffs.hiq))


#Stopped here

hurstexp(diff(master$SAP500))
hurstexp(diff(master_monthly$SAP500))
hurstexp(diff(master_yearly$SAP500))

ks.test(diff(master$SAP500), rnorm(length(diff(master$SAP500)), mean(diff(master$SAP500)), sd(diff(master$SAP500))))
ks.test(diff(master_monthly$SAP500), rnorm(length(diff(master_monthly$SAP500)), mean(diff(master_monthly$SAP500)), sd(diff(master_monthly$SAP500))))
ks.test(diff(master_yearly$SAP500), rnorm(length(diff(master_yearly$SAP500)), mean(diff(master_yearly$SAP500)), sd(diff(master_yearly$SAP500))))

hurstexp(diff(master$NASDAQ))
hurstexp(diff(master_monthly$NASDAQ))
hurstexp(diff(master_yearly$NASDAQ))

ks.test(diff(master$NASDAQ), rnorm(length(diff(master$NASDAQ)), mean(diff(master$NASDAQ)), sd(diff(master$NASDAQ))))
ks.test(diff(master_monthly$NASDAQ), rnorm(length(diff(master_monthly$NASDAQ)), mean(diff(master_monthly$NASDAQ)), sd(diff(master_monthly$NASDAQ))))
ks.test(diff(master_yearly$NASDAQ), rnorm(length(diff(master_yearly$NASDAQ)), mean(diff(master_yearly$NASDAQ)), sd(diff(master_yearly$NASDAQ))))


