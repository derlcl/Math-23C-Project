DJI <- read.csv('DJI.csv')
GSPC <- read.csv('GSPC.csv')
NASDAQ <- read.csv('IXIC.csv')

master <- data.frame(as.Date(DJI$Date), DJI$Open, GSPC$Open, NASDAQ$Open); colnames(master) <- c("Date", "DJI", "S&P500", "NASDAQ"); head(master)

{master_yearly <- master; master_yearly$Date <- format(as.Date(master$Date, format="%d/%m/%Y"),"%Y") 
  master_yearly <-  aggregate(master_yearly[,2:4], list(master_yearly$Date), mean, drop = TRUE) 
  colnames(master_yearly)[1] <- "Date"; head(master_yearly)
}

{master_monthly <- master; master_monthly$Date <- format(as.Date(master$Date, format="%d/%m/%Y"),"%Y-%m") 
  master_monthly <-  aggregate(master_monthly[,2:4], list(master_monthly$Date), mean, drop = TRUE) 
  colnames(master_monthly)[1] <- "Date"; head(master_monthly)
}
