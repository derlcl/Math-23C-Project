#This is where our main code for our project will go!
DJI <- read.csv("DJI.csv"); head(DJI)

#Create 'Regime' Categories (Regan, Bush, etc.)
DJI["Regime"] <- "None" #Set blank column

#Set the dates for each regime
#Each 'term' is denoted from inaugaration day 
#to the day before the next regime's inaugaration day
dates <- sapply(DJI["Date"], function(x) as.Date(x)) #Convert CSV date character objects to date objects
DJI[dates < as.Date("1989-01-20"), "Regime"] <- "RR"
DJI[dates >= as.Date("1989-01-20") & dates < as.Date("1993-01-20"), "Regime"] <- "GHWB"
DJI[dates >= as.Date("1993-01-20") & dates < as.Date("2001-01-20"), "Regime"] <- "BC"
DJI[dates >= as.Date("2001-01-20") & dates < as.Date("2009-01-20"), "Regime"] <- "GWB"
DJI[dates >= as.Date("2009-01-20") & dates < as.Date("2017-01-20"), "Regime"] <- "BO"
DJI[dates >= as.Date("2017-01-20"), "Regime"] <- "DJT"

#Create a binary category for the political party in control of the White House (Republican = 1, Democrat = 0)
Republican <- (DJI$Regime == "RR")|(DJI$Regime == "GHWB")|(DJI$Regime=="GWB")|(DJI$Regime=="DJT")
DJI <- data.frame(DJI, Republican); head(DJI); sum(Republican); mean(Republican)

