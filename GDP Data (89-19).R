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


