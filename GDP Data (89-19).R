#Note that we got this dataset from the FRED website, and this is real GDP growth (which had larger fluxes)
#Note 2: I did this based on year over year fluxes, but the .csv file has data by quarter
GDPb <- read.csv("GDPC1.csv"); head(GDP) 

BOY <- subset(GDPb, format.Date(GDP$DATE, "%m")=="01"); head(BOY) #selecting only dates for beginning of year
BOY$DATE <- as.Date(BOY$DATE); BOY <- data.frame(BOY) #formating dates and the data
N <- length(BOY$GDPC1)
GDPg <- numeric(N) #now creating a vector to store the flux in GDP for each year
for(i in 2:N){ 
  #since we're not interested the data until 1989, we can start from the second row
  GDPg[i] <- BOY$GDPC1[i] - BOY$GDPC1[i-1]
}
BOY <- data.frame(BOY,GDPg) 

GOM <- numeric(length(BOY$GDPg)) #GOA here stands for Growth Over Mean
for(i in 1:length(BOY$GDPg)){
  GOM[i] <- BOY$GDPg[i] > mean(GDPg) #using mean GDP growth here because 
  #there have been very few YOY recessions, so this will tell us
  #whether the growth in a given year was better than the average year from 1947-2019)
}
GDPf <- data.frame(BOY,GOM); head(GDPf)

#selecting only the relevant years for our analysis
GDP <- subset(GDPf, DATE >= "1989-01-01" & DATE <= "2019-01-01"); head(GDP) 


