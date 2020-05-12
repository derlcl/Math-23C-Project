# 10. Professional-looking software engineering (e.g defining and using your own functions).

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
