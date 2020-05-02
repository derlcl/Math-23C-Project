# 10. Professional-looking software engineering (e.g defining and using your own functions).

#Define the Chi-Squared function
ChiSq <- function(Obs,Exp){
  sum((Obs - Exp)^2/Exp)
}

RunFourier <- function(ncoeff, data){
  myCos <- function(m) cos((1:ncoeff)*m*2*pi/ncoeff)
  mySin <- function(m) sin((1:ncoeff)*m*2*pi/ncoeff)
  
  #This function gives the Fourier coefficient a_m
  coeffA <- function(m){
    sum(data*2*myCos(m)/ncoeff)
  }
  #This function gives the Fourier coefficient b_m
  coeffB <- function(m){
    sum(data*2*mySin(m)/ncoeff)
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
