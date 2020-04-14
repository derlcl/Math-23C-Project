#Define the Chi-Squared function
ChiSq <-function(Obs,Exp){
  sum((Obs-Exp)^2/Exp)
}