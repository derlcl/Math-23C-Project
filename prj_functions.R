# 10. Professional-looking software engineering (e.g defining and using your own functions).

#Define the Chi-Squared function
ChiSq <- function(Obs,Exp){
  sum((Obs - Exp)^2/Exp)
}

