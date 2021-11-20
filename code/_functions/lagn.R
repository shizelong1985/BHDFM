#---------------------------------------------------------------------------------------------------
#  
#  Returns lagged matrix/vector based on n. The rows are shifted down by the value of n and zeros 
#  are filled in the blank rows.
#  
#  Inputs:
#  --------------------------------------------- 
#      x  : M by N input matrix
#      n  : number of lags
#    
#  Outputs:
#  --------------------------------------------- 
#      y : M by N lagged matrix
#      
#---------------------------------------------------------------------------------------------------

 
lagn <- function (x, n) {
    
    # Preliminaries
    x <- as.matrix(x)
    nr <- nrow(x)
    nc <- ncol(x)
    y <- matrix(0, nrow = nr, ncol = nc)
    
    if (n > 0) {
        z <- trimr(x, 0, n)
        y[(n+1):nr, ] <- z
        
    } else if (n < 0) {
        z <- trimr(x, abs(n), 0)
        y[1:(nr-abs(n)), ] <- z
    } else {
        y <- x
    }
    
    return(y)
}
