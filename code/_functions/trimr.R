#---------------------------------------------------------------------------------------------------
#
#  Description: Returns matrix/vector stripped of specified rows.
#  
#  Inputs:
#      x : M by N input matrix
#      a : first a rows to strip
#      b : last b rows to strip
#    
#  Outputs:
#      z : (M-a) by (N-b) stripped matrix
#
#---------------------------------------------------------------------------------------------------

trimr <- function (x, a, b) {
    # Preliminaries
    x <- as.matrix(x)
    m <- nrow(x)
    n <- ncol(x)
    
    h1 <- a + 1
    h2 <- m - b
    z <- matrix(x[h1:h2, ], ncol = n, byrow = F)
    
    return(z)
}
