
flsa <- function(x, z, 
                 lambda1, lambda2, n, 
                 maxStep, tol, tau, flag, d, gamma1, gamma2) {
 # print(paste("Hey I'm here in flsa within flsa.R with d = ", d, sep=""))                   #### MODIFY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  infor <- matrix(as.double(rep(0, 4)), ncol=1)
  res <- .C('R_flsa', 
            x = as.matrix(as.double(rep(0, n)), ncol = 1),
            z = as.matrix(as.double(rep(0, (n-1) )), ncol = 1),
            infor=infor,
            as.double(lambda1), as.double(lambda2),
            as.integer(n), as.integer(maxStep),
            as.double(tol), as.integer(tau), 
            as.integer(flag),
            as.matrix(as.double(x), ncol = 1),
            as.matrix(as.double(z), ncol = 1),
            as.integer(d),
            as.integer(gamma1),
            as.integer(gamma2), 
            PACKAGE = "fusedlasso")
  res[1:3]
}
