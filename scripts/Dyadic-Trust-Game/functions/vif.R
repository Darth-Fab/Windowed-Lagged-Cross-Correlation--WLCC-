# Variance Inflation Factor ------
#to check if there are variable that are covariating a lot 

vif.mer <- function (fit) { # info about the VIF: https://onlinecourses.science.psu.edu/stat501/node/347/ 
  ## adapted from rms::vif ## function: https://github.com/aufrank/R-hacks/blob/master/mer-utils.R 
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}
