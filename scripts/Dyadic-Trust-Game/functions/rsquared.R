# Effect size f squared --------- 

fsquared.f <- function(model.incl, model.excl) { # Aiken and West (1991) as cited in Lorah (2018) - Large-Scale Assessments in Education 
  
  R.model.incl <- rsquared(model.incl)$Conditional
  R.model.excl <- rsquared(model.excl)$Conditional
  fsquared <- (R.model.incl - R.model.excl)/(1 - R.model.incl)
  fsquared 
} 

# small at a value of 0.02, medium at a value of 0.15, and large at a value of 0.35 (Cohen 1992).