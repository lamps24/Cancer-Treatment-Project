########################################################################
# Inverse Probability Weights (IPW)
# Compute value using the IPW estimator
#
# inputs:
#   - df = data
#   - eta = vector of eta values that will form policy
#   - trt = column number of treatment being considered 
#         Probably:
#           5: chemo
#           6: hormonal
#           7: amputation
# outputs:
#   - value, a scalar
#   - pi.d, a 272x1 vector
#   - c.d, a 272x1 vector
#
########################################################################

library(nnet)
library(dplyr)

ipw = function(df, eta, trt)
{
  
  # convert eta values to policies - can adjust form of policy
  i = ifelse(df$age < eta[1] & df$diam > eta[2] & df$grade > eta[3], 1, 0) 
  
  # extract only the variables needed for fitting so we can tell it to use all vars  
  x = df[, 8:ncol(df)]
  y = df[, trt]
  df_for_fit = as.data.frame(cbind(y, x))
  
  # fit the selected model
  pi.h = glm(y ~ ., data=df_for_fit, family=binomial)
  
  # compute the fitted values (propensity scores)
  pscores = predict(pi.h, type="response")
  
  # compute pi.d, c.d
  pi.d = (pscores * i) + ((1 - pscores) * (1 - i))
  c.d = (df[, trt] * i) + ((1 - df[, trt]) * (1 - i))
  
  # treatment effect
  value = mean((df$eventdeath * c.d) / pi.d)
  return(list(value=value, pi.d=pi.d, c.d=c.d))  
}
