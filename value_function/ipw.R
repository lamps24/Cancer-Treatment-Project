########################################################################
# Inverse Probability Weights (IPW)
# Compute value using the IPW estimator
#
# inputs:
#   - data = data
#   - i = treatment policy that's being tested
#         272x1 vector of 4 levels,which is then turned into:
#         272x4 matrix of dummies indicating which treatment each patient will receive
# outputs:
#   - value, a scalar
#   - pi.d, a 272x1 vector
#   - c.d, a 272x1 vector
#
########################################################################

library(nnet)
library(fastDummies)

ipw = function(df, i)
{
  
  # extract only the variables needed for fitting so we can tell it to use all vars  
  df_for_fit = df[, 8:ncol(df)]
  
  # fit model
  pi.h = multinom(treatment ~ ., data=df_for_fit)
  
  # compute the fitted values (propensity scores) - each patient has 8 probs
  pscores = predict(pi.h, type="probs")
  
  # turn treatment and regimen factors into 8 dummies - needs to be numeric
  original_treatment = as.numeric(as.matrix(dummy_cols(df$treatment))[, 2:5])
  original_treatment = matrix(original_treatment, nrow=272, ncol=4)
  
  new_treatment = as.numeric(as.matrix(dummy_cols(i))[, 2:5])
  new_treatment = matrix(new_treatment, nrow=272, ncol=4)
  
  # compute pi.d - probability original treatment = new treatment
  # pscore corresponding to new treatment
  pi.d = rowSums(pscores * new_treatment)
  
  # compute c.d - indicator for whether treatment = regimen
  c.d = rowSums(original_treatment * new_treatment)
  
  # treatment effect
  value = mean((df$eventdeath * c.d) / pi.d)
  return(list(value=value, pi.d=pi.d, c.d=c.d))
}
