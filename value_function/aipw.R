########################################################################
# Augmented Inverse Probability Weights (AIPW)
# Compute value using the IPW estimator
#
# inputs:
#   - data = data
#   - i = treatment policy that's being tested
#         272x1 vector of 4 levels
# outputs:
#   - value, a scalar
#
########################################################################

aipw = function(df, i)
{
  # ipw estimate
  ipw = ipw(df, i)
  ipw.value = ipw$value
  pi.d = ipw$pi.d
  c.d = ipw$c.d
  
  # or estimates
  q.or = glm(eventdeath ~ age + diam, data=df, family=binomial)
  
  # create new dataframe to be used for generating predictions based on new treatment decision i
  data.policy = df
  data.policy$treatment = i # overwrites the treatment variable with new treatment decision
  or.preds = predict(q.or, newdata=data.policy, type="response")
  
  # adjustment
  adjustment = mean((c.d - pi.d) / (pi.d) * or.preds)
  
  # compile estimate
  value = ipw.value - adjustment
  return(list(value=value))
}