########################################################################
# Augmented Inverse Probability Weights (AIPW)
# Compute value using the IPW estimator
#
# inputs:
#   - data = data
#   - eta = vector of eta values that will form policy
# outputs:
#   - value, a scalar
#
########################################################################

aipw = function(df, eta)
{
  
  # convert eta values to policies - can adjust form of policy
  chem = ifelse(df$var1 < eta[1], 1, 0) 
  amp = ifelse(df$timerecurrence < eta[2], 1, 0)
  policy = rep(0,272)
  policy = ifelse(chem == 0 & amp == 0, 1, policy)
  policy = ifelse(chem == 0 & amp == 1, 2, policy)
  policy = ifelse(chem == 1 & amp == 0, 3, policy)
  policy = ifelse(chem == 1 & amp == 1, 4, policy)
  i = policy
  
  # ipw estimate
  ipw = ipw(df, i)
  ipw.value = ipw$value
  pi.d = ipw$pi.d
  c.d = ipw$c.d
  
  # or estimates - fitting with all vars selected by LASSO
  df_for_fit = df[, c(2, 8:ncol(df))]
  q.or = glm(eventdeath ~ ., data=df_for_fit, family=binomial)
  
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
