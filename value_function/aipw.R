########################################################################
# Augmented Inverse Probability Weights (AIPW)
# Compute value using the IPW estimator
#
# inputs:
#   - data = data
#   - eta = vector of eta values that will form policy
#   - trt = column number of treatment being considered 
#         Probably:
#           5: chemo
#           6: hormonal
#           7: amputation
# outputs:
#   - value, a scalar
#
########################################################################

aipw = function(df, eta, trt)
{
  
  # convert eta values to policies - can adjust form of policy
  i = ifelse(df$age < eta[1] & df$diam > eta[2] & df$grade > eta[3], 1, 0) 
  
  # ipw estimate
  ipw = ipw(df, eta, trt)
  ipw.value = ipw$value
  pi.d = ipw$pi.d
  c.d = ipw$c.d
  
  # or estimates - fitting with all vars selected by LASSO
  df_for_fit = df[, c(2, 8:ncol(df))]
  q.or = glm(eventdeath ~ ., data=df_for_fit, family=binomial)
  
  # predictions with new policy (replace original treatment with new)
  data.policy = df
  data.policy[, trt] = i
  or.preds = predict(q.or, newdata=data.policy, type="response")
  
  # adjustment
  adjustment = mean((c.d - pi.d) / (pi.d) * or.preds)
  
  # compile estimate
  value = ipw.value - adjustment
  return(list(value=value, ipw.value=ipw.value, adjustment=adjustment))
  
}
