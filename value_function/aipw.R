########################################################################
# Augmented Inverse Probability Weights (AIPW)
# Compute value using the IPW estimator
#
# inputs:
#   - df = data
#   - trt = column number of treatment being considered
#         Probably:
#           5 = chemo
#           6 = hormonal
#           7 = amputation
#   - i = treatment policy (automatically created in the policy search function)
#
# outputs:
#   - value, a scalar
#   - ipw.value = the ipw value (as a QC)
#   - adjustment = the adjustment that was made to the ipw
#
########################################################################

aipw = function(df, trt, i)
{
  
  # ipw estimate
  ipw = ipw(df, trt, i)
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
