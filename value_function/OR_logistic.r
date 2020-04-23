########################################################################
# Outcome Regression (OR): Logistic 
# Computes the value using OR, usinglogistic regression to find the form of the Q function
#
# inputs:
#   - data = dataset
#   - i = treatment policy that's being tested
#         272x1 vector of 4 levels,which is then turned into:
# outputs:
#   - deltaOR_log, the resulting value (scalar)
#
########################################################################

library("dplyr")

OR_log = function(df, i){
  # Assign inputs to different var names
  newdat = df
  newdat = select(newdat, -c(chemo, amputation, hormonal)) # remove the treatment columns
  A = i
  
  # Create model for the betas, can be modified in the future.
  mod = glm(eventdeath ~ ., data = newdat, family = "binomial")
  
  # Instantiate new matrices for each treatment type
  X1 = newdat
  X1$treatment = as.factor(1)
  X2 = newdat
  X2$treatment = as.factor(2)
  X3 = newdat
  X3$treatment = as.factor(3)
  X4 = newdat
  X4$treatment = as.factor(4)
  
  # Will multiply the Q vector by 1 depending on the treatment policy, 0 if not
  vec1 = ifelse(A == 1, 1, 0)
  Q1 = predict(mod, X1) * vec1
  vec2 = ifelse(A == 2, 1, 0)
  Q2 = predict(mod, X2) * vec2 
  vec3 = ifelse(A == 3, 1, 0)
  Q3 = predict(mod, X3) * vec3
  vec4 = ifelse(A == 4, 1, 0)
  Q4 = predict(mod, X4) * vec4
  
  deltaOR_log = sum(Q1 + Q2 + Q3 + Q4) / N
  
  return(deltaOR_log)
}
