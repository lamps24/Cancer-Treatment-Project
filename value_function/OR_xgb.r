########################################################################
# Outcome Regression (OR): XGBoost
# Computes the value using OR, using XGBoost to find the form of the Q function
#
# inputs:
#   - data = dataset
#   - eta = values of eta for the policy for chemo and amputation
#      - currently has length 2, but can be changed depending on how we structure policy
# outputs:
#   - deltaOR_log, the resulting value (scalar)
#
########################################################################

library("dplyr")
library("xgboost")

OR_xgb = function(df, eta){
  xdat = df
  N = length(xdat$age)
  
  ##########################################################
  # Convert eta values to policies
  # Can be modified to accept several values for eta
  # Right now there is only 2, one for predicting chemo and one for amputation
  
  chem = ifelse(xdat$age < eta[1], 1, 0) 
  amp = ifelse(xdat$timerecurrence < eta[2], 1, 0)
  
  policy = rep(0,272)
  
  policy = ifelse(chem == 0 & amp == 0, 1, policy)
  policy = ifelse(chem == 0 & amp == 1, 2, policy)
  policy = ifelse(chem == 1 & amp == 0, 3, policy)
  policy = ifelse(chem == 1 & amp == 1, 4, policy)
  A = policy
  
  #########################################################
  
  # separating dataset into covariates and outcome (label)
  xdat = select(xdat, -c(chemo, amputation, hormonal))
  xdat$treatment = as.numeric(xdat$treatment)
  xdat = as.matrix(xdat)
  labels = xdat[, 2]
  xdat = xdat[, -2]
  
  # split the data into training and testing set
  index = sample(seq(1,272), 100)
  train = xdat[-index, ]
  test = xdat[index, ]
  train_labels = labels[-index]
  test_labels = labels[index]
  
  # convert to XGBoost appropriate matrices
  dtrain = xgb.DMatrix(train, label = train_labels)
  dtest = xgb.DMatrix(test, label = test_labels)
  watchlist = list(train = dtrain, eval = dtest)
  
  # create list of parameters, and begin training the model
  param = list(max_depth = 2, eta = 1, verbose = 0, nthread = 2, objective = "binary:logistic", eval_metric = "auc")
  bst = xgb.train(param, dtrain, nrounds = 10, watchlist)
  
  # remove the treatment columns from input data
  df = select(df, -c(chemo, amputation, hormonal))
  
  X1 = df[, -2]
  X1$treatment = 1
  X1 = as.matrix(X1)
  X1 = xgb.DMatrix(X1)
  
  X2 = df[, -2]
  X2$treatment = 2
  X2 = as.matrix(X2)
  X2 = xgb.DMatrix(X2)
  
  X3 = df[, -2]
  X3$treatment = 3
  X3 = as.matrix(X3)
  X3 = xgb.DMatrix(X3)
  
  X4 = df[, -2]
  X4$treatment = 4
  X4 = as.matrix(X4)
  X4 = xgb.DMatrix(X4)
  
  vec1 = ifelse(A == 1, 1, 0)
  Q1 = predict(bst, X1) 
  vec2 = ifelse(A == 2, 1, 0)
  Q2 = predict(bst, X2) * vec2 
  vec3 = ifelse(A == 3, 1, 0)
  Q3 = predict(bst, X3) * vec3
  vec4 = ifelse(A == 4, 1, 0)
  Q4 = predict(bst, X4) * vec4
  
  deltaXGB = sum(Q1 + Q2 + Q3 + Q4) / N
  
  return(deltaXGB)
  
}
