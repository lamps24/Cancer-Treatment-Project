########################################################################
# Outcome Regression (OR): XGBoost
# Computes the value using OR, using XGBoost to find the form of the Q function
#
# inputs:
#   - data = dataset
#   - trt = a scalar for which treatment is being used 
#      - chemo is column 5, and amputation is column 7
#   - i = policy vector, supplied from genetic algorithm
# outputs:
#   - deltaOR_log, the resulting value (scalar)
#
########################################################################

library("dplyr")
library("xgboost")

OR_xgb = function(df, trt, i){
  xdat = df
  N = length(xdat$age)
  
  # separating dataset into covariates and outcome (label)
  xdat = df[, c(trt, 8:ncol(df))]
  xdat = as.matrix(xdat)
  labels = df[, 2]
  
  # split the data into training and testing set
  index = sample(seq(1,272), 30) #can adjust the train/validation ratio
  train = xdat[-index, ]
  test = xdat[index, ]
  train_labels = labels[-index]
  test_labels = labels[index]
  
  # convert to XGBoost appropriate matrices
  dtrain = xgb.DMatrix(train, label = train_labels)
  dtest = xgb.DMatrix(test, label = test_labels)
  watchlist = list(train = dtrain, eval = dtest)
  
  # create list of parameters, and begin training the model
  param = list(max_depth = 2, eta = 1, verbose = 1, nthread = 2, objective = "binary:logistic", eval_metric = "auc")
  bst = xgb.train(param, dtrain, nrounds = 10, watchlist)
  
  # Create X0 matrix, where trt column is all 0's
  X0 = df[, c(trt, 8:ncol(df))]
  X0[trt] = 0
  X0 = as.matrix(X0)
  X0 = xgb.DMatrix(X0)
  
  # Create X1 matrix, where trt column is all 1's
  X1 = df[, c(trt, 8:ncol(df))]
  X1[trt] = 1
  X1 = as.matrix(X1)
  X1 = xgb.DMatrix(X1)
  
  # Make predictions with X0/X1 matrices
  Q0 = predict(bst, X0, type = 'response') * (1-i)
  Q1 = predict(bst, X1, type = 'response') * i
  
  value = sum(Q0 + Q1) / N
  
  return(list(value=value))
  
}
