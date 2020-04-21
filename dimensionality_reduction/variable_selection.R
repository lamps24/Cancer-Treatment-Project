########################################################################
# Variable Selection
# LASSO
# Initial dataset has n = 272 and p = 1570, so need to reduce size.
# Selection is done based on the variables impact on death.
# Metadata is automatically kept - LASSO only done on imaging variables
# Also turns three treatment indicators into one treatment factor with 8 levels
#
# inputs:
#   - original dataframe
# outputs:
#   - dataframe with only selected variables
#
########################################################################

library(glmnet)

# function
variable_select = function(original_df)
{
  
  # replace with own path
  original_df = read.csv("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Data/NKI_cleaned.csv")
  
  # subset data as matrices for cv.glmnet
  x = as.matrix(original_df[, 17:1570]) # independent vars excluding metadata
  y = as.factor(as.matrix(original_df[, 4])) # dependent var is death
  
  # lambda grid to search
  lambdas = 10^seq(3, -2, by = -.1)
  
  # cross-validation to search for best lambda
  cv = cv.glmnet(x, y, alpha = 1, family="binomial", lambda = lambdas) # alpha = 1 is LASSO
  
  # fit the model with the largest lambda w/in 1 se of minimized cv error
  m = glmnet(x, y, alpha = 1, family="binomial", lambda = cv$lambda.1se)
  
  # index of selected variables
  vars = m[["beta"]]@i + 1 # add one to the index if extracting vars from the original x matrix
  selected_vars = as.data.frame(x[, vars])
  
  # metadata
  metadata = original_df[, 3:14]
  
  # add an 8-level treatment factor
  metadata$treatment = ifelse(metadata$chemo==0 & metadata$amputation==0, 1, 0)
  metadata$treatment = ifelse(metadata$chemo==0 & metadata$amputation==1, 2, metadata$treatment)
  metadata$treatment = ifelse(metadata$chemo==1 & metadata$amputation==0, 3, metadata$treatment)
  metadata$treatment = ifelse(metadata$chemo==1 & metadata$amputation==1, 4, metadata$treatment)
  metadata$treatment = as.factor(metadata$treatment)
  
  # join metadata with selected vars
  newdata = cbind(metadata, selected_vars)
  
  return(list(newdata=newdata))

}
