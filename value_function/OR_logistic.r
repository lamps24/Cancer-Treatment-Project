########################################################################
# Outcome Regression (OR): Logistic 
# Computes the value using OR, usinglogistic regression to find the form of the Q function
#
# inputs:
#   - data = dataset
#   - trt = a scalar for which treatment is being used 
#      - chemo is column 5, and amputation is column 7
#   - i = policy vector, supplied from genetic algorithm
#
# outputs:
#   - deltaOR_log, the resulting value (scalar)
#
########################################################################

library("dplyr")

OR_log = function(df, trt, i){
  
  # Assign inputs to different var names
  newdat = df[, c(2, trt, 8:ncol(df))] # remove the treatment columns, extra response columns
  N = nrow(newdat)
  
  # Q function if treatment is chemo 
  if (trt == 5){
    mod = glm(formula = eventdeath ~ chemo + age + diam + grade + lymphinfil + 
                NM_000926 + NM_012067 + NM_003430 + Contig23211_RC + NM_016109 + 
                NM_016359 + NM_001109, family = "binomial", data = newdat)
  }
  
  # Q function if treatment is amputation 
  if (trt ==7){
    mod = glm(formula = eventdeath ~ amputation + age + diam + posnodes + 
                grade + lymphinfil + NM_000926 + NM_012067 + NM_003430 + 
                NM_006096 + Contig23211_RC + AL049265 + NM_016359 + NM_001109 + 
                NM_001333 + amputation:posnodes + amputation:lymphinfil + 
                amputation:NM_012067 + amputation:NM_003430 + amputation:AL049265 + 
                amputation:NM_001333, family = "binomial", data = newdat)
  }

  # Instantiate new matrices for each treatment
  X0 = df
  X0[trt] = 0
  X1 = df
  X1[trt] = 1

  # Will multiply the predictions with the policy vector i
  Q0 = predict(mod, X0, type = 'response') * (1-i)
  Q1 = predict(mod, X1, type = 'response') * i 

  value = sum(Q0 + Q1) / N
  
  return(list(value=value))
}
