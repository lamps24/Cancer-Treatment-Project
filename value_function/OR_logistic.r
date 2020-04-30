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
  # Currently there are a few different options of models for how to best include interactions with treatment
  if (trt == 5){
    #mod = glm(eventdeath ~ .*chemo, data = newdat, family = "binomial")
    #mod = glm(eventdeath ~ chemo + (histtype + grade + angioinv + lymphinfil + age + posnodes + diam)*chemo, data = newdat, family = "binomial")
    mod = glm(eventdeath ~ chemo + (histtype + grade + angioinv + lymphinfil + age + posnodes + diam+NM_012067
                                    + NM_021069 +NM_005132 +AL117638 + NM_006398+NM_004950 +NM_006461 +Contig23211_RC+NM_016109
                                    +Contig42011_RC+NM_001085 +NM_001109+NM_002274)*chemo , data = newdat, family = 'binomial')
    #mod = glm(eventdeath ~ chemo+(histtype + grade + angioinv + lymphinfil + age + posnodes + diam+ NM_001627+NM_012067+Contig21225_RC+NM_003430+NM_021069+NM_005132+AL117638+NM_006398+NM_004950+Contig23211_RC+NM_016359+L27560+NM_001109+M24895)*chemo, data = newdat, family = "binomial")
  }
  
  # Q function if treatment is amputation 
  if (trt ==7){
    mod = glm(eventdeath ~ .*amputation, data = newdat, family = "binomial", control = list(maxit = 50)) # without control, says algorithm doesn't converge
    #mod = glm(eventdeath ~ amputation + (histtype + grade + angioinv + lymphinfil + age + posnodes + diam+ NM_001627+NM_012067+Contig21225_RC+NM_003430+NM_021069+NM_005132+AL117638+NM_006398+NM_004950+Contig23211_RC+NM_016359+L27560+NM_001109+M24895+NM_002274)*amputation, data = newdat, family = "binomial")
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
  
  return(deltaOR_log)
}
