########################################################################
# Genetic Algorithm
# Searches for the best values of eta that minimizes the value function
#
# inputs:
#   - df = dataset
#   - trt = column number of treatment being considered
#         Probably:
#           5 = chemo
#           6 = hormonal
#           7 = amputation
#   - M = population size, default is 100
#   - u = mutation rate, default is 0.3
#   - lam = number of offspring, default is 3
#   - gen = number of generations, default is 1
#   - val = name of value function being used
#   - var_list = index of vars that contribute to policy decision
#          9 = age
#          10 = diam
#          12 = grade
#          13 = lymphinfil
#          18 = NM_003430
#          21 = Contig_23211_RC
#          25 = NM_016359
#
# outputs:
#   - eta_best: best values for eta
#   - min_value: corresponding value for the best etas
#
########################################################################

library(data.table)
library(MASS)

gen_alg3 = function(df, trt, M=100, u=0.3, lam=3, gen=1, val_fun=OR_log, var_list=c(18, 21, 25)) {
  
  # initialize values at means
  p = length(var_list)
  mu = rep(0, p)
  i = 1
  for (var in var_list) {
    mu[i] = mean(df[, var_list[i]])
    i = i + 1
  }
  Sigma = matrix(0, nrow=p, ncol=p)
  diag(Sigma) = 1
  Sigma[1, 1] = 0.75^2 
  Sigma[2, 2] = 0.04^2 
  Sigma[3, 3] = 0.75^2
  
  eta = cbind(mvrnorm(M, mu, Sigma))
  V = rep(0,M)
  min_value = Inf
  for (j in 1:M){
    
    # convert eta values to policies - can adjust form of policy
    i = ifelse(df[, var_list[1]] < eta[j, 1] & 
                 df[, var_list[2]] > eta[j, 2] & 
                 df[, var_list[3]] > eta[j, 3], 1, 0)
    
    V[j] = OR_log(df, trt, i)$value
    if (V[j] < min_value){
      min_value = V[j]
      eta_best = eta[j,]
    }
  }
  
  
  offspring = 1 + rpois(n = M, lambda = lam) # can change how to generate num offsprings
  L = sum(offspring)
  Z = matrix(rep(0, L*p), ncol = p)
  Vtemp = matrix(rep(0, L), ncol = 1)
  l = 1
  for (r in 1:M) {
    for(s in 1:offspring[r]) {
      if (rnorm(1) < u & s>1) {
        Z[l, ] = eta[r, ] * 0.98 # changed from addition to multiplication (mutation rate = 0.98)
      }
      else {
        Z[l, ] = eta[r, ]
      }
      
      # convert eta values to policies - can adjust form of policy
      i = ifelse(df[, var_list[1]] < Z[l, 1] & 
                   df[, var_list[2]] > Z[l, 2] & 
                   df[, var_list[3]] > Z[l, 3], 1, 0)
      
      Vtemp[l] = OR_log(df, trt, i)$value # can change which value search function 
      l = l + 1
    }
  }
  tempvec = data.table(cbind(Vtemp, Z))
  tempvec = tempvec[order(tempvec$V1, decreasing = FALSE)] # want to sort by lowest values
  tempvec = tempvec[1:M]
  V = as.matrix(tempvec[, 1])
  num_col = p + 1
  eta = as.matrix(tempvec[, 2:num_col])
  if (V[1] < min_value){
    min_value = V[1]
    eta_best = eta[1,]
  }
  return(list(min_value=min_value, eta_best=eta_best))
}
