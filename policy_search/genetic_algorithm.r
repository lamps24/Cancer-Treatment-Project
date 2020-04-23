########################################################################
# Genetic Algorithm
# Searches for the best values of eta that maximizes the value function (currently OR using logstic regression)
#
# inputs:
#   - data = dataset
#   - M = population size, default is 100
#   - u = mutation rate, default is 0.3
#   - lam = number of offspring, default is 3
#   - gen = number of generations, default is 1
#
# outputs:
#   - eta_best: best values for eta
#   - max_value: corresponding value for the best etas
#
########################################################################

library(data.table)

gen_alg = function(df, M = 100, u = 0.3, lam = 3, gen = 1){
  for (k in 1:gen){
    set.seed(5080)
    X = matrix(rnorm(M*2), nrow = M) #right now we only need 2 columns for eta
    X[,1] = abs(30*X[, 1]) #all of this is just arbitrary and can be modified
    X[,2] = abs(3*X[, 2])
    V = rep(0,M)
    max_value = -Inf
    for (j in 1:M){
      V[j] = OR_log(dat, X[j,])
      if (V[j] > max_value){
        max_value = V[j]
        eta_best = X[j,]
      }
    }
    
    offspring = 1 + rpois(n = M, lambda = lam) # can change how to generate num offsprings
    L = sum(offspring)
    Z = matrix(rep(0, L*2), ncol = 2)
    Vtemp = matrix(rep(0, L), ncol = 1)
    l = 1
    for (r in 1:M){
      for(s in 1:offspring[r]){
        if (rnorm(1) < u & s>1){
          Z[l, ] = X[r,] + 20*rnorm(2) # can change how to mutate offspring
        }
        else{
          Z[l, ]= X[r,]
        }
        Vtemp[l] = OR_log(dat, Z[l, ]) # can change which value search function 
        l = l + 1
      }
    }
    tempvec = data.table(cbind(Vtemp, Z))
    tempvec = tempvec[order(tempvec$V1, decreasing = TRUE)]
    tempvec = tempvec[1:M]
    V = as.matrix(tempvec[, 1])
    X = as.matrix(tempvec[, 2:3])
    if (V[1] > max_value){
      max_value = V[1]
      eta_best = X[1,]
    }
  }
  return(c(max_value, eta_best))
}
