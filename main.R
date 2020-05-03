# modify as needed
source("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Code/variable_selection.R")
source("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Code/ipw.R")
source("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Code/aipw.R")
source("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Code/or_logistic.R")
source("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Code/or_xgb.R")
source("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Code/policy_search1.R")
source("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Code/policy_search2.R")
source("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Code/policy_search3.R")
original_df = read.csv("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Data/NKI_cleaned.csv")

# set seed
set.seed(1234) #1234 results in 29 vars, 5080 results in 43

# select variables
df = variable_select(original_df)$newdata

# genetic algorithm
# var_list:
#   9 = age
#   10 = diam
#   12 = grade
#   13 = lymphinfil
#   18 = NM_003430
#   21 = Contig_23211_RC
#   25 = NM_016359
#   * note that the three imaging vars can change column indices based on seed - these are based on seed 1234
#test = gen_alg(df, trt=5, M=100, u=0.3, lam=3, gen=1, val_fun=OR_log, var_list=c(18, 21, 25))


# results matrix
results = array(NA, dim = c(24, 11))
colnames(results) = c('Treatment', 'Value Function', ' Policy Form', 'Best Value', 'Age <', 'Diam >', 'Grade >', 'Lymphinfil >',
                      'NM_003430 <', 'Contig23211_RC >', 'NM_016359 >')

# loop parameters
value_functions = c(OR_log, OR_xgb, ipw, aipw)
vars1 = c(9, 10, 12, 13)
vars2 = c(9, 10, 12, 13, 18, 21, 25)
vars3 = c(18, 21, 25)

i = 1
for (fun in value_functions) { 
  for (treatment in c(5, 7)) {
    for (policy_form in c(1, 2, 3)) {
      
      results[i, 1] = treatment
      results[i, 2] = 'fun'
      results[i, 3] = policy_form

      if (policy_form == 1) {
        temp = gen_alg1(df, trt=treatment, M=100, u=0.3, lam=3, gen=1, val_fun=fun, var_list=vars1)
        results[i, 4] = temp$min_value
        results[i, 5:8] = temp$eta_best
      }

      if (policy_form == 2) {
        temp = gen_alg2(df, trt=treatment, M=100, u=0.3, lam=3, gen=1, val_fun=fun, var_list=vars2)
        results[i, 4] = temp$min_value
        results[i, 5:11] = temp$eta_best
      }
      
      if (policy_form == 3) {
        temp = gen_alg3(df, trt=treatment, M=100, u=0.3, lam=3, gen=1, val_fun=fun, var_list=vars3)
        results[i, 4] = temp$min_value
        results[i, 9:11] = temp$eta_best
      }
      
      i = i + 1
      
    }
  }
}


