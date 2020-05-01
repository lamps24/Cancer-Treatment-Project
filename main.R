# modify as needed
source("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Code/variable_selection.R")
source("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Code/ipw.R")
source("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Code/aipw.R")
source("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Code/or_logistic.R")
source("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Code/or_xgb.R")
source("C:/Users/lamps/Documents/Class/IE 5080 - Personalized Medicine/Project/Code/policy_search2.R")
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
test = gen_alg(df, trt=7, M=100, u=0.3, lam=3, gen=1, val_fun=ipw, var_list=c(9, 10, 12, 13, 18, 21, 25))
