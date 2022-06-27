load("./data/chapron_treves_2016.Rdata")
library(rstan)
library(bridgesampling)

exp_null <- stan_model("./stan//exp_null.stan", 
                       model_name = "exp_null")

exp_policy <- stan_model("./stan//exp_policy.stan", 
                         model_name = "exp_policy")

logistic_null <- stan_model("./stan//logistic_null.stan", 
                            model_name = "logistic_null")

logistic_policy <- stan_model("./stan//logistic_policy.stan", 
                            model_name = "logistic_policy")


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

set.seed(42)

WI_data <- list(N_year = 18,
                min_count = wisconsin$Nobs.Min,
                max_count = wisconsin$Nobs.Max,
                H_t = wisconsin$H,
                D_t = treatment$days_cull_WI/365)



initial_value <- list(N_t = (wisconsin$Nobs.Max),
                      psi_t_min = wisconsin$Nobs.Min,
                      psi_t_max = wisconsin$Nobs.Max,
                      sigma_proc = 0.06,
                      gamma = 1.06,
                      beta_0 = 0.16,
                      beta_1 = -.5,
                      K = (max(wisconsin$Nobs.Max)),
                      o_min = 0.97,
                      o_max = 1.03,
                      sigma_obs = 3.82
                      )



## Take samples

sample_exp_null <- sampling(exp_null, 
                            data = WI_data,
                            control = list(adapt_delta = 0.999, max_treedepth=30), 
                            iter = 15000,chains = 4, thin = 5, 
                            init = list(initial_value,initial_value,initial_value,initial_value),
                            )#pars = c("sigma_proc","gamma","beta_0","o_min","o_max","sigma_obs"), 
                            #include = T)


sample_exp_policy <- sampling(exp_policy, 
                            data = WI_data,
                            control = list(adapt_delta = 0.999, max_treedepth=30), 
                            iter = 15000,chains = 4,thin = 5,
                            init = list(initial_value,initial_value,initial_value,initial_value),
                            )#pars = c("sigma_proc","gamma","beta_0","beta_1","o_min","o_max","sigma_obs"), 
                            #include = T)


sample_logistic_null <- sampling(logistic_null, 
                                 data = WI_data,
                                 control = list(adapt_delta = 0.999, max_treedepth=30), 
                                 iter = 15000,chains = 4,thin = 5,
                                 init = list(initial_value,initial_value,initial_value,initial_value))#,
#pars = c("sigma_proc","gamma","beta_0","o_min","o_max","sigma_obs","K"), 
#include = T)

sample_logistic_policy <- sampling(logistic_policy, 
                                 data = WI_data,
                                 control = list(adapt_delta = 0.999, max_treedepth=30), 
                                 iter = 15000,chains = 4,thin = 5,
                                 init = list(initial_value,initial_value,initial_value,initial_value)
)#pars = c("sigma_proc","gamma","beta_0","o_min","o_max","sigma_obs"), 
#include = T)


## BF for inference 

normalizing_exp_policy <- bridge_sampler(sample_exp_policy)
normalizing_exp_null <- bridge_sampler(sample_exp_null)

normalizing_logistic_policy <- bridge_sampler(sample_logistic_policy)
normalizing_logistic_null <- bridge_sampler(sample_logistic_null)

bayes_factor_exp_policy <- bridgesampling::bf(normalizing_exp_policy, normalizing_exp_null)
bayes_factor_logistic_null <- bridgesampling::bf(normalizing_logistic_null, normalizing_exp_null)
bayes_factor_logistic_policy <- bridgesampling::bf(normalizing_logistic_policy, normalizing_exp_null)


bayes_factor_logistic_vs_exp_policy <- bridgesampling::bf(normalizing_logistic_null, normalizing_exp_policy)

save.image("./model_res.RData")




