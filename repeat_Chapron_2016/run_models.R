load("./data/chapron_treves_2016.Rdata")
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

WI_data <- list(N_year = 18,
                min_count = wisconsin$Nobs.Min,
                max_count = wisconsin$Nobs.Max,
                H_t = wisconsin$H,
                D_t = treatment$days_cull_WI)

initial_value <- list(N_t = (wisconsin$Nobs.Max),
                      psi_t_min = wisconsin$Nobs.Min,
                      psi_t_max = wisconsin$Nobs.Max,
                      sigma_proc = 0.06,
                      gamma = 1.06,
                      beta_0 = 0.16,
                      beta_1 = 0,
                      linvK = log(max(wisconsin$Nobs.Max)),
                      o_min = 0.97,
                      o_max = 1.03,
                      sigma_obs_min = 3.82,
                      sigma_obs_max = 4.72,
                      mu_0 = log(wisconsin$Nobs.Max[1]))

exp_null <- stan_model("./stan//exp_null.stan", 
            model_name = "exp_null")


sample_exp_null <- sampling(exp_null, 
                            data = WI_data,
                            control = list(adapt_delta = 0.9995, max_treedepth=30), 
                            iter = 15000,chains = 2,
                            init = list(chain1 = initial_value,chain2 = initial_value),
                            pars = c("sigma_proc","gamma","beta_0","o_min","o_max","sigma_obs_min","sigma_obs_max"), 
                            include = T)

exp_policy <- stan_model("./stan//exp_policy.stan", 
                       model_name = "exp_policy")


sample_exp_policy <- sampling(exp_policy, 
                            data = WI_data,
                            control = list(adapt_delta = 0.9995, max_treedepth=30), 
                            iter = 15000,chains = 2,
                            init = list(chain1 = initial_value, chain2 = initial_value),
                            pars = c("sigma_proc","gamma","beta_0","beta_1","o_min","o_max","sigma_obs_min","sigma_obs_max"), 
                            include = T)


