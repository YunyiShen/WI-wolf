data {
    int<lower = 0> T; // number of years
    int<lower = 0> min_count[T]; //winter minimum count
    int<lower = 0> max_count[T]; //winter maximum count
    int<lower = 0> H_t[T]; //harvest
    int<lower = 0> D_t[T]; //days culled
}

parameters {
    real<lower = 0> sigma_proc;
    real<lower = 0> gamma;
    real beta_0;
    //real beta_1; //the "policy signal", useful in non-null's 
    //real linvK; // inverse carrying capacity
    real<lower = 0, upper = 1> o_min; //detection rate for min
    real<lower = 1, upper = 10> o_max; //detection multiplier for max
    real<lower = 0> sigma_obs; // observ error
    real mu_0; //baseline year 
}

model {
    // some useful local variables
    int<lower = 0> N_t[T]; // true population
    real<lower = 0> mu_t; // expected population
    real r_t;// population growth rate
    real<lower = 0> psi_t; // expected observed pop

    // priors


    // likelihood

    N_t[1] ~ lognormal(mu_0,sigma_proc);

    for(tt in 1:T){
        // this years observation 
        target += gamma_lpdf(psi_t | (N_t[tt]/sigma_obs)^2, N_t[tt]/((sigma_obs)^2) );
        target += poisson_lpmf(min_count[tt] | o_min * psi_t);
        target += poisson_lpmf(max_count[tt] | o_max * psi_t);
        // things about grwoth 
        if(tt <= T){
            r_t = beta_0; // growth rate
            //r_t = beta_0 + beta_1 * D_t[tt];
            mu_t = log(N_t[tt] * exp(r_t) - gamma * H_t[tt]); // exponential grwoth
            // mu_t = log(N_t[tt] * exp(r_t) * (1-exp(linvK) * N_t[tt]) - gamma * H_t[tt]) // logistic
            N_t[tt+1] ~ lognormal(mu_t, sigma_proc);
        }

    }

}


