##  This file fits a logistic growth to the whole population counts and perform likelihood ratio test on whether there is a K

# logistic growth function, reduce to exponential when invK=0, i.e. K=\infty
logistic_growth <- function(t,r,p0,invK=0){
  p0*exp(r*t)/(1+(p0*invK)*(exp(r*t)-1))
}

# with Poisson noise
negloglik_pois <- function(pars, y,t,dd = TRUE){
  pars <- exp(pars)
  if(dd){
    pred <- logistic_growth(t,pars[1],pars[2],pars[3] )
  }
  else{
    pred <- logistic_growth(t,pars[1],pars[2],0)
  }
  -sum(dpois(y,pred,log = TRUE))
}

# fit the Pois model
fit_pois_mod <- function(y,t,t_new,t0 = 1987,dd = TRUE, se = T, pred_inter = T){
  #opt <- optim(c(log(max(y)/min(y))/length(y),log(min(y)),-log(max(y))),negloglik_pois,y = y,t = t-t0, dd = dd, hessian = T)
  #pars <- exp(opt$par)
  if(dd){
    opt <- optim(c(log(max(y)/min(y))/length(y),log(min(y)),-log(max(y))),negloglik_pois,y = y,t = t-t0, dd = dd, hessian = T)
    pars <- exp(opt$par)
    pred <- logistic_growth(t_new-t0,pars[1],pars[2],pars[3] )
    if(se){
      the_se <- lapply(1:3, function(i, pars, t_new, pred){
        par1 <- pars
        par1[i] <- par1[i] + min(1e-8*abs(par1[i]),1e-8)
        pred1 <- logistic_growth(t_new,exp(par1[1]),exp(par1[2]),exp(par1[3]))
        (pred1-pred)/(par1[i]-pars[i])
      }, opt$par, t_new-t0, pred) |>
        Reduce(f = cbind) |> 
        apply(1,function(vec, hess){
          (as.vector(t(vec) %*% solve(hess, (vec))))
        }, opt$hessian) 
        the_se <- sqrt(the_se+(pred*pred_inter) )
    }
    else the_se <- NA
  }
  else{
    opt <- optim(c(log(max(y)/min(y))/length(y),log(min(y))),negloglik_pois,y = y,t = t-t0, dd = dd, hessian = T)
    pars <- exp(opt$par)
    pred <- logistic_growth(t_new-t0,pars[1],pars[2],0)
    # calculate Wald se using brutal force 
    if(se){
      the_se <- lapply(1:2, function(i, pars, t_new, pred){
        par1 <- pars
        par1[i] <- par1[i] + min(1e-8*abs(par1[i]),1e-8)
        pred1 <- logistic_growth(t_new,exp(par1[1]),exp(par1[2]),0)
        (pred1-pred)/(par1[i]-pars[i])
      }, opt$par, t_new-t0, pred) |>
        Reduce(f = cbind) |> 
        apply(1,function(vec, hess){
          (as.vector(t(vec) %*% solve(hess, (vec))))
        }, opt$hessian) 
      the_se <- sqrt(the_se+(pred*pred_inter) )
    }
    else the_se <- NA
  }
  
  return(list(opt = opt, pred = pred, se = the_se))
}

negloglik_gaus <- function(pars, y,t,dd = TRUE){
  pars <- exp(pars)
  if(dd){
    pred <- logistic_growth(t,pars[1],pars[2],pars[3] )
  }
  else{
    pred <- logistic_growth(t,pars[1],pars[2],0)
  }
  -sum(dnorm(y,pred,pars[3+dd],log = TRUE))
}

fit_norm_mod <- function(y,t,t_new,t0 = 1987,dd = TRUE,se = T, pred_inter = T){
  
  if(dd){
    opt <- optim(c(log(max(y)/min(y))/length(y),log(min(y)),-log(max(y)),log(var(y))/2),negloglik_gaus,y = y,t = t-t0, dd = dd, hessian = T)
    pars <- exp(opt$par)
    pred <- logistic_growth(t_new-t0,pars[1],pars[2],pars[3] )
    
    if(se){
      the_se <- lapply(1:4, function(i, pars, t_new, pred){
        par1 <- pars
        par1[i] <- par1[i] + min(1e-8*abs(par1[i]),1e-8)
        pred1 <- logistic_growth(t_new,exp(par1[1]),exp(par1[2]),exp(par1[3]))
        (pred1-pred)/(par1[i]-pars[i])
      }, opt$par, t_new-t0, pred) |>
        Reduce(f = cbind) |> 
        apply(1,function(vec, hess,samp_var,pred_inter){
          sqrt((t(vec) %*% solve(hess, (vec)))+samp_var * pred_inter)
        }, opt$hessian, pars[4]^2, pred_inter)
    }
    else the_se <- NA
    
  }
  else{
    opt <- optim(c(log(max(y)/min(y))/length(y),log(min(y)),log(1000)),negloglik_gaus,y = y,t = t-t0, dd = dd, hessian = T, control = list(reltol = 1e-10))
    pars <- exp(opt$par)
    pred <- logistic_growth(t_new-t0,pars[1],pars[2],0)
    
    if(se){
      the_se <- lapply(1:3, function(i, pars, t_new, pred){
        par1 <- pars
        par1[i] <- par1[i] + min(1e-8*abs(par1[i]),1e-8)
        pred1 <- logistic_growth(t_new,exp(par1[1]),exp(par1[2]),0)
        (pred1-pred)/(par1[i]-pars[i])
      }, opt$par, t_new-t0, pred) |>
        Reduce(f = cbind) |> 
        apply(1,function(vec, hess,samp_var,pred_inter){
          sqrt(t(vec) %*% solve(hess, (vec))+samp_var * pred_inter)
        }, opt$hessian, pars[3]^2, pred_inter)
    }
    else the_se <- NA
    
  }
  
  return(list(opt = opt, pred = pred, se = the_se))
}


# logistic vs exponential, for sensitivity, we can exclude the harvest year or all year after them
logistic_mod <- fit_pois_mod(wolf$Winter.Minimum.Count[-c(34:41)]
                             , wolf$year[-c(34:41)]
                             , wolf$year
                             )
exponential_mod <- fit_pois_mod(wolf$Winter.Minimum.Count[-c(34:41)]
                                , wolf$year[-c(34:41)]
                                , wolf$year
                                , dd = FALSE)

# see fitting


#pdf("./figs/population_overall.pdf", width = 6, height = 3.5)
png("./figs/population_overall.png", width = 6, height = 3.5, res = 500, unit = "in")
par(mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))
plot(wolf$year
     ,logistic_mod$pred, type = "l", ylab = "Population", xlab = 'Year', ylim = c(0,1300))


polygon(x = c(wolf$year, rev(wolf$year)),
        y = c(logistic_mod$pred - qnorm(0.975)*logistic_mod$se, 
              rev(logistic_mod$pred + qnorm(0.975)*logistic_mod$se)),
        col =  adjustcolor("black", alpha.f = 0.10), border = NA)


points(wolf$year[-c(34:41)]
       ,wolf$Winter.Minimum.Count[-c(34:41)]
       )
points(wolf$year[c(34:41)]
       ,wolf$Winter.Minimum.Count[c(34:41)]
       ,pch = 7
)

points(wolf$year[c(36:41)]-4
       ,wolf$Winter.Minimum.Count[c(36:41)]
       ,pch = 7,col = "red"
)

lines(wolf$year
      ,exponential_mod$pred
      , lty = 2)
#points(wolf$year,wolf$Winter.Minimum.Count)
legend("topleft",legend = c("observed:pre-hunting","observed:post-hunting",
                            "move 2015-2020 to 2011",
                            "logistic","exponential"), 
       lty = c(NA,NA,NA,1,2), pch = c(1,7,7,NA,NA), col = c("black","black","red","black","black"))

dev.off()

# likelihood ratio test
likelihood_ratio <- 2*(exponential_mod$opt$value-logistic_mod$opt$value)
pchisq(likelihood_ratio,1, lower.tail = F, log.p = T)

# is range expansion also logistic?
logistic_range <- fit_norm_mod(wolf_range$Winter.Minimum.Count[-c(34:41)]
                               , wolf_range$year[-c(34:41)]
                               , wolf_range$year
                               )
exponential_range <- fit_norm_mod(wolf_range$Winter.Minimum.Count[-c(34:41)]
                                  , wolf_range$year[-c(34:41)]
                                  , wolf_range$year
                                  ,dd = F, se = F)

png("./figs/range_overall.png", width = 6, height = 3.5, res = 500, unit = "in")
par(mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))
plot(wolf_range$year
     ,logistic_range$pred, type = "l", ylab = "Range", xlab = 'Year', ylim = c(0,1.3 * max(wolf_range$Winter.Minimum.Count)))

polygon(x = c(wolf$year, rev(wolf$year)),
        y = c(logistic_range$pred - qnorm(0.975)*logistic_range$se, 
              rev(logistic_range$pred + qnorm(0.975)*logistic_range$se)),
        col =  adjustcolor("black", alpha.f = 0.10), border = NA)

points(wolf_range$year[-c(34:41)]
       ,wolf_range$Winter.Minimum.Count[-c(34:41)]
       )
points(wolf_range$year[c(34:41)]
       ,wolf_range$Winter.Minimum.Count[c(34:41)]
       ,pch = 7
)
points(wolf_range$year[c(36:41)]-4
       ,wolf_range$Winter.Minimum.Count[c(36:41)]
       ,pch = 7
       ,col = "red"
)
lines(wolf_range$year
      ,exponential_range$pred
      ,lty = 2
     
      )
legend("topleft",legend = c("observed:pre-hunting","observed:post-hunting",
                            "move 2015-2020 to 2011",
                            "logistic","exponential"), 
       lty = c(NA,NA,NA,1,2), pch = c(1,7,7,NA,NA), col = c("black","black","red","black","black"))

dev.off()

likelihood_ratio_range <- 2*(exponential_range$opt$value-logistic_range$opt$value)
pchisq(likelihood_ratio_range,1, lower.tail = F, log.p = T)


# what's the predicted range capacity? and CI
1/exp(logistic_range$opt$par[3])
1/exp(logistic_range$opt$par[3]+1.96*sqrt((solve(logistic_range$opt$hessian))[3,3]))
1/exp(logistic_range$opt$par[3]-1.96*sqrt((solve(logistic_range$opt$hessian))[3,3]))

# is it coincide with the population level K and the observed equalibrium density?
last_decade <- data.frame(year = 2000:2020, density = (wolf$Winter.Minimum.Count/wolf_range$Winter.Minimum.Count)[year>=2000])

## population K
1/exp(logistic_mod$opt$par[3])
1/exp(logistic_mod$opt$par[3]+1.96*sqrt((solve(logistic_mod$opt$hessian))[3,3]))
1/exp(logistic_mod$opt$par[3]-1.96*sqrt((solve(logistic_mod$opt$hessian))[3,3]))

1/exp(logistic_mod$opt$par[3])/mean(last_decade$density)
1/exp(logistic_mod$opt$par[3]+1.96*sqrt((solve(logistic_mod$opt$hessian))[3,3]))/mean(last_decade$density)
1/exp(logistic_mod$opt$par[3]-1.96*sqrt((solve(logistic_mod$opt$hessian))[3,3]))/mean(last_decade$density)


