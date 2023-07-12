# least square loss for the growth curve
get_loss_fun <- function(t_obs, pop_obs, t_pred, pop_pred){
  loss_fun <- function(param){
    t_scale <- param[1]
    #pop_scale <- param[2]
    t_offset <- param[2]
    interp <- approxfun(t_pred/t_scale, pop_pred)
    #mean_ <- interp(t_obs)
    #negloglik <- -dpois(pop_obs, mean_, log = TRUE)
    #negloglik[is.na(negloglik)] <- Inf
    #return(sum(negloglik))
    
    residual_sq <- ((pop_obs-interp(t_obs-t_offset))^2)
    residual_sq[is.na(residual_sq)] <- Inf
    residual_sq[t_obs-t_offset < 0] <- (pop_obs[t_obs-t_offset < 0])^2
    return(sqrt(mean(residual_sq)))
  }
  return(loss_fun)
}

library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(raster)
source("./R/rotate.R")
wi_landscape <- raster("./Wolf_rasters/Wolf_WI.tif")
sourceCpp("src/reaction_diffusion.cpp")

# get empirical data, without harvest
load("../new_Analysis/wolf.RData")
pop_obs <- wolf$Winter.Minimum.Count[-c(34:41)]
t_obs <- wolf$year[-c(34:41)] - 1979

wi_landscape_mat <- matrix( values( wi_landscape), 
                            wi_landscape@nrows, 
                            wi_landscape@ncols, byrow = T)

wi_landscape_mat[is.na(wi_landscape_mat)] <- 0
wi_landscape_mat <- wi_landscape_mat/max(wi_landscape_mat)
wi_landscape_mat <- rotate(wi_landscape_mat, 270)
image(wi_landscape_mat)
bigK <- sum(wi_landscape_mat>1e-1) * 0.0254


#### a grid search ####
exponents <- seq(0.75,2,0.25)
rs <- 10^(-exponents)
thins <- (exponents * 2 - 1) * 10
thins[thins<=0] <- 10
n_save <- 300
losses <- exponents + Inf

for(i in 1:length(exponents)){
  den <- 0 * wi_landscape_mat
  den[floor(0.11*nrow(den)), floor(0.83*ncol(den))] <- .01
  
  r <- rs[i]
  n <- n_save * thins[i]
  wi_land_res <- simulate_reaction_diffusion(den, wi_landscape_mat, r, n, thins[i])
  
  population <- apply(wi_land_res$density,3,
                      function(w,wi_landscape_mat){
                        sum(w)/sum(wi_landscape_mat)
                        },wi_landscape_mat)
  population <- population * bigK
  t_pred <- wi_land_res$time
  
  rm(wi_land_res)
  gc()
  
  loss_fun <- get_loss_fun(t_obs, pop_obs, t_pred, population)
  scale_est <- optim(c(1, 1),loss_fun )
  losses[i] <- scale_est$value
  plot(t_pred/scale_est$par[1] + scale_est$par[2]+ 1979, 
       population, type = "l", 
       xlab = "Year", ylab = "Population", xlim = c(1979, 2043))
  points(t_obs + 1979, pop_obs)
  points(wolf$year[c(34:41)]
         ,wolf$Winter.Minimum.Count[c(34:41)]
         ,pch = 7
  )
  points(wolf$year[c(36:41)]-4
         ,wolf$Winter.Minimum.Count[c(36:41)]
         ,pch = 7,col = "red"
  )
  text(x = wolf$year[34:35], 
       y = wolf$Winter.Minimum.Count[c(34:35)], labels = c("*","*"), 
       pos = 1, offset = 0.3)
  
  #legend("topleft",legend = c("observed:pre-hunting","observed:post-hunting",
  #                            "shift 2015-2020 to 2011",
  #                            "reaction-diffusion"), 
  #       lty = c(NA,NA,NA,1), pch = c(1,7,7,NA), col = c("black","black","red","black"),
  #       bg = "white")
  cat(i, losses[i])
}


#### the best fit ###
r <- 10^(-1.25)
thin <- 15
n <- n_save * thin
wi_land_res <- simulate_reaction_diffusion(den, wi_landscape_mat, r, n, thins[i])

population <- apply(wi_land_res$density,3,
                    function(w,wi_landscape_mat){
                      sum(w)/sum(wi_landscape_mat)
                    },wi_landscape_mat)
population <- population * bigK
t_pred <- wi_land_res$time

loss_fun <- get_loss_fun(t_obs, pop_obs, t_pred, population)
scale_est <- optim(c(1, 1),loss_fun )

t_pred_scale <- t_pred/scale_est$par[1] + scale_est$par[2]+ 1979


range_time <- apply(wi_land_res$density,3,function(w,wi_landscape_mat){sum(w>1e-1)},wi_landscape_mat)
area_pop<- data.frame(Population=wolf$Winter.Minimum.Count, Range=wolf_range$Winter.Minimum.Count, area = "WI", Year = wolf$year)



pdf("Figs/wi_with_fit.pdf",width = 11/1.2,height = 6/1.2)
par(mfcol = c(2,3),mar = c(2,2,2,2)+.3, mgp = c(1.3, 0.5, 0))
plot(t_pred/scale_est$par[1] + scale_est$par[2]+ 1979, 
     population, type = "l", 
     xlab = "Year", ylab = "Population", xlim = c(1979, 2043), main = "a)", adj = 0)
points(t_obs + 1979, pop_obs)
points(wolf$year[c(34:41)]
       ,wolf$Winter.Minimum.Count[c(34:41)]
       ,pch = 7
)
points(wolf$year[c(36:41)]-4
       ,wolf$Winter.Minimum.Count[c(36:41)]
       ,pch = 7,col = "red"
)
text(x = wolf$year[34:35], 
     y = wolf$Winter.Minimum.Count[c(34:35)], labels = c("*","*"), 
     pos = 1, offset = 0.3)

legend("topleft",legend = c("observed:pre-hunting","observed:post-hunting",
                            "shift 2015-2020 to 2011",
                            "Fisher-KPP"), 
       lty = c(NA,NA,NA,1), pch = c(1,7,7,NA), col = c("black","black","red","black"),
       bg = "white")

plot(range_time, population, type= "l", xlab = "Range", 
     ylab = "Population", main = "b)", adj = 0)
points(Population~Range,area_pop)
legend("topleft",legend = c("observed",
                            "Fisher-KPP"), 
       lty = c(NA,1), pch = c(1,NA), col = c("black","black"),
       bg = "white")

image(wi_land_res$density[,,9],main = "c) relative density 2000", adj = 0)
image(wi_land_res$density[,,20],main = "d) relative density 2010", adj = 0)
image(wi_land_res$density[,,28],main = "e) relative density 2020", adj = 0)
image(wi_land_res$density[,,37],main = "f) relative density 2030", adj = 0)
dev.off()
