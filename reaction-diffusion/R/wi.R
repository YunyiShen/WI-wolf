library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(raster)
source("./R/rotate.R")
wi_landscape <- raster("./Wolf_rasters/Wolf_WI.tif")
sourceCpp("src/reaction_diffusion.cpp")

wi_landscape_mat <- matrix( values( wi_landscape), 
                            wi_landscape@nrows, 
                            wi_landscape@ncols, byrow = T)

wi_landscape_mat[is.na(wi_landscape_mat)] <- 0
wi_landscape_mat <- wi_landscape_mat/max(wi_landscape_mat)
wi_landscape_mat <- rotate(wi_landscape_mat, 270)

image(wi_landscape_mat)

den <- 0 * wi_landscape_mat
den[floor(0.11*nrow(den)), floor(0.83*ncol(den))] <- 1
image(den)

#r <- .1
r <- 1e-3
n <- 1200 * 3

wi_land_res <- simulate_reaction_diffusion(den, wi_landscape_mat, r, n, 30)
wi_land <- wi_land_res$density
wi_land_frontier <- wi_land_res$frontier
wi_land_growth <- wi_land_res$growth
range_time <- apply(wi_land,3,function(w,wi_landscape_mat){sum(w>1e-1)/sum(wi_landscape_mat>1e-1)},wi_landscape_mat)
population <- apply(wi_land,3,function(w,wi_landscape_mat){sum(w)/sum(wi_landscape_mat)},wi_landscape_mat)
frontier <- 1:n
for(i in 1:n){
  frontier[i] <- sum(((wi_land_growth[,,i]+wi_land_frontier[,,i])>=1e-5) * (wi_land_growth[,,i]>=1e-1), na.rm = T)/(sum(wi_landscape_mat>1e-1))
}

pdf("Figs/wi.pdf",width = 8,height = 6)
#png("Figs/wi.png",width = 8,height = 6)
par(mfrow = c(3,3),mar = c(2,2,2,2)+.3, mgp = c(1.3, 0.5, 0))
image(wi_land[,,50],main = "a) density t=50", adj = 0)
image(wi_land[,,100],main = "b) density t=100", adj = 0)
plot(range_time, type = "l", xlab = "time", ylab = "relative range and population", 
     main = "c)", adj = 0)
lines(population, lty = 2)
legend("bottomright",legend = c("range","population"),lty = c(1,2))

image(wi_land[,,500],main = "d) t=500", adj = 0)
image(wi_land[,,1000],main = "e) t=1000", adj = 0)
plot(range_time,population,
     xlab = "relative range", 
     ylab = "relative population", main = "f) population-range curve", adj = 0)



image(((wi_land_growth[,,500]+wi_land_frontier[,,500])>=1e-5) * (wi_land_growth[,,500]>=1e-5),
      main = "g) frontier t=500", adj = 0)
image(((wi_land_growth[,,1000]+wi_land_frontier[,,1000])>=1e-5) * (wi_land_growth[,,1000]>=1e-5),
      main = "h) frontier t=1000", adj = 0)
plot(frontier, type = "l", xlab = "time", ylab = "relative range in frontier", 
     main = "i) range in frontier", adj = 0)


dev.off()
