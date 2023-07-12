library(Matrix)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("src/reaction_diffusion.cpp")

Landscape = matrix(1, 100,100)
X_init <- matrix(0,100,100)
X_init[51,1] <- .1
r <- .02
n <- 500

### tunnel induce Allee effect
allee_landscape <- Landscape
allee_landscape[c(1:46,101-1:46), 1:50] <- 1e-7

allee_land_res <- simulate_reaction_diffusion(X_init, allee_landscape, r, n)
allee_land <- allee_land_res$density
allee_land_frontier <- allee_land_res$frontier
allee_land_growth <- allee_land_res$growth
range_time <- apply(allee_land,3,function(w,allee_landscape){sum(w>1e-5)/sum(allee_landscape>1e-5)},allee_landscape)
population <- apply(allee_land,3,function(w,allee_landscape){sum(w)/sum(allee_landscape)},allee_landscape)
frontier <- 1:n
for(i in 1:n){
  frontier[i] <- sum(((allee_land_growth[,,i]+allee_land_frontier[,,i])>=1e-5) * (allee_land_growth[,,i]>=1e-5))/(sum(allee_landscape>1e-5))
}

pdf("Figs/allee.pdf",width = 8,height = 6)
par(mfrow = c(3,3),mar = c(2,2,2,2)+.3, mgp = c(1.3, 0.5, 0))
image(allee_land[,,1],main = "a) density t=1", adj = 0)
image(allee_land[,,100],main = "b) density t=100", adj = 0)
plot(range_time, type = "l", xlab = "time", ylab = "relative range and population", 
     main = "c)", adj = 0)
lines(population, lty = 2)
legend("bottomright",legend = c("range","population"),lty = c(1,2))

image(allee_land[,,200],main = "d) t=200", adj = 0)
image(allee_land[,,300],main = "e) t=300", adj = 0)
plot(range_time,population,
     xlab = "relative range", 
     ylab = "relative population", main = "f) population-range curve", adj = 0)



image(((allee_land_growth[,,200]+allee_land_frontier[,,200])>=1e-5) * (allee_land_growth[,,200]>=1e-5),
      main = "g) frontier t=200", adj = 0)
image(((allee_land_growth[,,300]+allee_land_frontier[,,300])>=1e-5) * (allee_land_growth[,,300]>=1e-5),
      main = "h) frontier t=300", adj = 0)
plot(frontier, type = "l", xlab = "time", ylab = "relative range in frontier", 
     main = "i) range in frontier", adj = 0)


dev.off()

