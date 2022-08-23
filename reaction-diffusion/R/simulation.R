library(Matrix)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("src/reaction_diffusion.cpp")

Landscape = matrix(1, 100,100)
X_init <- matrix(0,100,100)
X_init[51,1] <- .1
r <- .02
n <- 500

## uniform landscape
uniform_land_res <- simulate_reaction_diffusion(X_init, Landscape, r, n)
uniform_land <- uniform_land_res$density
uniform_land_frontier <- uniform_land_res$frontier
uniform_land_growth <- uniform_land_res$growth
range_time <- apply(uniform_land,3,function(w,Landscape){sum(w>1e-5)/sum(Landscape>1e-5)},Landscape)
population <- apply(uniform_land,3,function(w,Landscape){sum(w)/sum(Landscape)},Landscape)
frontier <- 1:n
for(i in 1:n){
  frontier[i] <- sum(((uniform_land_growth[,,i]+uniform_land_frontier[,,i])>=1e-5) * (uniform_land_growth[,,i]>=1e-5))/(sum(Landscape>1e-5))
}

png("Figs/hardboundry.png",width = 8,height = 6, units = "in", res = 500)
par(mfrow = c(3,3),mar = c(2,2,2,2)+.3, mgp = c(1.3, 0.5, 0))
image(uniform_land[,,1],main = "t=1")
image(uniform_land[,,100],main = "t=100")
plot(range_time, type = "l", xlab = "time", ylab = "relative range and population")
lines(population, lty = 2)
legend("bottomright",legend = c("range","population"),lty = c(1,2))


image(((uniform_land_growth[,,200]+uniform_land_frontier[,,200])>=1e-5) * (uniform_land_growth[,,200]>=1e-5),main = "frontier t=200")
image(((uniform_land_growth[,,300]+uniform_land_frontier[,,300])>=1e-5) * (uniform_land_growth[,,300]>=1e-5),main = "frontier t=300")
plot(frontier, type = "l", xlab = "time", ylab = "relative range in frontier")

image(uniform_land[,,200],main = "t=200")
image(uniform_land[,,300],main = "t=300")
plot(range_time,population,xlab = "relative range", ylab = "relative population")
dev.off()

## with gradient

gradient_1d <- c(seq(0.1,1,length.out = 10),rep(1,80),seq(1,.1,length.out = 10))
Landscape2 <- as.matrix(gradient_1d)%*%t(as.vector(gradient_1d))
image(Landscape2)

core_land_res <- simulate_reaction_diffusion(X_init, Landscape2, r, n)
core_land <- core_land_res$density
core_land_frontier <- core_land_res$frontier
core_land_growth <- core_land_res$growth

range_time2 <- apply(core_land,3,function(w,Landscape){sum(w>1e-5)/sum(Landscape>1e-5)},Landscape2)
population2 <- apply(core_land,3,function(w,Landscape){sum(w)/sum(Landscape)},Landscape2)
frontier2 <- 1:n #apply(core_land_frontier,3,function(w){sum(w>=0.1*max(max(w),1e-10))/sum(w>-10)})
for(i in 1:n){
  frontier2[i] <- sum(((core_land_growth[,,i]+core_land_growth[,,i])>=1e-5) * (core_land_growth[,,i]>=1e-5))/(sum(Landscape>1e-5))
}

png("Figs/softboundry.png",width = 8,height = 6, units = "in", res = 500)
par(mfrow = c(3,3),mar = c(2,2,2,2)+.3, mgp = c(1.3, 0.5, 0))
image(core_land[,,1],main = "t=1")
image(core_land[,,100],main = "t=100")
plot(range_time2, main = "relative range/population", type = "l", xlab = "t", ylab = "relative range and population")
lines(population2, lty = 2)
abline(v = 100, lty = 3)
abline(v = 225, lty = 3)
abline(v = 350, lty = 3)
legend("bottomright",legend = c("range","population"),lty = c(1,2))

image(((core_land_growth[,,225]+core_land_frontier[,,225])>=1e-5) * (core_land_growth[,,225]>=1e-5),main = "frontier t=225")
image(((core_land_growth[,,350]+core_land_frontier[,,350])>=1e-5) * (core_land_growth[,,350]>=1e-5),main = "frontier t=350")
plot(frontier2, type = "l", xlab = "time", ylab = "relative range in frontier")

image(core_land[,,225],main = "t=225")
image(core_land[,,350],main = "t=350")
plot(range_time2,population2,xlab = "relative range", ylab = "relative population")
dev.off()

# patchy landscape
rmvn <- function(n, mu = 0, V = matrix(1)) {
  p <- length(mu)
  if (any(is.na(match(dim(V), p)))) 
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}

# Set up a square lattice region
simgrid <- expand.grid(1:100, 1:100)
nn <- nrow(simgrid)

# Set up distance matrix
distance <- as.matrix(dist(simgrid))
# Generate random variable
set.seed(1234)
phi <- 0.05
X <- rmvn(1, rep(0, nn), exp(-phi * distance))
Landscape31 <- matrix(X,100,100)
Landscape31 <- (Landscape31 * (Landscape31>=0))
Landscape31 <- apply(Landscape31,c(1,2),min,1)
image(Landscape31)
Landscape3 <- Landscape31
Landscape3[Landscape3==0] <- 0
X_init2 <- 0 * Landscape3
X_init2[1,100] <- .1
r2 <- .02
n <- 1000

patchy_land_res <- simulate_reaction_diffusion(X_init2, Landscape3, r2, n)
patchy_land <- patchy_land_res$density
patchy_land_frontier <- patchy_land_res$frontier
patchy_land_growth <- patchy_land_res$growth

image(((patchy_land_growth[,,500]+patchy_land_frontier[,,500])>=1e-5) * (patchy_land_growth[,,500]>=1e-5))

range_time3 <- apply(patchy_land,3,function(w,Landscape){sum(w>1e-5)/sum(Landscape>1e-5)},Landscape3)
population3 <- apply(patchy_land,3,function(w,Landscape){sum(w)/sum(Landscape)},Landscape3)
#frontier3 <- apply(patchy_land_frontier,3,function(w){sum(w>=0.1*max(max(w),1e-10))/sum(w>-10)})

frontier3 <- 1:n #apply(core_land_frontier,3,function(w){sum(w>=0.1*max(max(w),1e-10))/sum(w>-10)})
for(i in 1:n){
  frontier3[i] <- sum(((patchy_land_frontier[,,i]+patchy_land_growth[,,i])>=1e-5) * (patchy_land_growth[,,i]>=1e-5),na.rm = T)/(sum(Landscape3>1e-5,na.rm=T))
}


png("Figs/patchy.png",width = 8,height = 6, units = "in", res = 500)
par(mfrow = c(3,3),mar = c(2,2,2,2)+.3, mgp = c(1.3, 0.5, 0))
image(patchy_land[,,100],main = "t=100")
image(patchy_land[,,300],main = "t=300")
plot(range_time3, main = "relative range/population", type = "l", xlab = "t", ylab = "relative range and population")
lines(population3, lty = 2)
abline(v = 100, lty = 3)
abline(v = 300, lty = 3)
abline(v = 500, lty = 3)
legend("bottomright",legend = c("range","population"),lty = c(1,2))

image(((patchy_land_growth[,,100]+patchy_land_frontier[,,100])>=1e-5) * (patchy_land_growth[,,100]>=1e-5),main = "frontier t=100")
image(((patchy_land_growth[,,300]+patchy_land_frontier[,,300])>=1e-5) * (patchy_land_growth[,,300]>=1e-5),main = "frontier t=300")
plot(frontier3,type = "l", xlab = "time", ylab = "relative range in frontier")


image(patchy_land[,,500],main = "t=500")
image(patchy_land[,,1000],main = "t=1000")
plot(range_time3,population3,xlab = "relative range", ylab = "relative population")
dev.off()

