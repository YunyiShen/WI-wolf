library(Matrix)
library(Rcpp)
library(RcppArmadillo)
sourceCpp("src/reaction_diffusion.cpp")

Landscape = matrix(1, 100,100)
X_init <- matrix(0,100,100)
X_init[50,1] <- .1
r <- .02
n <- 500

## uniform landscape

uniform_land <- simulate_reaction_diffusion(X_init, Landscape, r, n)
range_time <- apply(uniform_land,3,function(w,Landscape){sum(w>1e-5)/sum(Landscape>1e-5)},Landscape)
population <- apply(uniform_land,3,function(w,Landscape){sum(w)/sum(Landscape)},Landscape)


png("Figs/hardboundry.png",width = 8,height = 4, units = "in", res = 500)
par(mfrow = c(2,3),mar = c(2,2,2,2)+.3, mgp = c(1.3, 0.5, 0))
image(uniform_land[,,1],main = "t=1")
image(uniform_land[,,100],main = "t=100")
plot(range_time, type = "l", xlab = "time", ylab = "relative range and population")
lines(population, lty = 2)
legend("bottomright",legend = c("range","population"),lty = c(1,2))
image(uniform_land[,,200],main = "t=200")
image(uniform_land[,,300],main = "t=300")
plot(range_time,population,xlab = "relative range", ylab = "relative population")
dev.off()

## with gradient

gradient_1d <- c(seq(0.1,1,length.out = 20),rep(1,60),seq(1,.1,length.out = 20))
Landscape2 <- as.matrix(gradient_1d)%*%t(as.vector(gradient_1d))
image(Landscape2)

core_land <- simulate_reaction_diffusion(X_init, Landscape2, r, n)
range_time2 <- apply(core_land,3,function(w,Landscape){sum(w>1e-5)/sum(Landscape>1e-5)},Landscape2)
population2 <- apply(core_land,3,function(w,Landscape){sum(w)/sum(Landscape)},Landscape2)


png("Figs/softboundry.png",width = 8,height = 4, units = "in", res = 500)
par(mfrow = c(2,3),mar = c(2,2,2,2)+.3, mgp = c(1.3, 0.5, 0))
image(core_land[,,1],main = "t=1")
image(core_land[,,100],main = "t=100")
plot(range_time2, main = "relative range/population", type = "l", xlab = "t", ylab = "relative range and population")
lines(population2, lty = 2)
legend("bottomright",legend = c("range","population"),lty = c(1,2))
image(core_land[,,200],main = "t=200")
image(core_land[,,300],main = "t=300")
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
Landscape31 <- (Landscape31 * (Landscape31>=0.0))
Landscape31 <- apply(Landscape31,c(1,2),min,1)
image(Landscape31)
Landscape3 <- Landscape31
Landscape3[Landscape3==0] <- 0
X_init2 <- 0 * Landscape3
X_init2[1,100] <- .1
n <- 1000

patchy_land <- simulate_reaction_diffusion(X_init2, Landscape3, r, n)

range_time3 <- apply(patchy_land,3,function(w,Landscape){sum(w>1e-5)/sum(Landscape>1e-5)},Landscape3)
population3 <- apply(patchy_land,3,function(w,Landscape){sum(w)/sum(Landscape)},Landscape3)

png("Figs/patchy.png",width = 8,height = 4, units = "in", res = 500)
par(mfrow = c(2,3),mar = c(2,2,2,2)+.3, mgp = c(1.3, 0.5, 0))
image(patchy_land[,,1],main = "t=1")
image(patchy_land[,,100],main = "t=100")
plot(range_time3, main = "relative range/population", type = "l", xlab = "t", ylab = "relative range and population")
lines(population3, lty = 2)
legend("bottomright",legend = c("range","population"),lty = c(1,2))
image(patchy_land[,,300],main = "t=300")
image(patchy_land[,,1000],main = "t=1000")
plot(range_time3,population3,xlab = "relative range", ylab = "relative population")
dev.off()

