library(deSolve)

toy_model <- function(t, x, parms = 10){
  A <- x[2]
  D <- x[1]
  return(list(
    c(D * (1-D) * A,
    parms * A * (1-A) * D)
  ))
}

ts <- seq(0,200,length.out = 500)
xstart <- c(.05,.05)
ode(
  func=toy_model,
  y=xstart,
  times=ts,
  parms = 1/20
) |>
  as.data.frame() -> out_0.05


ts <- seq(0,30,length.out = 500)
xstart <- c(.05,.05)
ode(
  func=toy_model,
  y=xstart,
  times=ts,
  parms = 1.
) |>
  as.data.frame() -> out_1.1

ts <- seq(0,4.5,length.out = 500)
xstart <- c(.05,.05)
ode(
  func=toy_model,
  y=xstart,
  times=ts,
  parms = 20
) |>
  as.data.frame() -> out_20


png("./figs/time_scale_ode.png", width = 5,height = 5, res = 500, units = "in")
par(mfrow = c(3,2),mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))
plot(out_0.05[,1],out_0.05[,2], xlab = "time", ylab = "rA/rD=0.05",type = "l", ylim = c(0,1))
lines(out_0.05[,1],out_0.05[,3], lty = 2)
lines(out_0.05[,1],out_0.05[,3] * out_0.05[,2], lty = 3)
legend("bottomright", legend = c("Density","Range","Population"),lty = c(1,2, 3))
plot(out_0.05[,3],out_0.05[,3] * out_0.05[,2], xlab = "range", ylab = "population")
abline(0,1)

plot(out_1.1[,1],out_1.1[,2], xlab = "time", ylab = "rA/rD=1",type = "l", ylim = c(0,1))
lines(out_1.1[,1],out_1.1[,3], lty = 2)
lines(out_1.1[,1],out_1.1[,3] * out_1.1[,2], lty = 3)
legend("topleft", legend = c("Density","Range", "Population"),lty = c(1,2, 3))
plot(out_1.1[,3],out_1.1[,3] * out_1.1[,2], xlab = "range", ylab = "population")
abline(0,.05)
abline(0,1)

plot(out_20[,1],out_20[,2], xlab = "time", ylab = "rA/rD=20",type = "l", ylim = c(0,1))
lines(out_20[,1],out_20[,3], lty = 2)
lines(out_20[,1],out_20[,3] * out_20[,2], lty = 3)
legend("topleft", legend = c("Density","Range", "Population "),lty = c(1,2,3))
plot(out_20[,3],out_20[,3] * out_20[,2], xlab = "range", ylab = "population")
abline(0,.05)
dev.off()



