poland <- read.csv("../data/Poland_wolf_Nowak_and_Mysłajek_2016.csv")
plot(poland$area_min, poland$pop_size_min, xlab = "Range", ylab = "Population")
pop_range_lm <- lm((pop_size_max+pop_size_min)/2~area_min-1, data = poland)
abline(pop_range_lm)
text(x = 2e3, y = 100, # Coordinates
     label = expression("Population = 0.0124 * Range\n p<2e-16, R^2=0.99"))


par(mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))
plot(poland$year, poland$pop_size_max, xlab = "Year", ylab = "Population")
exponential_model <- nls(pop_size_max~p0*exp(r*(year-2000)),
    data = poland, start = list(p0=10, r = 0.1))
pred_exp <- predict(exponential_model, newdata = poland)
lines(poland$year, pred_exp)
quadratic_model <- lm(pop_size_max~year + I(year^2), data = poland)
pred_quad <- predict(quadratic_model, newdata = poland)
lines(poland$year, pred_quad, col = "red")
legend("topleft",legend = c("Observed","exponential","quadratic"), 
       lty = c(NA,1,1), pch = c(1,NA,NA),col = c("black","black","red")
)

AIC(exponential_model)
AIC(quadratic_model)
