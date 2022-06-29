# this file analyzed the last decades which seems stationary, we determine the effect of harvest
last_decade <- data.frame(year = 2000:2020-2010, density = (wolf$Winter.Minimum.Count/wolf_range$Winter.Minimum.Count)[year>=2000])
lmtrend <- lm(density~year,last_decade[-c(13:21),]) # 2012-14 is legal harvest
lmtrend <- lm(density~year,last_decade) # 2012-14 is legal harvest
tseries::kpss.test(last_decade$density) # kpss test for stationary, seems stationary
lm_time <- lm(density~year, last_decade)
pred_lm <- predict(lm_time,se = T)
png("./figs/den_last20.png", width = 6, height = 3.5, res = 500, unit = "in")
par(mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))
plot(2000:2020, last_decade$density, xlab = "Year", ylab = "Density")
lines(2000:2020,pred_lm$fit)
polygon(x = c(last_decade$year+2010, rev(last_decade$year+2010)),
        y = c(pred_lm$fit - qt(0.975,pred_lm$df)*pred_lm$se.fit, 
              rev(pred_lm$fit + qt(0.975,pred_lm$df)*pred_lm$se.fit)),
        col =  adjustcolor("black", alpha.f = 0.10), border = NA)
text(x = 2014, y = 0.0285, # Coordinates
     label = expression("Density = 0.0259-9.7e-5 * (Year-2010)\n p=0.13, R^2=0.067"))
dev.off()
