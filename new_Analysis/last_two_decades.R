# this file analyzed the last decades which seems stationary, we determine the effect of harvest
last_decade <- data.frame(year = 2000:2020, density = (wolf$Winter.Minimum.Count/wolf_range$Winter.Minimum.Count)[year>=2000])
lmtrend <- lm(density~year,last_decade[-c(13:21),]) # 2012-14 is legal harvest
lmtrend <- lm(density~year,last_decade) # 2012-14 is legal harvest
tseries::kpss.test(last_decade$density) # kpss test for stationary, seems stationary
lm_time <- lm(density~year, last_decade)
pred_lm <- predict(lm_time,se = T)
png("./figs/den_last20.png", width = 6, height = 3.5, res = 500, unit = "in")
par(mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))
plot(last_decade, xlab = "Year", ylab = "Density")
abline(lm(density~year, last_decade))
polygon(x = c(last_decade$year, rev(last_decade$year)),
        y = c(pred_lm$fit - qt(0.975,pred_lm$df)*pred_lm$se.fit, 
              rev(pred_lm$fit + qt(0.975,pred_lm$df)*pred_lm$se.fit)),
        col =  adjustcolor("black", alpha.f = 0.10), border = NA)
text(x = 2015, y = 0.0285, # Coordinates
     label = expression("Density = 0.212-0.000097 * Year\n p=0.13, R^2=0.067"))
dev.off()
