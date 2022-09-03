## see if we really have expansion dominated growth 
area_pop<- data.frame(Population=wolf$Winter.Minimum.Count, Range=wolf_range$Winter.Minimum.Count, area = "WI", Year = wolf$year)
area_pop_lm <- lm(Population~Range-1, area_pop)
area_pop_pred <- predict(area_pop_lm, se = T)

png("./figs/Pop_vs_range.png", width = 6, height = 2.5, res = 500, unit = "in")

par(mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))
plot(Population~Range,area_pop)
abline(area_pop_lm)
polygon(x = c(area_pop$Range, rev(area_pop$Range)),
        y = c(area_pop_pred$fit - qt(0.975,area_pop_pred$df)*(area_pop_pred$se.fit), 
              rev(area_pop_pred$fit + qt(0.975,area_pop_pred$df)*area_pop_pred$se.fit)),
        col =  adjustcolor("black", alpha.f = 0.10), border = NA)
text(x = 1e4, y = 800, # Coordinates
     label = expression("Population = 0.0254 * Range\n p<2e-16, R^2=0.99"))
dev.off()

## direct test for stationary of the last two decades
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

## pack size 

pack_pop<- data.frame(Population=wolf$Winter.Minimum.Count, Pack=pack$Winter.Minimum.Count)
pack_pop_lm <- lm(Population~Pack-1, pack_pop)
pack_pop_pred <- predict(pack_pop_lm, se = T)

png("./figs/Pop_vs_pack.png", width = 6, height = 2.5, res = 500, unit = "in")
par(mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))
plot(Population~Pack,pack_pop, xlab = "Pack count")
abline(pack_pop_lm)
polygon(x = c(pack_pop$Pack, rev(pack_pop$Pack)),
        y = c(pack_pop_pred$fit - qt(0.975,pack_pop_pred$df)*(pack_pop_pred$se.fit), 
              rev(pack_pop_pred$fit + qt(0.975,pack_pop_pred$df)*pack_pop_pred$se.fit)),
        col =  adjustcolor("black", alpha.f = 0.10), border = NA)
text(x = 70, y = 800, # Coordinates
     label = expression("Population = 3.78 * Pack\n p<2e-16, R^2=0.99"))
dev.off()


## Also in UP of Michigan
MI_data <- read.csv("../data/MI_UP_Wolf_Density_overtime_1995_2013.csv")
plot(N~Area.occupied..km2., MI_data, xlab = "Range", ylab = "Population", ylim = c(0,800))
MI_lm <- lm(N~Area.occupied..km2., data = MI_data)
abline(MI_lm)
summary(MI_lm)

## combine dataset

MI_area_pop <- data.frame(Population = MI_data$N, Range = MI_data$Area.occupied..km2., area = "MI")
WI_MI_area_pop <- rbind(area_pop, MI_area_pop) |> na.omit()

WI_MI_area_pop <- WI_MI_area_pop[order(WI_MI_area_pop$Range),]

wimi_area_pop_lm <- lm(Population~Range-1, WI_MI_area_pop)
wimi_area_pop_pred <- predict(wimi_area_pop_lm, se = T)

png("./figs/Pop_vs_range_WI_MI.png", width = 6, height = 3.5, res = 500, unit = "in")

par(mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))
plot(Population~Range,WI_MI_area_pop)
points(Population~Range,WI_MI_area_pop[WI_MI_area_pop$area=="MI",], pch = 7)
abline(wimi_area_pop_lm)
polygon(x = c(WI_MI_area_pop$Range, rev(WI_MI_area_pop$Range)),
        y = c(wimi_area_pop_pred$fit - qt(0.975,wimi_area_pop_pred$df)*(wimi_area_pop_pred$se.fit), 
              rev(wimi_area_pop_pred$fit + qt(0.975,wimi_area_pop_pred$df)*wimi_area_pop_pred$se.fit)),
        col =  adjustcolor("black", alpha.f = 0.10), border = NA)
text(x = 1e4, y = 800, # Coordinates
     label = expression("Population = 0.0243 * Range\n p<2e-16, R^2=0.98"))
legend("bottomright", legend = c("WI","MI-UP"),pch = c(1,7))
dev.off()

wimi_area_pop_lm2 <- lm(Population~Range:area+Range, WI_MI_area_pop)

## add MI to WI
MI_data2 <- na.omit(MI_data)
mi_wi_merge <- merge(area_pop, MI_data2)
mi_wi_merge$Range <- mi_wi_merge$Range + mi_wi_merge$Area.occupied..km2.
mi_wi_merge$Population <- mi_wi_merge$Population + mi_wi_merge$N
mi_wi_merge <- mi_wi_merge[,c("Year","Range","Population")]

png("./figs/Pop_vs_range_WI_MI_added.png", width = 6, height = 3.5, res = 500, unit = "in")

par(mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))

plot(Population~Range, data = mi_wi_merge)
wimi_area_pop_merge_lm <- lm(Population~Range-1, mi_wi_merge)
wimi_area_pop_merge_lm_pred <- predict(wimi_area_pop_merge_lm, se = T)
abline(wimi_area_pop_merge_lm)
polygon(x = c(mi_wi_merge$Range[order(mi_wi_merge$Range)], rev(mi_wi_merge$Range[order(mi_wi_merge$Range)])),
        y = c(wimi_area_pop_merge_lm_pred$fit[order(mi_wi_merge$Range)] - qt(0.975,wimi_area_pop_merge_lm_pred$df)*(wimi_area_pop_merge_lm_pred$se.fit[order(mi_wi_merge$Range)]), 
              rev(wimi_area_pop_merge_lm_pred$fit[order(mi_wi_merge$Range)] + qt(0.975,wimi_area_pop_merge_lm_pred$df)*wimi_area_pop_merge_lm_pred$se.fit[order(mi_wi_merge$Range)])),
        col =  adjustcolor("black", alpha.f = 0.10), border = NA)
text(x = 2e4, y = 1200, # Coordinates
     label = expression("Population = 0.0242 * Range\n p<2e-16, R^2=0.98"))
dev.off()
