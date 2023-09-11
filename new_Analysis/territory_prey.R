deer_by_zone <- read.csv("../data/DeerInfo1980_2012.csv")
wolf_by_zone <- read.csv("../data/WI_zone_wolfnumbers_07202022.csv")

deer_den <- deer_by_zone[,c("Year","Zone","MeanSAK_Den")]
wolf_by_zone$meanden <- wolf_by_zone$Wolves/(wolf_by_zone$Pack.Area.in.Zone.m2/(1000*1000))

deer_wolf <- merge(deer_den, wolf_by_zone, 
                   by = c("Year","Zone"))[,c("Year","Zone","meanden","MeanSAK_Den")]
deer_wolf$MeanSAK_Den <- deer_wolf$MeanSAK_Den / 2.58999 # to square km

deer_wolf_lm <- lm(formula = meanden ~ MeanSAK_Den, data = deer_wolf)
deer_wolf_pred <- predict(deer_wolf_lm, se = T, 
                          newdata = data.frame(MeanSAK_Den = 20:130))


pdf("./figs/wolfden_vs_deer.pdf", width = 6, height = 2.5)

par(mar = c(3,4,2,1), mgp = c(1.8, 0.5, 0))
plot(meanden ~ MeanSAK_Den, deer_wolf, xlab = "Estimated deer density", ylab = "Average local density\n in zone and year")
abline(deer_wolf_lm)

polygon(x = c(20:130, rev(20:130)),
        y = c(deer_wolf_pred$fit - qt(0.975,deer_wolf_pred$df)*(deer_wolf_pred$se.fit), 
              rev(deer_wolf_pred$fit + qt(0.975,deer_wolf_pred$df)*deer_wolf_pred$se.fit)),
        col =  adjustcolor("black", alpha.f = 0.10), border = NA)
text(x = 1.1e2, y = 0.07, # Coordinates
     label = expression("p=0.106, R^2=0.011"))
dev.off()
