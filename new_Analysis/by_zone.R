source("interesting_time_points.R")
wi_zone <- read.csv("../data/WI_zone_wolfnumbers_07202022.csv")
habitat_zone <- read.csv("../data/ZoneSummaries_HabitatSuitability.csv")
wi_zone$Pack.Area.in.Zone.m2 <- wi_zone$Pack.Area.in.Zone.m2 * (1e-6)
focused_zone <- c("1","2","3","4","5","6")

second_order_time <- c()
third_order_time <- c()

for(zone in focused_zone){
  temp <- wi_zone[wi_zone$Zone==zone, ]
  habitat <- habitat_zone[habitat_zone$Geo==zone,c("Belant2022HabitatReclass_Area","Mladenoff09_HabitatProbability_Area")]
  logistic_fit <- nls(Pack.Area.in.Zone.m2~p0*exp(r*(Year-1980))/(1+(p0*invK)*(exp(r*(Year-1980))-1)),
                      data = temp, start = list(p0=100, invK=1e-5, r = 0.1))
  
  #logistic_fit <- nls(Wolves~p0*exp(r*(Year-1980))/(1+(p0*invK)*(exp(r*(Year-1980))-1)),
  #                    data = temp, start = list(p0=1, invK=0.01, r = 0.1))
                      
  
  logistic_predict <- predict(logistic_fit, newdata = temp, se = T)
  p0 <- logistic_fit$m$getPars()["p0"]
  invK <- logistic_fit$m$getPars()["invK"]
  r <- logistic_fit$m$getPars()["r"]
  
  
  pdf(paste0("./figs/range_zone_",zone,".pdf"), width = 4, height = 3)
  par(mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))
  
  plot(temp$Y
       ,logistic_predict, type = "l", ylab = "Range", xlab = 'Year', ylim = c(0,1.1 * max(1/logistic_fit$m$getPars()["invK"],max(temp$Pack.Area.in.Zone.m2),as.numeric(habitat))))
  
  points(temp$Y[-c(34:41)]
         ,temp$Pack.Area.in.Zone.m2[-c(34:41)]
  )
  points(temp$Y[c(34:41)]
         ,temp$Pack.Area.in.Zone.m2[c(34:41)]
         ,pch = 7)
  abline(v = 1980+(-sum(log(logistic_fit$m$getPars()[c("p0","invK")])))/logistic_fit$m$getPars()[c("r")], lty = 2, col = "#ffaf0f")
  second_order_time[zone] <- 1980+(-sum(log(logistic_fit$m$getPars()[c("p0","invK")])))/logistic_fit$m$getPars()[c("r")]
  #abline(h = .5/logistic_fit$m$getPars()["invK"], lty = 2)
  abline(v = (root_3rd_der(r,p0,invK))[1], lty = 2, col = "#55ff00")
  abline(v = (root_3rd_der(r,p0,invK))[2], lty = 2, col = "#ff0000")
  third_order_time[zone] <- max(root_3rd_der(r,p0,invK))
  abline(h = 1/logistic_fit$m$getPars()["invK"], lty = 3)
  abline(h=habitat[1],lty = 3)
  abline(h=habitat[2],lty = 3)
  text(1990,habitat[1],"Gantchoff et al. 2022")
  text(1990,habitat[2],"Mladenoff et al. 2009")
  #text(1980+(-sum(log(logistic_fit$m$getPars()[c("p0","invK")])))/logistic_fit$m$getPars()[c("r")],1/logistic_fit$m$getPars()["invK"],"'K'")
  text(max(root_3rd_der(r,p0,invK)),1/logistic_fit$m$getPars()["invK"],"'K'" )
  #legend("topleft",legend = c("observed:pre-hunting","observed:post-hunting",
                              
  #                            "logistic"), 
  #       lty = c(NA,NA,1), pch = c(1,7,NA), col = c("black","black","black"))
  dev.off()
}

area_pop<- data.frame(Population=wolf$Winter.Minimum.Count, Range=wolf_range$Winter.Minimum.Count, area = "WI")
area_pop_lm <- lm(Population~Range-1, area_pop)
area_pop_pred <- predict(area_pop_lm, se = T)


pdf("./figs/pop_vs_range_by_zone_all.pdf", width = 6, height = 3.5)

par(mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))

plot(wolf_range$Winter.Minimum.Count, area_pop_pred$fit, xlab = "Range", ylab = "Population", type = "l", xlim = c(0,15000), ylim = c(0,450))
abline(area_pop_lm)

focused_zone2 <- c("1","2","3","4","5","6")
for(zone in focused_zone2){
  temp <- wi_zone[wi_zone$Zone==zone, ]
  #png(paste0("./figs/population_zone_",zone,".png"), width = 4, height = 3, res = 500, unit = "in")
  points(temp$Pack.Area.in.Zone.m2,temp$Wolves, pch = as.numeric(zone))
  
}

legend("topleft",legend = c(paste("zone",focused_zone2),
                            "overall fit"), 
       lty = c(rep(NA,length(focused_zone2)),1), pch = c(as.numeric(focused_zone2),NA),
       bg = "white"
       )
dev.off()


## an overall range fit

logistic_fit_range <- nls(Winter.Minimum.Count~p0*exp(r*(year-1980))/(1+(p0*invK)*(exp(r*(year-1980))-1)),
                    data = wolf_range[-c(34:41),], start = list(p0=1000, invK=1e-5, r = 0.1))

second_order_time["WI"] <- 1980+(-sum(log(logistic_fit_range$m$getPars()[c("p0","invK")])))/logistic_fit_range$m$getPars()[c("r")]
p0 <- logistic_fit_range$m$getPars()["p0"]
invK <- logistic_fit_range$m$getPars()["invK"]
r <- logistic_fit_range$m$getPars()["r"]
third_order_time["WI"] <- max(root_3rd_der(r,p0,invK))

second_order_time
third_order_time
times_interesting <- data.frame(zone = c(1:6,"WI"), second = second_order_time, third = third_order_time)
write.csv(times_interesting,"./figs/interesting_time_for_range.csv", row.names = F)


# draw the mortality risk fig
mor_risk <- read.csv("../data/mortality_risk.csv")
wlm_mor <- lm(Mean~Year, data = mor_risk, weights = (1/mor_risk$Standard.deviation)^2)


ReLU <- function(x) x * (x>=0)
png(paste0("./figs/mortality.png"), width = 4, height = 3, res = 500, unit = "in")
par(mar = c(3,4,2,2), mgp = c(1.8, 0.5, 0))
plot(mor_risk$Year, (mor_risk$Mean), xlab = "Year", ylab = "Mortality risk in wolf range \n (Stenglein et al. 2018)", pch = 15, ylim = c(-.0,max(mor_risk$Mean+mor_risk$Standard.deviation)))
lines(mor_risk$Year, (mor_risk$Mean))
arrows(x0=mor_risk$Year, y0=(mor_risk$Mean-mor_risk$Standard.deviation), x1=mor_risk$Year, y1=mor_risk$Mean+mor_risk$Standard.deviation, code=3, angle=90, length=0.05, col="black", lwd=1)
abline(wlm_mor, lwd = 1.5)

dev.off()
