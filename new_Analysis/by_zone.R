wi_zone <- read.csv("../data/WI_zone_wolfnumbers_07202022.csv")
focused_zone <- c("1","2","3","4","5","6")

for(zone in focused_zone){
  temp <- wi_zone[wi_zone$Zone==zone, ]
  
  logistic_fit <- nls(Wolves~p0*exp(r*(Year-1980))/(1+(p0*invK)*(exp(r*(Year-1980))-1)),
                      data = temp, start = list(p0=1, invK=0.01, r = 0.1))
  
  logistic_predict <- predict(logistic_fit, newdata = temp, se = T)
  
  png(paste0("./figs/population_zone_",zone,".png"), width = 4, height = 3, res = 500, unit = "in")
  par(mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))
  
  plot(temp$Y
       ,logistic_predict, type = "l", ylab = "Population", xlab = 'Year', ylim = c(0,1.1 * max(temp$Wolves)))
  
  points(temp$Y[-c(34:41)]
         ,temp$Wolves[-c(34:41)]
  )
  points(temp$Y[c(34:41)]
         ,temp$Wolves[c(34:41)]
         ,pch = 7)
  abline(v = 1980+(-sum(log(logistic_fit$m$getPars()[c("p0","invK")])))/logistic_fit$m$getPars()[c("r")], lty = 2)
  abline(h = .5/logistic_fit$m$getPars()["invK"], lty = 2)
  #legend("topleft",legend = c("observed:pre-hunting","observed:post-hunting",
                              
  #                            "logistic"), 
  #       lty = c(NA,NA,1), pch = c(1,7,NA), col = c("black","black","black"))
  dev.off()
}

area_pop<- data.frame(Population=wolf$Winter.Minimum.Count, Range=wolf_range$Winter.Minimum.Count, area = "WI")
area_pop_lm <- lm(Population~Range-1, area_pop)
area_pop_pred <- predict(area_pop_lm, se = T)


png("./figs/pop_vs_range_by_zone_all.png", width = 6, height = 3.5, res = 500, unit = "in")

par(mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))

plot(wolf_range$Winter.Minimum.Count, area_pop_pred$fit, xlab = "Range", ylab = "Population", type = "l", xlim = c(0,15000), ylim = c(0,450))
abline(area_pop_lm)

focused_zone2 <- c("1","2","3","4","5","6")
for(zone in focused_zone2){
  temp <- wi_zone[wi_zone$Zone==zone, ]
  #png(paste0("./figs/population_zone_",zone,".png"), width = 4, height = 3, res = 500, unit = "in")
  points(temp$Pack.Area.in.Zone.m2/1e6,temp$Wolves, pch = as.numeric(zone))
  
}

legend("topleft",legend = c(paste("zone",focused_zone2),
                            "overall fit"), 
       lty = c(rep(NA,length(focused_zone2)),1), pch = c(as.numeric(focused_zone2),NA),
       )
dev.off()


## an overall range fit

logistic_fit_range <- nls(Winter.Minimum.Count~p0*exp(r*(year-1980))/(1+(p0*invK)*(exp(r*(year-1980))-1)),
                    data = wolf_range, start = list(p0=1000, invK=1e-5, r = 0.1))


