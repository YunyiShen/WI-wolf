library(rdd)
source("./plot.RD") # hack the rdd package plot function to allow x and y labs

# take the last 20 years again 
rdd_wolf <- wolf[21:41,]
rdd_den <- wolf_den[21:41,]
rdd_range <- wolf_range[21:41, ]
rdd_pack <- pack[21:41,]

rdd_wolf_fit <- RDestimate(Winter.Minimum.Count~year, data = rdd_wolf, cutpoint = 2013)
plot(rdd_wolf_fit, ylab = "Population",xlab = "Year")
abline(v = 2013)
summary(rdd_wolf_fit)


rdd_range_fit <- RDestimate(Winter.Minimum.Count~year, data = rdd_range, cutpoint = 2013)
plot(rdd_range_fit, xlab = "Year",ylab = "Range")
abline(v = 2013)
summary(rdd_range_fit)

rdd_pack_fit <- RDestimate(Winter.Minimum.Count~year, data = rdd_pack, cutpoint = 2013)
plot(rdd_pack_fit, xlab = "Year",ylab = "Pack")
abline(v = 2013)
summary(rdd_range_fit)
