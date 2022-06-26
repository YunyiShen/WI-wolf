## This file load the data and perfomed some quick Exploratory data analysis (EDA)

rawdata <- read.csv("../data/Wolf_Time_Series.csv")
rawdata <- rawdata[rawdata$year!=2021,]
wolf <- rawdata[rawdata$X=="Wolves",]
wolf_range <- rawdata[rawdata$X=="Range",]
pack <- rawdata[rawdata$X=="Packs",]
wolf_den <- wolf
wolf_den$Winter.Minimum.Count <- wolf$Winter.Minimum.Count/wolf_range$Winter.Minimum.Count

year <- wolf$year

save.image("wolf.RData")

# plot all data
par(mfrow = c(2,2))
plot(year,wolf$Winter.Minimum.Count, ylab = "counts")
plot(year,wolf_range$Winter.Minimum.Count, ylab = "range")
plot(year,pack$Winter.Minimum.Count,ylab = "pack counts")
plot(year,wolf$Winter.Minimum.Count/wolf_range$Winter.Minimum.Count, ylab = "pack density")
dev.off()


densities <- (wolf$Winter.Minimum.Count/wolf_range$Winter.Minimum.Count)[-41] # remove 2021, which uses new method
lambdas <- wolf$Winter.Minimum.Count[-41]/wolf$Winter.Minimum.Count[-1]
plot(densities, lambdas)
abline(h = 1, col = "red")

DD_data <- data.frame(lambda = lambdas,densities = densities, year = as.factor( year[-41]) )

library(ggplot2)
ggplot(DD_data, aes(x = densities, y = lambda, col = year)) + geom_point() + geom_hline(yintercept = 1)


