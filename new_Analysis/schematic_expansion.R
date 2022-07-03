area_pop<- data.frame(Population=wolf$Winter.Minimum.Count, Range=wolf_range$Winter.Minimum.Count)
area_pop_lm <- lm(Population~Range-1, area_pop)
## draw a schematic for different growth when there are expansion

# I need a y=x^2 (resp. y=sqrt(x)) type of function pass (20000, 507.6) and y'(20000)=0.02538 
super_linear <- function(x,a){
  a*(x-20000+0.02538/(2*a))^2+507.6-a*(0.02538/(2*a))^2
}

sub_linear <- function(x,a){
  b <- 20000-(a/(2*0.02538))^2
  c <- 507.6-a^2/(2*0.02538)
  a * sqrt(x-b) + c
}

x_prime <- 20000:50000

png("./figs/Pop_vs_range_schematic.png", width = 6*1.2, height = 3.5*1.2, res = 500, unit = "in")

par(mar = c(3,3,2,2), mgp = c(1.8, 0.5, 0))
plot(Population~Range,area_pop)
abline(area_pop_lm)
lines(x_prime, super_linear(x_prime, 5e-6), lty = 2)
lines(x_prime, sub_linear(x_prime, 1), lty = 4)
lines(rep(20000,100),seq(507.6, 1500, length.out = 100), lty = 5 )
lines(x_prime, 507.6 + 0 * x_prime, lty = 3)

legend("topleft",legend = c("expand+density unchange","expand+density increase",
                            "expand+density decrease",
                            "density increase", "expand","WI-wolf"), 
       lty = c(1,2,4,5,3,NA), pch = c(NA,NA,NA,NA,NA,1))
dev.off()
