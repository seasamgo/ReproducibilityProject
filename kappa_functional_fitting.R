require(reshape2)
require(akima)
require(ks)
require(fields)
require(lattice)

#Reshaping the data to feed into heat map function interp()
critval = read.csv(paste(getwd(),"critval.csv", sep=""))
critval2 = critval
colnames(critval2) = c("kappa", seq(0.0050, 0.9950, by = 0.0050))

critval3 = melt(critval2, id = "kappa")
colnames(critval3) = c("kappa","time","gamma")
critval3$time = as.numeric(as.character(critval3$time))

#Plot contours
contour(a, pch = 22, lty = 2, lwd = 2, labcex = 1.4, levels = seq(0.90, 0.99, by = 0.01), 
        add=FALSE, xlab = "L", ylab = expression(kappa))

#Functional Fitting
fit = lm(gamma - 1 ~ -1 + I(kappa + kappa^2) + time:kappa, data = critval3)
summary(fit)

kappa.fun = function(L, conf.level) {
  a = -0.4272; b = 0.2848
  kappa = (-(a+b*L)-sqrt((a+b*L)^2-4*a*(1-conf.level)))/(2*a)
  return(kappa)
}

#Plot contours from functional fit
L.grid = seq(0, 1, by = 0.001)
kappa90 = kappa.fun(seq(0, 1, by = 0.001), 0.90)
kappa91 = kappa.fun(seq(0, 1, by = 0.001), 0.91)
kappa92 = kappa.fun(seq(0, 1, by = 0.001), 0.92)
kappa93 = kappa.fun(seq(0, 1, by = 0.001), 0.93)
kappa94 = kappa.fun(seq(0, 1, by = 0.001), 0.94)
kappa95 = kappa.fun(seq(0, 1, by = 0.001), 0.95)
kappa96 = kappa.fun(seq(0, 1, by = 0.001), 0.96)
kappa97 = kappa.fun(seq(0, 1, by = 0.001), 0.97)
kappa98 = kappa.fun(seq(0, 1, by = 0.001), 0.98)
kappa99 = kappa.fun(seq(0, 1, by = 0.001), 0.99)
lines(L.grid, kappa90, col = "blue", lwd = 2)
lines(L.grid, kappa91, col = "blue", lwd = 2)
lines(L.grid, kappa92, col = "blue", lwd = 2)
lines(L.grid, kappa93, col = "blue", lwd = 2)
lines(L.grid, kappa94, col = "blue", lwd = 2)
lines(L.grid, kappa95, col = "blue", lwd = 2)
lines(L.grid, kappa96, col = "blue", lwd = 2)
lines(L.grid, kappa97, col = "blue", lwd = 2)
lines(L.grid, kappa98, col = "blue", lwd = 2)
lines(L.grid, kappa99, col = "blue", lwd = 2)
