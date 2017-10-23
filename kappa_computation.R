library(LambertW)

# evaluate the psi function

psi = function(x) {
  sqrt(-suppressWarnings(W(-x^2, branch = -1)))
}

# calculate optimal band with known kappa

opt.ci <- function(time, kappa){
  c_hat = psi(kappa%*%t(time))*sqrt(t(replicate(nrow(kappa), c(time))))
  return(c_hat)
}

args<-commandArgs(TRUE)

# simulate coverage table

coveragetable = function(L, U, outfile) {
  K = 4*10^5
  N = 10^5
  dt = seq(1/K, 1, by = 1/K)

  kappa = t(t(seq(L, U, by = 0.0005)))
  time = t(t(seq(0.0025, 1, by = 0.0025)))

  lower = -opt.ci(dt, kappa)
  upper = -lower

  kappa_time = matrix(0, nrow = length(kappa), ncol = nrow(time)-1)

  lower_tltu = vector("list", length(kappa))
  upper_tltu = vector("list", length(kappa))
  
  m = length(time)

  for (k in 1:length(kappa)){
    lower_tltu[[k]] = mapply(function(x, y) lower[k, (dt > x & dt < y)], time[1:(m-1)], time[2:m])
    upper_tltu[[k]] = mapply(function(x, y) upper[k, (dt > x & dt < y)], time[1:(m-1)], time[2:m])
    print(k)
  }

  for (n in 1:N){
    W = rnorm(K, 0, 1/sqrt(K)); W = cumsum(W)
    W_tltu = mapply(function(x, y) W[(dt > x & dt < y)], time[1:(m-1)], time[2:m])
  
    for (k in 1:length(kappa)){
      intervalprod = mapply(function(x, y, z) prod((x >= y)*(x <= z)), W_tltu, lower_tltu[[k]], upper_tltu[[k]])
      kappa_time[k,] = kappa_time[k,] + rev(cumprod(rev(intervalprod)))
    }
    print(n)
  }

  crit.val = data.frame(kappa_time/N)
  rownames(crit.val) = seq(L, U, by = 0.0005)
  colnames(crit.val) = seq(0.0025, 0.9975, by = 0.0025)
  write.csv(crit.val, paste(getwd(),"/", outfile, sep =""))
}

L = as.numeric(args[1])
U = as.numeric(args[2])
outfile = args[3]
coveragetable(L,U,outfile)

