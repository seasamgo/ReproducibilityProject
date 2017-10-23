require(optband)

gen.data <- function(n, mu, nu=NA, strata=F){
  
  ## generate surival data for simulation
  ## n = sample size, mu = mean failure rate
  ## censoring rate dependent on mu::nu
  ## strata=T indicates 2-sample data
  
  x <- rexp(n, rate=mu)
  if(is.na(nu)==T) ## no censoring
    return(as.data.frame(cbind(x,numeric(n)+1)))
  else {
    x2 <- rexp(n, rate=nu)
    d <- x < x2; x <- pmin(x,x2)
    if(strata==T){
      s <- rbinom(n,1,0.5)
      return(as.data.frame(cbind(x,d,s)))
    }
    return(as.data.frame(cbind(x,d)))
  }
}

band.sim <- function(n, R, mu, nu, method, fun, restr, c.l){
  
  ## band area sim
  ## R = iteration #, method = opt. method, fun = estimated function
  ## restr = truncation, c.l = coverage level
  
  area <- vector(length = R)
  for(r in 1:R){
    
    ## generate data and fit
    dat <- gen.data(n=n, mu=mu, nu=nu)
    S_KM <- survfit(Surv(dat[[1]], dat[[2]]) ~ 1, type="kaplan-meier")
    L = quantile(S_KM$time[S_KM$n.event>0], restr[1])
    U = quantile(S_KM$time[S_KM$n.event>0], restr[2])
    
    if(method == "hall-wellner" || method == "epband"){
      fit <- km.ci(S_KM, conf.level=c.l, method = method, tl = L, tu = U)
      time <- fit$time; time_shift <- shift(time); time_shift[length(time)] <- time_shift[length(time)-1]
      area[r] = riemsum(time, fit$upper - fit$lower)
      
    } else if(method == "optband"){
      fit <- opt.ci(S_KM, conf.level=c.l, fun = fun, tl = L, tu = U)
      time <- fit$time; time_shift <- shift(time); time_shift[length(time)] <- time_shift[length(time)-1]
      area[r] = riemsum(time, fit$upper - fit$lower)
      
    } else
      stop("Method not recognized")
  }
  return(mean(area))
}

coverage.sim <- function(n, R, mu, nu, method, fun, restr, c.l, t=NA){
  
  ## band coverage sim
  ## R = iteration #, method = opt. method, fun = estimated function
  ## restr = truncation, c.l = coverage level
  
  prob = vector(length = R)
  
  for (r in 1:R){
    
    ## generate data and fit
    
    dat <- gen.data(n=n, mu=mu, nu=nu)
    if(is.na(t)!=T){
      dat[dat[,1] >= t,1]=t
      dat[dat[,1] >= t,2]=0
    }
    S_KM <- survfit(Surv(dat[[1]], dat[[2]]) ~ 1, type = "kaplan-meier")
    L <- quantile(S_KM$time[S_KM$n.event>0], restr[1])
    U <- quantile(S_KM$time[S_KM$n.event>0], restr[2])
    
    if(method == "hall-wellner" || method == "epband"){
      fit <- km.ci(S_KM, conf.level=c.l, method = method, tl = L, tu = U)
      time <- fit$time; time_shift <- shift(time); time_shift[length(time)] <- time_shift[length(time)-1]
      success <- cbind(exp(-mu*time_shift)<=fit$upper, exp(-mu*time)>=fit$lower)
      prob[r] <- prod(success[,1]*success[,2], na.rm = T)*n
      
    } else if(method == "optband"){
      fit <- opt.ci(S_KM, conf.level=c.l, fun = "surv", tl = L, tu = U)
      time <- fit$time; time_shift <- shift(time); time_shift[length(time)] <- time_shift[length(time)-1]
      success <- cbind(exp(-mu*time_shift)<=fit$upper, exp(-mu*time)>=fit$lower)
      prob[r] <- prod(success[,1]*success[,2], na.rm = T)*n
      
    } else
      stop("Method not recognized")
    
  }
  return(sum(prob)/(R*n))
}

compare.sim <- function(n, R, mu, nu, fun, restr, c.l){
  
  ## comparison simulations based on band area and coverage
  ## for a particular N, R, censor %, coverage and truncation
  
  a <- c(); c <- c(); methods=c("hall-wellner", "epband", "optband")
  
  ## band area and coverage sims
  for(i in 1:3){
    a <- c(a, band.sim(n=n, R=R, mu=mu, nu=nu,  method=methods[i], fun=fun, restr=restr, c.l=c.l))
    c <- c(c, coverage.sim(n=n, R=R, mu=mu, nu=nu,  method=methods[i], fun=fun, restr=restr, c.l=c.l))
  }
  
  ## calculate area relative to Hall-Wellner
  tab <- cbind(c,a/a[1])
  return(tab)
}

table.sim <- function(n, R=10000, mu=1, nu, fun="surv", truncation=NA, cov.lev=NA){
  
  ## table of simulations based on band area and coverage
  ## for a particular N, R, censor %
  ## coverage at (90%, 95%, 99%) and truncation at various levels
  
  if(is.na(truncation)==T)
    truncation <- rbind(c(0,.999999), c(.01,.99), c(.2,.8), c(.01,0.8), c(.2,.99))
  if(is.na(cov.lev)==T)
    cov.lev <- c(.9,.95,.99)
  tab <- c()
  
  for(i in 1:nrow(truncation)){
    section <- c()
    for(j in 1:length(cov.lev))
      section <- rbind(section, compare.sim(n=n, R=R,  mu=mu, nu=nu, fun=fun, restr=truncation[i,], c.l=cov.lev[j]))
    tab <- cbind(tab,section)
  }
  return(round(tab,2))
}

plot.survival <- function(survi, methods="optband", fun="surv",
                          bands=T, legend=T, color=NULL, position="topright", yl=c(0,1), main="", samples=1){
  ## Plots
  ## survi is a survival object or data frame, methods may be one of (optband, epband, hall-wellner)
  ## fun indicates survival or cumulative-hazard function, bands/legend=T indicates plot bands/legend
  ## optional color, position, main arguments alter defaults for plot
  ## yl determines ylim in plot, samples indicates 1 or 2 sample cases
  
  ## Create surv object if data input
  if (data.class(survi) != "survfit")
    survi <- survfit(Surv(survi[[1]], survi[[2]])~1, type="kaplan-meier")
  
  L = quantile(survi$time[survi$n.event>0], .01)
  U = quantile(survi$time[survi$n.event>0], .99)
  
  if(is.null(color))
    color <- c("grey", "darkblue", "red", "green")
  
  par(pty="s")
  
  m=length(methods)
  if(bands==T & fun=="surv"){
    
    plot(survi, mark.time=FALSE, xlab="",
         ylab="", main=main, col = "white", ylim=yl)
    
    for(i in 1:m){
      
      method=methods[i]
      
      if(method == "optband"){
        c <- opt.ci(survi, fun=fun, conf.level = 0.95)
        lines(c, col=color[i+1], lty=2, lwd=2, mark.time=F)
      } else if(method == "hall-wellner" || method == "epband"){
        c <- km.ci(survi, conf.level = 0.95, method = method)
        lines(c, col=color[i+1], lty=3, lwd=2, mark.time=F)
      } else {
        stop("Method not recognized")
      }
    }
    lines(survi,col="grey", lwd=2, mark.time=F)
    
  } else if (length(methods)>1 || methods != "optband")
    stop("All methods available for survival curve, but only optband for cumulative hazard")
  else if (bands==T & fun == "cumhaz" & samples==1){
    position="topleft"
    plot(survi, fun=fun, mark.time=FALSE, xlab="",
         ylab="", main=main, col = color[1], pch=20)
    c <- opt.ci(survi, fun=fun, conf.level = 0.95)
    lines(-log(c$upper)~c$time, col=color[2], lty=3, lwd=2)
    lines(-log(c$lower)~c$time, col=color[2], lty=3, lwd=2)
  } else if(bands==T & fun == "cumhaz" & samples==2){
    position="topleft"
    c <- opt.ci(survi, fun=fun, conf.level = 0.95, samples=2)
    plot(c$difference~c$time, xlab="",
         ylab="", main=main, col = "white", ylim=yl, pch=20)
    lines(c$difference~c$time, col = "grey", lty=1, lwd=2)
    lines(-log(c$upper)~c$time, col=color[2], lty=3, lwd=2)
    lines(-log(c$lower)~c$time, col=color[2], lty=3, lwd=2)
  }
  if(legend==TRUE)
    legend(position, c(fun, methods), lwd=2, lty=c(1,rep(2,m)), col=color)
  
}
