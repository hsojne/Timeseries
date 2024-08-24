rm(list = ls(all=T)) # This code clears all
graphics.off()

###################### Taylor Rule ##########################
install.packages("quantmod")
library(quantmod) # Load the package

#download data from the Federal Reserve Bank of St. Louis
ffr<-getSymbols("FEDFUNDS",src="FRED",from="1948-01-01", to="2023-03-01", auto.assign = FALSE)
unemployment<-getSymbols("UNRATE",src="FRED",from="1948-01-01", to="2023-03-01", auto.assign = FALSE) 
#unemployment trend
unemployment.n<-getSymbols("NROU",src="FRED",from="1948-01-01", to="2023-03-01", auto.assign = FALSE) 
pce <- getSymbols("PCEPI",src="FRED",from="1948-01-01", to="2023-03-01",auto.assign = FALSE) 
deflator <- getSymbols("GDPDEF",src="FRED",from="1948-01-01", to="2023-03-01",auto.assign = FALSE) 


date.ffr <-index(ffr) #monthly 
ffr.m <- ffr[(date.ffr >= "1960-01-01" & date.ffr <= "2008-12-01"), ]

date.un <-index(unemployment) #monthly
unemployment.m <- unemployment[(date.un >= "1960-01-01" & date.un <= "2008-12-01"), ]
tail(unemployment.m)

date.n <-index(unemployment.n) #quarterly
unemployment.n.q <- unemployment.n[(date.n >= "1960-01-01" & date.n <= "2008-12-01"), ]
tail(unemployment.n.q)


date.pce <- index(pce) #monthly
inflation.pce <- 1200*diff(log(pce)) #annualize the data
inflation.pce.m <- inflation.pce[(date.pce >= "1960-01-01" & date.pce <= "2008-12-01"), ] 
head(inflation.pce.m)

date.deflator <- index(deflator) #quarterly
inflation.def <- 400*diff(log(deflator)) #annualize the data
inflation.def.q <- inflation.def[(date.deflator >= "1960-01-01" & date.deflator <= "2008-12-01"), ]
tail(inflation.def.q)


#----------Convert monthly data into quarterly data------------------------
n <- length(ffr.m)
i<-1
j<-1
ffr.q <- rep(0,n/3)
unemployment.q <- rep(0,n/3)
inflation.pce.q <- rep(0,n/3)

while (i <= n) {
  ffr.q[j] <-mean(ffr.m[i:(i+2),1])
  inflation.pce.q[j] <-mean(inflation.pce.m[i:(i+2),1])
  unemployment.q[j] <-mean(unemployment.m[i:(i+2),1])
  i <- i+3
  j <- j+1
}

#unemployment gap
unemploymentgap.q<-unemployment.q-unemployment.n.q #unemployment gap
tail(unemploymentgap.q)

#--------------------plot the data-----------------------------------------
dev.new()
par(mfrow=c(2,2))#2by1 graph
date2 <- seq(1960.0, 2008.75, 0.25)
plot(date2, ffr.q, type="l", ylab="Federal Funds Rate", xlab="date", col="blue")
plot(date2, unemploymentgap.q, type="l",ylab="Unemployment Rate", xlab="date", col="blue")
plot(date2, inflation.pce.q, type="l",ylab="PCE Inflation", xlab="date", col="blue")
plot(date2, inflation.def.q, type="l",ylab="GDP Deflator Inflation", xlab="date", col="blue")
par(mfrow=c(1,1))#2by1 graph


###########################################################################
# Taylor Rule: Estimation 
###########################################################################
fn_mle=function(p){
  
  rho=p[1]
  phi=p[2]
  theta=p[3]
  sigma=p[4]
  
  u=rep(0,length(ffr.q))
  
  for(t in 2:length(ffr.q)) {
    u[t]=ffr.q[t]-rho*ffr.q[t-1]-(1-rho)*(phi*inflation.def.q[t]+theta*unemploymentgap.q[t])
  }
  
  log_likihood <-0
  for (i in 2:length(ffr.q)) {
    log_likihood = log_likihood + (-(1/2)*log(2*pi)-(1/2)*log(sigma^2) -(u[i]^2)/(2*sigma^2)) 
  }
  minus_log_likihood = -log_likihood
  return(minus_log_likihood)
}

#------------------estimation and t-values------------------------------------
output=optim(par=c(.8,1.5,-0.5, 1),fn_mle, method="BFGS", hessian=TRUE)

parameters<-output$par
parameters

information <- output$hessian           
var_cov <- solve(information) #variance-covariance matrix for the estimates
var_cov

se <-sqrt(diag(var_cov))  #The function diag takes the diagnal elements. 
se

t <-output$par/se
t

#print the estimates of the model parameters, standard errors, and t-values
summary.output<-data.frame(output$par, se, t)
colnames(summary.output)<-c("parameter", "standard error", "t-value")
rownames(summary.output)<-c("rho", "phi", "theta", "sigma")
summary.output

#------------------fitted values and actual FFR-------------------------------
rho <- output$par[1]
phi <- output$par[2]
theta <- output$par[3]
n<-length(ffr.q)
fit<-rep(0,n)
for (t in 2:n) {
  fit[t] <- rho*ffr.q[t-1]+(1-rho)*(phi*inflation.def.q[t]+theta*unemploymentgap.q[t])
}

#------------------plot fitted values and FFR-----------------------------------------
dev.new()
date3 <- seq(1960.25, 2008.75, 0.25)
plot(date3, ffr.q[2:n], type="l", lty=1, lwd=2, col="blue", xlab="date", ylab="")
lines(date3, fit[2:n], type="l", lty=4, lwd=2, col="red")
legend("topleft", legend=c("FFR","fitted values"),col=c("blue", "red"), lty=1:2, cex=1.0, box.lty=0)

#residual
u.hat<-ffr.q[2:n]-fit[2:n]
plot(date3, u.hat, type="l", col="blue")