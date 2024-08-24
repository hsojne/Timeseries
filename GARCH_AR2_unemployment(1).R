rm(list = ls(all=T)) # This code clears all
graphics.off()
################################################
library(quantmod) # Load the package

unemployment<-getSymbols("UNRATE",src="FRED",auto.assign = FALSE) 
date <-index(unemployment)
unemployment <- unemployment[(date >= "1960-01-01" & date <= "2019-12-01"), ] #fix the data set 

head(unemployment)
tail(unemployment)

n<-length(unemployment)
n

#Eliminate dates from the data
date <-index(unemployment)
unemployment<-coredata(unemployment)
head(unemployment)


dev.new()
plot(date, unemployment, type="l",ylab="Unemployment Rate", xlab="date")


#-------------AR(2) model with constant variance---------------------------
#sigma가 시점별로 다른 경우
results <- arima(unemployment, order=c(2,0,0))  
results #print results

sigma2_ar2<-results$sigma2
coefficient<-results$coef
coefficient

var<-results$var.coef #variance-covariance matrix
var

t_ar2<-coefficient/sqrt(diag(var))
t_ar2 #t values

###########################################################################
# Unemployment: Estimation of an AR(2) model with time-varying variance
###########################################################################
fn_mle_g=function(p){
  alpha<-abs(p)
  n <- length(unemployment)

  sigma2 <- sigma2_ar2*rep(1, n) #sigma2_ar is the variance from the estimated AR(2) model above.
                                 #We will use this value as a starting point. 
  u=rep(0,n)
  
  #u=error
  for(t in 3:n) {
    u[t]=unemployment[t]-alpha[1]-alpha[2]*unemployment[t-1]-alpha[3]*unemployment[t-2]
  }
  

  for (i in 3:n) {
    sigma2[i] <- alpha[4] + alpha[5]*u[i-1]^2 + alpha[6]*sigma2[i-1] 
  }
  
  
  log_likihood <-0
  for (i in 3:n) {
    log_likihood = log_likihood + (-(1/2)*log(2*pi)-(1/2)*log(sigma2[i]) -(u[i]^2)/(2*sigma2[i])) 
  }
  
  minus_log_likihood=-log_likihood
  return(minus_log_likihood)
}

#------------------estimation and t-values------------------------------------
output=optim(par=c(0.5, 0.9, 0.13, 0.21, 0.2, 0.6),fn_mle_g, method="BFGS", hessian=TRUE)
output                                 #Algorithms: ??Nelder-Mead??, ??BFGS??, ??CG??, ??L-BFGS-B??, ??SANN??, ??Brent??

output$par

information <- output$hessian           
var_cov <- solve(information)
var_cov

se <-sqrt(diag(var_cov))  #The function diag takes the diagnal elements. 
se

t <-output$par/se
t

#------------------fitted values and actual unemployment--------------------------------------
phi0 <- output$par[1]
phi1 <- output$par[2]
phi2 <- output$par[3]
theta <- output$par[4:6]

fit<-rep(0,n)
for (t in 3:n) {
  fit[t] <-phi0+phi1*unemployment[t-1]+phi2*unemployment[t-2]
}

#------------------plot fitted values and unemployment-----------------------------------------
dev.new()
date2<-date[3:n]
plot(date2, unemployment[3:n], type="l", lty=1, lwd=1, col="blue", xlab="date", ylab="")
lines(date2, fit[3:n], type="l", lty="dotted", lwd=2, col="red")
legend("topleft", legend=c("unemployment","fitted values"),col=c("blue", "red"), lty=1:2, cex=1.0)

#------------------plot errors and volatility (variance)-----------------------------------------  

# Compute conditional variance at optimal parameters
u=rep(0,n)

for(t in 3:n) {
  u[t]=unemployment[t]-phi0-phi1*unemployment[t-1]-phi2*unemployment[t-2] #u=unemployment[3:n]-fit[3:n]
}

sigma2 <- sigma2_ar2*rep(1,n)

for (i in 3:n) {
  sigma2[i] <- theta[1] + theta[2]*u[i-1]^2 + theta[3]*sigma2[i-1]
}

#----------------graph of conditional variance-----------------------------
dev.new()  
par(mar=c(5.0,5.0,1.0,1.0), mfcol=c(2,1)) #sets the bottom, left, top and right margins respectively of the plot region in number of lines of text. 
n1=1
plot(date[n1:n], sigma2[n1:n], type="l", xlab="date", ylab=expression(sigma[t]^2), col="blue")
plot(date[n1:n], u[n1:n], type="l", xlab="date", ylab="error", col="blue")
par(mfcol=c(1,1))


