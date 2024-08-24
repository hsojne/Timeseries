rm(list = ls(all=T)) # This code clears all

###################### ARMA(1,1) ##########################

phi=.9
theta=.3

n=1000
y=rep(0,n)
set.seed(1234)
e=rnorm(n)

#Generate y series using the ARMA(1,1) model
for(t in 2:n) {
  y[t]=phi*y[t-1]+e[t]+theta*e[t-1]
}

#plot y against x 
x=seq(1, n, 1)
dev.new()
plot(x,y,type="l")


#################################################
# ARMA(1,1): Estimation of phi, theta, and sigma 
#################################################
fn_mle2=function(p){
  
  phi=p[1]
  theta=p[2]
  sigma=p[3]
  u=rep(0,length(y))
  
  for(t in 2:length(y)) {
    u[t]=y[t]-phi*y[t-1]-theta*u[t-1]
  }
  
  log_likihood <-0
  for (i in 2:length(y)) {
    log_likihood = log_likihood + (-(1/2)*log(2*pi)-(1/2)*log(sigma^2) -(u[i]^2)/(2*sigma^2)) 
  }
  
  minus_log_likihood = -log_likihood
  return(minus_log_likihood)
}

output=optim(par=c(.1,.1, 1),fn_mle2, method="BFGS", hessian=TRUE)
output$par

information <- output$hessian           
var_cov <- solve(information)
var_cov

se <-sqrt(diag(var_cov))  #The function diag takes the diagonal elements. 
se

t <-output$par/se
t

#print the estimates of the model parameters, standard errors, and t-values
summary.output<-data.frame(output$par, se, t)
colnames(summary.output)<-c("parameter", "standard error", "t-value")
rownames(summary.output)<-c("phi", "theta", "sigma")
summary.output


phi.hat <- output$par[1]
theta.hat <- output$par[2]
sigma.hat <- output$par[3]

#---------Do the estimated errors track the actual errors?---------
#u.hat for the estimated errors.
#e contains the actual errors. 
u.hat<-rep(0,n) 
for (t in 2:n) {
  u.hat[t]=y[t]-phi.hat*y[t-1]-theta.hat*u.hat[t-1]
}

dev.new()
plot(e, type = "l", col="blue")
lines(u.hat, type="l", col="red")

#---------fitted values-------------------------------------------------
fit<-rep(0,n)

for (t in 2:n) {
  fit[t]=phi.hat*y[t-1]+theta.hat*u.hat[t-1]
}

dev.new()
plot(y, type = "l", col="blue", ylim = c(-10, 10))
lines(fit, type="l", col="red")
legend("topleft", legend=c("actual data","fitted value"), 
       lty=c(1,1), col=c("blue", "red"), cex=1.5)


#-----estimation of the model using the built-in function, arima()------
#results <- arima(y, order=c(1,0,1)) 
#results #print results


#---------forecasting---------------------------------------------------
#install.packages("forecast")
library('forecast')
n=length(y)
output.arma = arima(y[1:(n-10)], order=c(1,0,1)) #estimation based on y[1:(n-10)]
fcast <- forecast(output.arma, h=10)

dev.new()
plot(as.numeric(fcast$mean), main="forecasting and 80% interval", type = "l", 
     col="blue", ylim = c(-3,6), ylab = "y")
lines(as.numeric(fcast$lower[,1]), type="l", col="grey") #fcast$lower[,2]) for 95% interval
lines(as.numeric(fcast$upper[,1]), type="l", col="grey")
lines(y[991:n], type="l", col="red")
legend("bottomleft", legend=c("prediction", "lower bound", "upper bound", "actual series"), 
       lty=c(1,1,1,1), col=c("blue","grey","grey", "red"), cex=1.5, box.lty=0)



#------Plot the forecasts along with the actual series from 980 to 990---------------------
# references
# http://www.alisonsinclair.ca/2011/03/shading-between-curves-in-r/
# http://www.sthda.com/english/wiki/add-legends-to-plots-in-r-software-the-easiest-way
# https://www.datanovia.com/en/blog/awesome-list-of-657-r-color-names/
# The function rev() is for reversing its arguements.

forecast<-c(y[981:990], as.numeric(fcast$mean)) #actual data + forecasts
lower.forecast<-c(y[981:990],as.numeric(fcast$lower[,1]))
upper.forecast<-c(y[981:990],as.numeric(fcast$upper[,1]))
actual.series<-y[981:n]

dev.new()
plot(981:n, forecast, main="forecasts and 80% interval", type = "l", 
     col="blue", ylim = c(-3,6), ylab = "y", xlab="period")
polygon(c(990:n,rev(990:n)),
        c(upper.forecast[10:20],rev(lower.forecast[10:20])),
        col="gray92",border="gray92")
lines(981:n,forecast, type = "l", col="blue", lwd=2)
lines(981:n,actual.series, type="l", col="red", lwd=2)
legend("bottomleft", legend=c("prediction", "Confidence Interval", "actual series"), 
       lty=c(1,1,1),lwd=c(2,6,2), col=c("blue","gray60", "red"), cex=1.5, box.lty=0)

