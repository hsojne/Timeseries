rm(list = ls(all=T)) # This code clears all
graphics.off()


################################################
install.packages("quantmod")
library(quantmod) # Load the package

unemployment<-getSymbols("UNRATE",src="FRED",from="1948-01-01", to="2024-02-01",auto.assign = FALSE) 
date <-index(unemployment)

unemployment <- unemployment[(date >= "1960-01-01" & date <= "2024-02-01"), ] #fix the data set 
date <-index(unemployment)

head(unemployment)
tail(unemployment)

n<-length(unemployment)
n

dev.new()
plot(date, unemployment, type="l",ylab="Unemployment Rate", xlab="date")

#Eliminate dates from the data
unemployment<-coredata(unemployment)
head(unemployment)
class(unemployment)

###########################################################################
# Unemployment: Estimation of an AR(2) model 
###########################################################################
fn_mle=function(p){
  
  phi0=p[1]
  phi1=p[2]
  phi2=p[3]
  sigma=p[4]
  
  u=rep(0,n)

  for(t in 3:n) {
    u[t]=unemployment[t]-phi0-phi1*unemployment[t-1]-phi2*unemployment[t-2]
  }
  
  log_likihood <-0
  for (i in 3:n) {
    log_likihood = log_likihood + (-(1/2)*log(2*pi)-(1/2)*log(sigma^2) -(u[i]^2)/(2*sigma^2)) 
  }
  minus_log_likihood <- -log_likihood
  
  return(minus_log_likihood)
}

#------------------estimation and t-values------------------------------------
output=optim(par=c(0.1,0.7,0.2,1.0),fn_mle, method="BFGS", hessian=TRUE)
output

likelihood <- output$val
likelihood <- -likelihood
likelihood 

information <- output$hessian #information matrix = -hessian           
var_cov <- solve(information) #variance-covariance matrix
var_cov

se <-sqrt(diag(var_cov))  #The function diag takes the diagnal elements. 
se

t <-output$par/se
t

summary.ar2<-data.frame(output$par, se, t)
colnames(summary.ar2)<-c("parameter", "standard error", "t-value")
rownames(summary.ar2)<-c("phi0", "phi1", "phi2", "sigma")
summary.ar2

#mean estimate of unemploymnet
mu_u<-output$par[1]/(1-output$par[2]-output$par[3])
mu_u

#mean
mean(unemployment)

# Estimate the model using the built-in function, arima. 
results <- arima(unemployment, order=c(2,0,0))
results #print results
#The function arma gives us the intercept, which is the mean value of unemployment (mu_u). 

#------------------fitted values and actual unemployment--------------------------------------
phi0 <- output$par[1]
phi1 <- output$par[2]
phi2 <- output$par[3]

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

#------------------plot errors------------------------------------------------------------------  
dev.new()
error<-unemployment[3:n]-fit[3:n]
plot(date2, error, type="l", lty=1, lwd=1, col="blue", xlab="date", ylab="Error")
