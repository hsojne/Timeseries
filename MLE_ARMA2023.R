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

#par=c(), include parameters
output=optim(par=c(.1,.1, 1),fn_mle2, method="BFGS", hessian=TRUE)
output

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
lines(fit, type="l", col="orangered")
legend("topleft", legend=c("actual data","fitted value"), 
       lty=c(1,1), col=c("blue", "red"), cex=1.5)


#-----estimation of the model using the built-in function, arima()------
#results <- arima(y, order=c(1,0,1)) 
#results #print results

