#============================================================================
#
#    State Space Model: Kalman Filter 
#
#============================================================================

rm (list = ls(all=TRUE))
graphics.off()

#install.packages("quantmod")
library(quantmod) # Load the package including getSymbols(), coredata(), and index()

rGDP<-getSymbols("GDPC1",src="FRED",auto.assign = FALSE) 
tail(rGDP)

RGDP<-coredata(rGDP)
y <- log(RGDP) #-mean(log(RGDP)) 
y <- y*100 #We later estimate the percentage deviation of output from its potential. 
t <- nrow(y)


##################################################################
#--------------------------- Estimation --------------------------
# Load quarterly US data for the period 1947:1 to 2023:1 (T =305)
##################################################################

#load functions for likelihood, Kalman filter, and smoothing
source("F:\\teaching\\Time_Series_Analysis_2023\\R_code_examples_2023\\kalman_filter\\likelihood_kalman.R")

# Estimate the model
initial.values <- c(1.1,-0.2, 0.5, 0.25,0.1)
output <- optim(initial.values, likelihood, y=y, method="Nelder-Mead")
#Algorithms: ¡°Nelder-Mead¡±, ¡°BFGS¡±, ¡°CG¡±, ¡°L-BFGS-B¡±, ¡°SANN¡±, ¡°Brent¡±
para <- output$par
para
likelihood.value <- output$val
likelihood.value


# Smoothing
f.smooth <- smoothing(y, para)
ss <- f.smooth$ss
trend <- ss[,1] #trend of real GDP


#--------------figure----------------------------------------
dev.new()
start<-1
nn<-nrow(y)
date<-seq(1947.0, 2023.0, by=0.25)
matplot( date,cbind(y[start:nn], ss[start:nn,1]), type='l', 
         main = '', ylab = 'Log(real GDP)*100', xlab = 'year', lwd=2)

#output gap
output.gap=y[,1]-f.smooth$ss[,1] #output gap = real GDP - potential output 

#figure for the output gap
dev.new()
start=53 
date2<-seq(1960.0, 2023.0, by=0.25)
plot(date2,output.gap[start:nn], lty=1, type="l", col="blue", ylab="output gap", xlab="", lwd=2) 
abline(h=0, lty=3, col="blue")


#growth rate and its trend
dev.new()
start=13 
date2<-seq(1950.0, 2023.0, by=0.25)
plot(date[start:nn], 4*f.smooth$ss[start:nn,4], type = "l", 
     col="blue", lwd=3, ylab = "growth rate", xlab = "year",
     ylim=c(min(4*diff(y)), max(4*diff(y))))
lines(date[start:nn], 4*diff(y[(start-1):nn]), type = "l", col="red")
text(date[length(date)-60], 12, labels="trend growth rate", col = "blue", cex=1.5)
text(date[length(date)-60], 10, labels="actual growth rate", col = "red", cex=1.5)

