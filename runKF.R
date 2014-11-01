#Run Kalman Filter
#CTS 25 July 2014


####################
#Run Kalman filter

#Set working directory
setwd('C:/Users/csolom7/Documents/Research projects/C model/KF code')

#Load function
source('kalmanStep.R')
source('kalmanState.R')

#Load test data
load('C:/Users/csolom7/Documents/Research projects/C model/Data/Test/testData')

#Initial guess of parameters
log(truePars["sProcDIC"]^2)
log(truePars["sProcDOC"]^2)
parGuess <- c(0.001,8e-07,0.3,9,13) #r, iota, k, log(DIC process error variance), log(DOC process error variance)

#Initial estimate of observation error covariance matrix
H <- matrix(c(10^2,0,0,100^2),nrow=2,byrow=TRUE)

#Check that kalmanStep function works
kalmanStep(p=parGuess,Y=lakeCObs[,2:3],H=H,drivers=data1)

#Run Kalman filter
fit <- optim(parGuess,kalmanStep,control=list(maxit=1e6),Y=lakeCObs[,2:3],H=H,drivers=data1)
fit

#Pull together true, guess, and ML parameter values and compare
#All parameters in natural units (no log-transformed values), and proc errors expressed as sd
guessPars <- parGuess
guessPars[4:5] <- sqrt(exp(guessPars[4:5]))
estPars <- fit$par
estPars[4:5] <- sqrt(exp(estPars[4:5]))
parsMatrix <- rbind(truePars,guessPars,estPars)
colnames(parsMatrix) <- c("r","iota","k","sProcDIC","sProcDOC")
rownames(parsMatrix) <- c("true","guess","ML")
print(signif(parsMatrix,3))

#Run KF with ML parameter estimates to pull out estimates of state
stateEst <- kalmanState(p=fit$par,Y=lakeCObs[,2:3],H=H,drivers=data1)

#Plot these state estimates along with the simulated (and maybe observed?) data
plot(DIC~time,data=lakeCSim,type='l',xlab="DOY",ylab="Lake C (mol)")
points(stateEst[1,]~lakeCSim$time,type='l',col='blue')
par(new=T)
plot(DOC~time,data=lakeCSim,type='l',lty=2,axes=F,xlab="",ylab="")
axis(4)
points(stateEst[2,]~lakeCSim$time,type='l',lty=2,col="blue")
legend("bottomright",legend=c("DIC (left axis)","DOC (right axis)"),lty=c(1,2),bty='n')


