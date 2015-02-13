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

#Keep only every 5th observation
keepRows <- seq(1,dim(lakeCObs)[1],5)
lakeCObs[-keepRows,2:3] <- NA

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

#Plot observed data, estimated states, and rough 95% CI for predicted state (+- 2 sd)
par(mfrow=c(2,1),mar=c(3.1,3.1,1.1,1.1))
plot(DIC~time,data=lakeCObs,xlab="DOY",ylab="Lake C (mol)")
points(stateEst$a[1,]~lakeCObs$time,type='l')
sd <- sqrt(stateEst$P[1,1,])
points(lakeCObs$time,(stateEst$a[1,]+2*sd),type='l',col='gray')
points(lakeCObs$time,(stateEst$a[1,]-2*sd),type='l',col='gray')
plot(DOC~time,data=lakeCObs,xlab="DOY",ylab="Lake C (mol)")
points(stateEst$a[2,]~lakeCObs$time,type='l')
sd <- sqrt(stateEst$P[2,2,])
points(lakeCObs$time,(stateEst$a[2,]+2*sd),type='l',col='gray')
points(lakeCObs$time,(stateEst$a[2,]-2*sd),type='l',col='gray')


#Plot state estimates along with the simulated data
plot(DIC~time,data=lakeCSim,type='l',xlab="DOY",ylab="Lake C (mol)")
points(stateEst[1,]~lakeCSim$time,type='l',col='blue')
par(new=T)
plot(DOC~time,data=lakeCSim,type='l',lty=2,axes=F,xlab="",ylab="")
axis(4)
points(stateEst[2,]~lakeCSim$time,type='l',lty=2,col="blue")
legend("bottomright",legend=c("DIC (left axis)","DOC (right axis)"),lty=c(1,2),bty='n')


