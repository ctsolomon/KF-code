#Test KF with missing obs in state variables
#CTS 12 Feb 2015

####################
#Run Kalman filter

#Set working directory
setwd('C:/Users/csolom7/Documents/Research projects/C model/KF code')

#Load function
source('kalmanStep.R')
source('kalmanState.R')

#Load test data
load('C:/Users/csolom7/Documents/Research projects/C model/Data/Test/testData')

#Introduce some NAs to the observed state variable time series
lakeCObs[10,"DIC"] <- NA
lakeCObs[20,c("DIC","DOC")] <- NA

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
