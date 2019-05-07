library(iotools)
library(backports, "/share/pkg/r/3.5.2/install/lib64/R/library")
library(devtools,"/share/pkg/r/3.5.2/install/lib64/R/library")
library(rjags,"/share/pkg/r/3.5.2/install/lib64/R/library" )
library(coda, "/share/pkg/r/3.5.2/install/lib64/R/library" )
library(ecoforecastR)


source('/usr3/graduate/tmccabe/mccabete/Fire_forecast_509/scripts/tess_geo_fork/fire_area_forecast/01.2.1_data__helper_functions.R')
load('/usr3/graduate/tmccabe/mccabete/Fire_forecast_509/output/mcmc/20190506.historical_calibration_with_modis_and_viirs_better_priors_and_sigma.convergence_passed_GBR_test.JAGS_run.Rdata')


### Make forward funciton
ensemble_number <- 700
ne <- ensemble_number
X <- matrix(NA, nrow = ne, ncol =1)


SSEM <- function(X,params,inputs, ensemble_number){ # X is previous step latent state
  mu <- matrix(NA, nrow = ensemble_number, ncol = 1)
  #X <- as.matrix(X)
  beta_precip <- params$beta_precip
  beta_temp <- params$beta_temp
  beta_sigma <- params$sigma
  
  mu <- X[,1] + beta_precip*inputs[1] + beta_temp*inputs[2]#+beta_IC
  print(mu)
  N <- rnorm(ensemble_number, mu, abs(beta_sigma))
  
 #N[is.nan(N)] <- 999999999999 #put max value, assuming that it is giving infinity
  return(N)
}

### Set parameters
tmp <- as.matrix(jags.burn)
beta_precip_orig <- tmp[,1]
beta_temp_orig <- tmp[,2]
sigma_orig <- tmp[,3]

params <- list()
params$beta_precip <-  rnorm(ne, mean(beta_precip_orig), sd(beta_precip_orig)) #sample(beta_precip_orig, ne)
params$beta_temp <-  rnorm(ne, mean(beta_temp_orig), sd(beta_temp_orig)) #sample(beta_temp_orig, ne)#
params$sigma_orig <- rnorm (ne, mean(sigma_orig, ne))# sample(sigma_orig, ne)

## Assign initial conditions
tmp <-as.matrix(jags.burn)
last_x <- 14 + 3 # 14th week is April, aproximatly when we will be starting
X.orig <- sample(tmp[, last_x], ne)
X.orig <- as.matrix(X.orig)


## Set up modis and GEFs data
modis <- read.csv("/usr3/graduate/tmccabe/mccabete/Fire_forecast_509/data/MOD14A2/2019/MOD14A2.csv")
modis_days <- get_days(first_date = "20190306", data_type = "MOD14A2")
dates <- format(as.POSIXct(modis$X1), "%Y%m%d")
modis <- modis$X2[dates %in% modis_days]

### Read in temp and precip ensembles
#met_days <- get_days(first_date = "20190306", data_type = "GEFS", day_interval = 8)
met_days <- modis_days # Issue with lubridate is breaking the helper function
temp_list <- matrix(ncol = 21, nrow = length(met_days))
precip_list <- matrix(ncol = 21, nrow = length(met_days))
nt <- length(met_days)

for (i in seq_along(met_days)){
  file_path <- paste("/usr3/graduate/tmccabe/mccabete/Fire_forecast_509/data/GEFS", met_days[i], sep = "/" )
  temp_name <- paste(met_days[i], "_full_ensemble_tempurature.csv", sep ="")
  precip_name <- paste(met_days[i], "_full_ensemble_precipitation.csv", sep ="")
  
  temp_path <- paste (file_path, temp_name, sep = "/")
  precip_path <- paste (file_path, precip_name, sep = "/")
  temp_object <- read.csv(temp_path)
  temp_object <- na.omit(temp_object[,-1])
  precip_object <- read.csv(precip_path)
  precip_object <- na.omit(precip_object[,-1])

  ### Sumamrize to 8-day agggrigate (super basic)
  for (j in 1:21){
    precip_list[i,j] <- mean(precip_object[,j])
    temp_list[i,j] <- max(temp_object[,j])
  }
  precip <- rowMeans(precip_list)
  temp <- rowMeans(temp_list)
  inputs <- cbind(precip, temp)
}

## Define storage variable for forward ensemble
output = array(0.0,c(nt,ne,1))

## Foreward ensemble simulation
nt <- length(modis_days)  #Number of steps we are forecasting
X <- X.orig
for(t in 1:nt){
  output[t,,] <- SSEM(X,params,inputs[t,], ne)
  X <- as.matrix(output[t,,1])
  print(t)
}

#output_na_omit <- output[!is.nan(output)]
## Plot Foreward ensemble sim
par(mfrow=c(1,1))
ci = apply(output[,,],1,quantile,c(0.025,0.5,0.975), na.rm = TRUE)
plot(ci[2,],main="Forward Ensemble",xlab="time",ylab="Burn Area",type='l',ylim=range(ci))
ciEnvelope(1:ncol(ci),ci[1,],ci[3,],col=col.alpha("lightGrey",0.5))
lines(ci[2,])


###### Non-resampling particle filter ######
#window <- rep(1, nt)
window <- 1:7
LAIm = t(apply(output[,,1],2,tapply,window,mean))

##temporary <- as.matrix(output[,,1])
#LAIm = tapply(output[,,1],2, INDEX = 1, FUN = mean)

#LAIm <- output
LAIlike = array(NA,dim(LAIm))
LAIm.ci  = apply(LAIm,2,quantile,c(0.025,0.5,0.975))
sel=1:ncol(LAIm.ci)
LAIr <- modis
LAIr.sd <- jitter(rep(sd(modis), 7))## MADE UP LAIr.sd FOR NOW NEED TO UPDATE WITH REAL
sel=1:ncol(LAIm.ci)
for(i in 1:ne){
  LAIlike[i,] = dnorm(LAIm[i,],LAIr[sel],LAIr.sd[sel], log = TRUE)  ## calculate log likelihoods
  print(LAIlike[i,])
  #LAIlike[i,is.na(LAIlike[i,])] = 0       ## missing data as weight 1; log(1)=0
  LAIlike[i,] <- LAIlike[i, ]
  LAIlike[i,] = exp(cumsum(LAIlike[i,]))  ## convert to cumulative likelihood (Had to remove exp due to Inf errors)
}
hist(LAIlike[,ncol(LAIlike)],main="Final Ensemble Weights")

## Non-resampling Particle Filter
## calculation of CI
nobs = ncol(LAIlike)                     ## number of observations
LAIpf = matrix(NA,3,nobs)
wbar = apply(LAIlike,2,mean)             ## mean weight at each time point
for(i in 1:nobs){
  LAIpf[,i] = wtd.quantile(LAIm[,i],LAIlike[,i]/wbar[i],c(0.025,0.5,0.975))  ## calculate weighted median and CI
}




## plot original ensemble and PF with data
Mtime <- as.Date(modis_days)
Msel <- 1:7
col.pf   = c(col.alpha("lightGrey",0.5),col.alpha("lightBlue",0.5),col.alpha("lightGreen",0.5)) ## color sequence
names.pf = c("ensemble","non-resamp PF","resamp PF")                         ## legend names
plot(Mtime[Msel],LAIm.ci[2,], ylim=range(c(range(LAIm.ci),range(LAIr,na.rm=TRUE))),
     type='n',ylab="Burn Area",xlab="Time")
ciEnvelope(Mtime,LAIm.ci[1,],LAIm.ci[3,],col=col.pf[1])                ## original ensemble
ciEnvelope(Mtime,LAIpf[1,],LAIpf[3,],col=col.pf[2])                    ## non-resampling Particle Filter
points(Mtime,LAIr)                                                           ## observations
for(i in 1:length(LAIr)){                                                    ## observation uncertainty
    lines(rep(Mtime[i],2),LAIr[i]+c(-1,1)*LAIr.sd[i])
}
legend("bottomright",legend=names.pf[1:2],col=col.pf[1:2],lwd=5)
