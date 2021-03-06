#### load libraries ####
library(nimble, warn.conflicts = FALSE, quietly = TRUE)
library(dplyr)
library(mcmcplots)
library(coda)

#### Data Preparation ####
df  <- read.csv('prior_data.csv',header=TRUE)
N   <- dim(df)[1]

prevData    <- list(Xi = df$reactive)
prevConst   <- list(N = N,ni = df$tested)

#### starting values ####
prevInits   <- list(Se = 0.88, Sp = 0.98,
                    Tpi = rep(0.5, prevConst$N))

########## Weakly Informative - Primary Analysis ##########

# Derive alpha and beta
nWI       <- 10
alphaWI   <- rep(0, N)
betaWI    <- rep(0, N)
x         <- df$Super.region3
for (i in 1:N){
  alphaWI[i]  <- x[i]/100 + 1
  betaWI[i]   <- nWI - x[i]/100 + 1
}
alphaWI
betaWI

muWI <- rep(0,N)
varWI <- rep(0,N)

# Derive the mean and variance
for (i in 1:N){
  muWI[i]  <- alphaWI[i]/(alphaWI[i]+betaWI[i])
  varWI[i] <- alphaWI[i] * betaWI[i] / ((alphaWI[i]+betaWI[i])^2 * (alphaWI[i]+betaWI[i]+1))
}

# Model
prevWI  <- nimbleCode({
  for (i in 1:N){
    Opi[i] <- Tpi[i]*Se + (1 - Tpi[i])*(1-Sp)
    Xi[i] ~ dbin(Opi[i], ni[i])
  }
  Tpi[1]  ~ dbeta(1.1930,10.8070)
  Tpi[2]  ~ dbeta(1.1930,10.8070)
  Tpi[3]  ~ dbeta(1.1930,10.8070)
  Tpi[4]  ~ dbeta(2.6641,9.3359)
  Tpi[5]  ~ dbeta(1.1930,10.8070)
  Tpi[6]  ~ dbeta(2.6641,9.3359)
  Tpi[7]  ~ dbeta(1.1930,10.8070)
  Tpi[8]  ~ dbeta(2.1792,9.8208)
  Tpi[9]  ~ dbeta(1.1930,10.8070)
  Tpi[10] ~ dbeta(1.1930,10.8070)
  Tpi[11] ~ dbeta(1.1930,10.8070)
  Se ~ T(dbeta(205, 29),0.8,0.95) # truncated beta to plausible range
  Sp ~ T(dbeta(288, 2),0.9,1) # truncated beta to plausible range
})

burnin    <- 100000
modelWI   <- nimbleMCMC(code = prevWI,
                        constants = prevConst,
                        data = prevData,
                        inits = prevInits,
                        monitors = c("Se", "Sp",
                                     paste('Tpi[',1:N,']', sep = ''),
                                     paste('Opi[',1:N,']', sep = '')),
                        nburnin = burnin, niter = 2*burnin, nchains = 4)

### (1) True Sero-prevalence: 
## a. Convergence Diagnostics.
# Visual inspection
region <- df$Region.
par(mfrow = c(3,4))
for(i in 1:N){
  Tp1   <- modelWI$chain1[,paste('Tpi[',i,']', sep = '')] 
  Tp2   <- modelWI$chain2[,paste('Tpi[',i,']', sep = '')]
  Tp3   <- modelWI$chain3[,paste('Tpi[',i,']', sep = '')] 
  Tp4   <- modelWI$chain4[,paste('Tpi[',i,']', sep = '')]
  plot(1:burnin, Tp1, type = 'l', xlab = 'b', ylab = region[i],
       xlim = c(0, 4*burnin))
  lines(burnin + 1:burnin, Tp2, col = 2)
  lines(2*burnin + 1:burnin, Tp3, col = 3)
  lines(3*burnin + 1:burnin, Tp4, col = 4)
}

# Gelman-Rubin diagnostic
TpWI <- vector(mode="list", length=11)
GRWI <- vector(mode="list", length=11)
for (i in 1:N){
  TpWI[[i]]   <- c(modelWI$chain1[,paste('Tpi[',i,']', sep = '')], 
                   modelWI$chain2[,paste('Tpi[',i,']', sep = '')],
                   modelWI$chain3[,paste('Tpi[',i,']', sep = '')],
                   modelWI$chain4[,paste('Tpi[',i,']', sep = '')])
  GRWI[[i]]   <- mcmc.list(as.mcmc(modelWI$chain1[,paste('Tpi[',i,']', sep = '')]),
                           as.mcmc(modelWI$chain2[,paste('Tpi[',i,']', sep = '')]), 
                           as.mcmc(modelWI$chain3[,paste('Tpi[',i,']', sep = '')]), 
                           as.mcmc(modelWI$chain4[,paste('Tpi[',i,']', sep = '')]))
  print(gelman.diag(GRWI[[i]]))
}

## b.Credible Intervals
Tp_q_WI <- vector(mode="list", length=N)

for (i in 1:N){
  Tp_q_WI[[i]]  <- quantile(TpWI[[i]], probs = c(0.5, 0.025, 0.975))
}

tp_q_WI             <- data.frame(matrix(unlist(Tp_q_WI), nrow=N, byrow=T), 
                                stringsAsFactors = FALSE)
tp_q_WI$Regions     <- df$Region.
new_tp_q_WI         <- tp_q_WI %>%
  select(Regions, everything())
names(new_tp_q_WI)  <- c("Regions","50%", "2.5%", "97.5%")
new_tp_q_WI


### (2) Observed Prevalence

## Convergence Diagnostics
#  Visual inspection
par(mfrow = c(3,4))
for(i in 1:N){
  Op1   <- modelWI$chain1[,paste('Opi[',i,']', sep = '')] 
  Op2   <- modelWI$chain2[,paste('Opi[',i,']', sep = '')]
  Op3   <- modelWI$chain3[,paste('Opi[',i,']', sep = '')] 
  Op4   <- modelWI$chain4[,paste('Opi[',i,']', sep = '')]
  plot(1:burnin, Op1, type = 'l', xlab = 'b', ylab = region[i],
       xlim = c(0, 4*burnin))
  lines(burnin + 1:burnin, Op2, col = 2)
  lines(2*burnin + 1:burnin, Op3, col = 3)
  lines(3*burnin + 1:burnin, Op4, col = 4)
}

# Gelman-Rubin diagnostic
GRWI_Op <- vector(mode="list", length=11)
for (i in 1:N){
  GRWI_Op[[i]]  <- mcmc.list(as.mcmc(modelWI$chain1[,paste('Opi[',i,']', sep = '')]), 
                             as.mcmc(modelWI$chain2[,paste('Opi[',i,']', sep = '')]), 
                             as.mcmc(modelWI$chain3[,paste('Opi[',i,']', sep = '')]), 
                             as.mcmc(modelWI$chain4[,paste('Opi[',i,']', sep = '')]))
  print(gelman.diag(GRWI_Op[[i]]))
}


### (3) Sensitivity and Specificity
SeWI  <- c(modelWI$chain1[,'Se'], modelWI$chain2[,'Se'],
           modelWI$chain3[,'Se'], modelWI$chain4[,'Se'])
SpWI  <- c(modelWI$chain1[,'Sp'], modelWI$chain2[,'Sp'],
           modelWI$chain3[,'Sp'], modelWI$chain4[,'Sp'])

par(mfrow = c(1,2))
plot(SeWI, type = 'l')
plot(SpWI, type = 'l')

# Gelman-Rubin diagnostic
GRWI_Se   <- mcmc.list(as.mcmc(modelWI$chain1[,'Se']), as.mcmc(modelWI$chain2[,'Se']), 
                       as.mcmc(modelWI$chain3[,'Se']), as.mcmc(modelWI$chain4[,'Se']))
print(gelman.diag(GRWI_Se))

GRWI_Sp   <- mcmc.list(as.mcmc(modelWI$chain1[,'Sp']), as.mcmc(modelWI$chain2[,'Sp']), 
                       as.mcmc(modelWI$chain3[,'Sp']), as.mcmc(modelWI$chain4[,'Sp']))
print(gelman.diag(GRWI_Sp))

quantile(SeWI, probs = c(0.5, 0.025, 0.975))
quantile(SpWI, probs = c(0.5, 0.025, 0.975))

##########  Non-Informative Model  ##########
prevFlat  <- nimbleCode({
  for (i in 1:N){
    Tpi[i]  ~ dbeta(1, 1)
    Opi[i]  <- Tpi[i]*Se + (1 - Tpi[i])*(1-Sp)
    Xi[i]   ~ dbin(Opi[i], ni[i])
  }
  
  Se ~ T(dbeta(205, 29),0.8,0.95) # truncated betas to plausible range
  Sp ~ T(dbeta(288, 2),0.9,1) # truncated betas to plausible range
})

burnin  <- 100000
modelF  <- nimbleMCMC(code = prevFlat,
                      constants = prevConst,
                      data = prevData,
                      inits = prevInits,
                      monitors = c("Se", "Sp",
                                   paste('Tpi[',1:N,']', sep = ''),
                                   paste('Opi[',1:N,']', sep = '')),
                      nburnin = burnin, niter = 2*burnin, nchains = 4)

### (1) True Sero-prevalence:

## a. Convergence Diagnostics.
# Visual inspection
par(mfrow = c(3,4))
for(i in 1:N){
  Tp1   <- modelF$chain1[,paste('Tpi[',i,']', sep = '')] 
  Tp2   <- modelF$chain2[,paste('Tpi[',i,']', sep = '')]
  Tp3   <- modelF$chain3[,paste('Tpi[',i,']', sep = '')] 
  Tp4   <- modelF$chain4[,paste('Tpi[',i,']', sep = '')]
  plot(1:burnin, Tp1, type = 'l', xlab = 'b', ylab = region[i],
       xlim = c(0, 4*burnin))
  lines(burnin + 1:burnin, Tp2, col = 2)
  lines(2*burnin + 1:burnin, Tp3, col = 3)
  lines(3*burnin + 1:burnin, Tp4, col = 4)
}

# Gelman-Rubin diagnostic
TpF <- vector(mode="list", length=11)
GRF <- vector(mode="list", length=11)
for (i in 1:N){
  TpF[[i]]  <- c(modelF$chain1[,paste('Tpi[',i,']', sep = '')], 
                 modelF$chain2[,paste('Tpi[',i,']', sep = '')],
                 modelF$chain3[,paste('Tpi[',i,']', sep = '')],
                 modelF$chain4[,paste('Tpi[',i,']', sep = '')])
  GRF[[i]]  <- mcmc.list(as.mcmc(modelF$chain1[,paste('Tpi[',i,']', sep = '')]), 
                         as.mcmc(modelF$chain2[,paste('Tpi[',i,']', sep = '')]), 
                         as.mcmc(modelF$chain3[,paste('Tpi[',i,']', sep = '')]), 
                         as.mcmc(modelF$chain4[,paste('Tpi[',i,']', sep = '')]))
  print(gelman.diag(GRF[[i]]))
}



## b.Credible Intervals for Each Region
Tp_q_F <- vector(mode="list", length=N)

for (i in 1:N){
  Tp_q_F[[i]]   <- quantile(TpF[[i]], probs = c(0.5, 0.025, 0.975))
}

tp_q_F            <- data.frame(matrix(unlist(Tp_q_F), nrow=11, byrow=T), stringsAsFactors = FALSE)
tp_q_F$Regions    <- df$Region.
new_tp_q_F        <- tp_q_F %>% 
  select(Regions, everything())
names(new_tp_q_F) <- c("Regions","50%", "2.5%", "97.5%")
new_tp_q_F

### (2) Observed Sero-prevalence

## Convergence Diagnostics
#  Visual inspection
par(mfrow = c(3,4))
for(i in 1:N){
  Op1   <- modelF$chain1[,paste('Opi[',i,']', sep = '')] 
  Op2   <- modelF$chain2[,paste('Opi[',i,']', sep = '')]
  Op3   <- modelF$chain3[,paste('Opi[',i,']', sep = '')] 
  Op4   <- modelF$chain4[,paste('Opi[',i,']', sep = '')]
  plot(1:burnin, Op1, type = 'l', xlab = 'b', ylab = region[i],
       xlim = c(0, 4*burnin))
  lines(burnin + 1:burnin, Op2, col = 2)
  lines(2*burnin + 1:burnin, Op3, col = 3)
  lines(3*burnin + 1:burnin, Op4, col = 4)
}

# Gelman-Rubin diagnostic
GRF_Op <- vector(mode="list", length=11)
for (i in 1:N){
  GRF_Op[[i]] <- mcmc.list(as.mcmc(modelF$chain1[,paste('Opi[',i,']', sep = '')]), 
                        as.mcmc(modelF$chain2[,paste('Opi[',i,']', sep = '')]), 
                        as.mcmc(modelF$chain3[,paste('Opi[',i,']', sep = '')]), 
                        as.mcmc(modelF$chain4[,paste('Opi[',i,']', sep = '')]))
  print(gelman.diag(GRF_Op[[i]]))
}


### (3) Sensitivity and Specificity

SeF  <- c(modelF$chain1[,'Se'], modelF$chain2[,'Se'],
          modelF$chain3[,'Se'], modelF$chain4[,'Se'])
SpF  <- c(modelF$chain1[,'Sp'], modelF$chain2[,'Sp'],
          modelF$chain3[,'Sp'], modelF$chain4[,'Sp'])

## Convergence Diagnostics
#  Visual inspection
par(mfrow = c(1,2))
plot(SeF, type = 'l')
plot(SpF, type = 'l')

# Gelman-Rubin diagnostic
GRF_Se  <- mcmc.list(as.mcmc(modelF$chain1[,'Se']), as.mcmc(modelF$chain2[,'Se']), 
                     as.mcmc(modelF$chain3[,'Se']), as.mcmc(modelF$chain4[,'Se']))
print(gelman.diag(GRF_Se))

GRF_Sp  <- mcmc.list(as.mcmc(modelF$chain1[,'Sp']), as.mcmc(modelF$chain2[,'Sp']), 
                     as.mcmc(modelF$chain3[,'Sp']), as.mcmc(modelF$chain4[,'Sp']))
print(gelman.diag(GRF_Sp))

## Credible Interval
quantile(SeF, probs = c(0.5, 0.025, 0.975))
quantile(SpF, probs = c(0.5, 0.025, 0.975))

########## Informative Model ##########

# Derive alpha and beta
nI      <- 100
alphaI  <- rep(0, N)
betaI   <- rep(0, N)
for (i in 1:N){
  alphaI[i]   <- x[i]/10 + 1
  betaI[i]    <- nI - x[i]/10 + 1
}
alphaI
betaI

muI   <- rep(0,N)
varI  <- rep(0,N)
# Derive the mean
for (i in 1:N){
  muI[i]    <- alphaI[i]/(alphaI[i]+betaI[i])
  varI[i]   <- alphaI[i] * betaI[i] / ((alphaI[i]+betaI[i])^2 * (alphaI[i]+betaI[i]+1))
}

muI   <- data.frame(muI)
varI  <- data.frame(varI)


#Model
prevI  <- nimbleCode({
  for (i in 1:N){
    Opi[i] <- Tpi[i]*Se + (1 - Tpi[i])*(1-Sp)
    Xi[i] ~ dbin(Opi[i], ni[i])
  }
  Tpi[1]  ~ dbeta(2.930,99.070)
  Tpi[2]  ~ dbeta(2.930,99.070)
  Tpi[3]  ~ dbeta(2.930,99.070)
  Tpi[4]  ~ dbeta(17.641,84.359)
  Tpi[5]  ~ dbeta(2.930,99.070)
  Tpi[6]  ~ dbeta(17.641,84.359)
  Tpi[7]  ~ dbeta(2.930,99.070)
  Tpi[8]  ~ dbeta(12.792,89.208)
  Tpi[9]  ~ dbeta(2.930,99.070)
  Tpi[10] ~ dbeta(2.930,99.070)
  Tpi[11] ~ dbeta(2.930,99.070)
  Se ~ T(dbeta(205, 29),0.8,0.95) # truncated betas to plausible range
  Sp ~ T(dbeta(288, 2),0.9,1) # truncated betas to plausible range
})

burnin  <- 100000
modelI  <- nimbleMCMC(code = prevI,
                      constants = prevConst,
                      data = prevData,
                      inits = prevInits,
                      monitors = c("Se", "Sp",
                                   paste('Tpi[',1:N,']', sep = ''),
                                   paste('Opi[',1:N,']', sep = '')),
                      nburnin = burnin, niter = 2*burnin, nchains = 4)

### (1) True Sero-prevalence: 
## a. Convergence Diagnostics.
# Visual inspection
par(mfrow = c(3,4))
for(i in 1:N){
  Tp1   <- modelI$chain1[,paste('Tpi[',i,']', sep = '')] 
  Tp2   <- modelI$chain2[,paste('Tpi[',i,']', sep = '')]
  Tp3   <- modelI$chain3[,paste('Tpi[',i,']', sep = '')] 
  Tp4   <- modelI$chain4[,paste('Tpi[',i,']', sep = '')]
  plot(1:burnin, Tp1, type = 'l', xlab = 'b', ylab = region[i],
       xlim = c(0, 4*burnin))
  lines(burnin + 1:burnin, Tp2, col = 2)
  lines(2*burnin + 1:burnin, Tp3, col = 3)
  lines(3*burnin + 1:burnin, Tp4, col = 4)
}

# Gelman-Rubin diagnostic
TpI <- vector(mode="list", length=11)
GRI <- vector(mode="list", length=11)
for (i in 1:N){
  TpI[[i]]  <- c(modelI$chain1[,paste('Tpi[',i,']', sep = '')], 
                 modelI$chain2[,paste('Tpi[',i,']', sep = '')],
                 modelI$chain3[,paste('Tpi[',i,']', sep = '')],
                 modelI$chain4[,paste('Tpi[',i,']', sep = '')])
  GRI[[i]]  <- mcmc.list(as.mcmc(modelI$chain1[,paste('Tpi[',i,']', sep = '')]), 
                         as.mcmc(modelI$chain2[,paste('Tpi[',i,']', sep = '')]), 
                         as.mcmc(modelI$chain3[,paste('Tpi[',i,']', sep = '')]), 
                         as.mcmc(modelI$chain4[,paste('Tpi[',i,']', sep = '')]))
  print(gelman.diag(GRI[[i]]))
}

## b.Credible Intervals
Tp_q_I  <- vector(mode="list", length=N)

for (i in 1:N){
  Tp_q_I[[i]]   <- quantile(TpI[[i]], probs = c(0.5, 0.025, 0.975))
}

tp_q_I              <- data.frame(matrix(unlist(Tp_q_I), nrow=N, byrow=T), stringsAsFactors = FALSE)
tp_q_I$Regions      <- df$Region.
new_tp_q_I          <- tp_q_I %>%
  select(Regions, everything())
names(new_tp_q_I)   <- c("Regions","50%", "2.5%", "97.5%")
new_tp_q_I


### (2) Observed Prevalence

## Convergence Diagnostics
#  Visual inspection - traceplot
par(mfrow = c(3,4))
for(i in 1:N){
  Op1   <- modelI$chain1[,paste('Opi[',i,']', sep = '')] 
  Op2   <- modelI$chain2[,paste('Opi[',i,']', sep = '')]
  Op3   <- modelI$chain3[,paste('Opi[',i,']', sep = '')] 
  Op4   <- modelI$chain4[,paste('Opi[',i,']', sep = '')]
  plot(1:burnin, Op1, type = 'l', xlab = 'b', ylab = region[i],
       xlim = c(0, 4*burnin))
  lines(burnin + 1:burnin, Op2, col = 2)
  lines(2*burnin + 1:burnin, Op3, col = 3)
  lines(3*burnin + 1:burnin, Op4, col = 4)
}

# Gelman-Rubin diagnostic
GRI_Op <- vector(mode="list", length=11)
for (i in 1:N){
  GRI_Op[[i]] <- mcmc.list(as.mcmc(modelI$chain1[,paste('Opi[',i,']', sep = '')]), 
                           as.mcmc(modelI$chain2[,paste('Opi[',i,']', sep = '')]), 
                           as.mcmc(modelI$chain3[,paste('Opi[',i,']', sep = '')]), 
                           as.mcmc(modelI$chain4[,paste('Opi[',i,']', sep = '')]))
  print(gelman.diag(GRI_Op[[i]]))
}

### (3) Sensitivity and Specificity
SeI  <- c(modelI$chain1[,'Se'], modelI$chain2[,'Se'],
          modelI$chain3[,'Se'], modelI$chain4[,'Se'])
SpI  <- c(modelI$chain1[,'Sp'], modelI$chain2[,'Sp'],
          modelI$chain3[,'Sp'], modelI$chain4[,'Sp'])

## Convergence Diagnostic
#  Visual inspection
par(mfrow = c(1,2))
plot(SeI, type = 'l')
plot(SpI, type = 'l')

# Gelman-Rubin diagnostic
GRI_Se  <- mcmc.list(as.mcmc(modelI$chain1[,'Se']), as.mcmc(modelI$chain2[,'Se']), 
                     as.mcmc(modelI$chain3[,'Se']), as.mcmc(modelI$chain4[,'Se']))
print(gelman.diag(GRF_Se))

GRI_Sp  <- mcmc.list(as.mcmc(modelI$chain1[,'Sp']), as.mcmc(modelI$chain2[,'Sp']), 
                     as.mcmc(modelI$chain3[,'Sp']), as.mcmc(modelI$chain4[,'Sp']))
print(gelman.diag(GRI_Sp))

quantile(SeI, probs = c(0.5, 0.025, 0.975))
quantile(SpI, probs = c(0.5, 0.025, 0.975))

