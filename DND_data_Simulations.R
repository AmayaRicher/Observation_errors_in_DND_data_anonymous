###No FPE with detection heterogeneity###

rm(list = ls())	# Clean the workspace


#install packages
library(dplyr)
library(ggplot2)
library(tidyverse)
library(igraph)
library(R6)
library(coda)
library(readxl)
library(readr)
library(numDeriv)
library(pracma)
library(posterior)
library(mcmcplots)
library(rjags)
library(bayestestR)
library(knitr)
library(kableExtra)


#To upload and print figures without running the simulations, go to end of script and import pre-saved tables

true_psi.list <- c(0.1, 0.3, 0.5)	# occupancy probability 
true_p11.list <- c(0.7, 0.9) #c(0.5, 0.7, 0.9)	# detection probability
true_p10.list <- c(0, 0.01, 0.02, 0.05)	# false-positive probability


### CREATE TABLE TO STORE SIMULATION RESULTS
TABLE <- matrix(NA, length(true_psi.list)*length(true_p11.list)*length(true_p10.list), 3 + 3 + 1 + 2 + 6*2) # 3 + 3 + 1 + 2 + 6*2 because 3 for psi, p11 and p10 + 3 for rel_psi, rel_p11, rel_p10 (realisated parameters in data set); 1 for convergence  and + 2 for estim_sd_lpsi and estim_sd_lp11;  and 2 parameters to estimate here (psi and p11) and 6 for biais, error, BCI_coverage, BCI_sidth, mse and rmse
colnames(TABLE) <- c("true_psi", "true_p11", "true_p10", "rel_psi", "rel_p11", "rel_p10", "converged", "estim_sd_lpsi", "estim_sd_lp11", "biais.psi", "error.psi", "BCI_cover.psi", "BCI_width.psi", "mse.psi", "rmse.psi", "biais.p11", "error.p11", "BCI_cover.p11", "BCI_width.p11", "mse.p11", "rmse.p11")




### FUNCTION TO SIMULATE DATA
sim.data <- function(psi = true_psi, psi_sd = true_psi/3, p11 = true_p11, p11_sd = true_p11/10, p10 = true_p10, nquadrats = 20, ncells = 100, k = 3) { #k is number of replicates per site
  lpsi <- rnorm(nquadrats, mean = psi, sd = psi_sd) #list of psi for each of the 20 quadrats
  lp11 <- rnorm(nquadrats, mean = p11, sd = p11_sd) #list of p11 for each of the 20 quadrats
  for (a in 1:nquadrats) {
    if (lpsi[a] < 0) {lpsi[a] <- 0}
    if (lpsi[a]>1) {lpsi[a] <- 1}
    if (lp11[a] < 0) {lp11[a] <- 0}
    if (lp11[a]>1) {lp11[a] <- 1}
  }
  z <- array(0,dim=c(nquadrats,ncells)) #presence (1) or absence (0)
  mu <- array(0,dim=c(nquadrats,ncells)) #'binomial' parameter, conditional on site occupancy status
  y <- array(0,dim=c(nquadrats,ncells)) #detection (1) non-detection (0)
  #z <- rbinom(nsites,1,psi)
  for (i in 1:nquadrats) {
    z[i,] <- rbinom(ncells,1,lpsi[i])
    for (j in 1:ncells) {
      mu[i,j] <- (1-z[i,j])*p10 + z[i,j]*lp11[i] #so we suppose here that p10 is constant across all quadrats and cells, but p11 varies between quadrats
      y[i,j] <- rbinom(1,k,mu[i,j]) #size = k for all 3 visits
    }
  }
  return(list(z=z, mu=mu, y=y, k=k, nquadrats=nquadrats, ncells=ncells))
}





### STAN MODELS

#NO FALSE POSITIVES MODEL
modelString = "
model{
  for(i in 1:nquadrats){  # Loop through quadrats
    # Likelihood
    for(j in 1:ncells) { # Loop through cells
      z[i, j] ~ dbern(psi[i])
      # Observation model
      mu[i,j] <- z[i, j]*p11[i]
      y[i, j] ~ dbin(mu[i,j], k)
    }
    
    # Priors for Psi (quadrat level)
    lpsi[i] ~ dnorm(mu.lpsi, tau.lpsi)
    psi[i] <- ilogit(lpsi[i])
    
    #Detection heterogeneity:
    # Priors for p11 (quadrat level)
    lp11[i] ~ dnorm(mu.lp11, tau.lp11)
    p11[i] <- ilogit(lp11[i])
    
  }
  
  # Hyperpriors for psi (mean across all quadrats)
  psi.mean ~ dunif(0.01,1) #dbeta(1, 1)
  mu.lpsi <- logit(psi.mean)
  sd.lpsi ~ dunif(0, 5)
  tau.lpsi <- 1/sd.lpsi^2
  
  # Hyperpriors for p11 (mean across all quadrats)
  p11.mean ~ dunif(0.01,1) #dbeta(1, 1)
  mu.lp11 <- logit(p11.mean)
  sd.lp11 ~ dunif(0, 5)
  tau.lp11 <- 1/sd.lp11^2
}
"

writeLines(modelString,con="model_No_FP_het.txt") 



reps=100 #number of simulations

### CREATE TABLE TO SAVE ALL of the reps (here 100) ESTIMATIONS OF PSI FOR EACH SCENARIO:
LIST_ESTIM_PSI <- matrix(NA,length(true_psi.list)*length(true_p11.list)*length(true_p10.list), 3 + reps)
colnames(LIST_ESTIM_PSI) <- c(c("true_psi", "true_p11", "true_p10"), as.character(c(1:reps)))


###------------------------------------------------------------------------###

ct <- 0 #counter to keep track of iterations (combination of parameters)
time.start <- Sys.time

### START ITERATIONS

#loop across all parameters:
for (a in 1:length(true_psi.list)) {
  for (b in 1:length(true_p11.list)) {
    for (c in 1:length(true_p10.list)) {
      
      
      ct <- ct + 1
      
      
      true_psi <- true_psi.list[a]
      true_p11 <- true_p11.list[b]
      true_p10 <- true_p10.list[c]
      
      
      ### CREATE OBJECTS TO SAVE RESULTS OF ESTIMATORS
      estim_psi <- NA
      estim_p11 <- NA 
      estim_sd_lpsi <- NA
      
      
      ### START SIMULATIONS
      
      if (is.na(TABLE[ct,'true_psi'])) {
        print(c('go for',true_psi,true_p11,true_p10))
        
        estim_psi_biais <- c() #list of psi biais (mean(psi) - true_psi) for all sites for one reps
        mean_estim_psi <- c()
        estim_psi_sd <- c()
        estim_lpsi_sd <- c()
        estim_p11_biais <- c() #same for p11
        estim_p11_sd <- c()
        estim_lp11_sd <- c()
        conv <- 0 #counts number of models that converged among the 100 reps
        BCI_cover.psi <- 0 #counts number of models for which true_psi is in BCI 95% among converged models
        BCI_width.psi <- c() #list of BCI width of all reps for psi 
        BCI_cover.p11 <- 0 #same for true_p11
        BCI_width.p11 <- c() #list of BCI width of all reps for p11
        
        
        for(r in 1:reps) {
          
          #SIMULATE THE DATA:
          x <- sim.data()
          data <- list(y=x$y, nquadrats=x$nquadrats, ncells = x$ncells, k=x$k)
          
          #obtain rel_psi, rel_p11 and rel_p10
          rel_psi <- mean(rowSums(as.data.frame(x$z))/x$ncells) #mean of portion of occupied cells across all quadrats
          n_detect_p11 <- 0
          n_detect_p10 <- 0
          for (i in 1:x$nquadrats) {
            for (j in 1:x$ncells) {
              n_detect_p11 <- n_detect_p11 + x$z[i,j]*x$y[i,j]
              n_detect_p10 <- n_detect_p10 + (1-x$z[i,j])*x$y[i,j]
            }
          }
          rel_p11 <- n_detect_p11/(sum(x$z)*x$k) #portion of detections among truly occupied cells (*k because 3 observations of each cell)
          rel_p10 <- n_detect_p10/((x$nquadrats*x$ncells - sum(x$z))*x$k)
          
          inits <- function()list(psi=runif(1,0,1), p11 = runif(1,0,1))
          inits1 <- inits()
          inits2 <- inits()
          inits3 <- inits()
          inits <- list(inits1, inits2, inits3)
          params <- c('psi.mean', 'p11.mean', 'psi','p11', 'sd.lpsi', 'sd.lp11')
          
          adaptSteps = 500 
          burnInSteps = 5000#1000
          nChains = 3 
          numSavedSteps=20000 
          thinSteps=3 
          nIter = 20000#ceiling((numSavedSteps*thinSteps) / nChains) 
          
          #Create, initialize and adapt the model
          jagsModel = jags.model(
            "model_No_FP_het.txt" , data=data , # inits=inits ,
            n.chains=nChains , n.adapt=adaptSteps )
          #check inits: jagsModel$state()
          
          #Burn-in
          cat( "Burning in the MCMC chain...\n" ) 
          update( jagsModel , n.iter=burnInSteps )
          
          cat( "Sampling final MCMC chain...\n" )
          codaSamples = coda.samples( jagsModel ,
                                      variable.names=params, n.iter=nIter,
                                      thin=thinSteps )
          
          
          #Results
          mcmcChain = as.matrix( codaSamples )
          
          estim_psi <- mcmcChain[,'psi.mean']
          estim_p11 <- mcmcChain[,'p11.mean']
          estim_sd_lpsi <- mcmcChain[,'sd.lpsi']
          estim_sd_lp11 <- mcmcChain[,'sd.lp11']
          
          estim_psi_biais <- c(estim_psi_biais, mean(estim_psi) - true_psi) #list with mean estimates of all sites per reps
          mean_estim_psi <- c(mean_estim_psi, mean(estim_psi))
          estim_psi_sd <- c(estim_psi_sd, sd(estim_psi))
          estim_lpsi_sd <- c(estim_lpsi_sd, mean(estim_sd_lpsi))
          estim_p11_biais <- c(estim_p11_biais, mean(estim_p11) - true_p11)
          estim_p11_sd <- c(estim_p11_sd, sd(estim_p11))
          estim_lp11_sd <- c(estim_lp11_sd, mean(estim_sd_lp11))
          
          #CONVERGENCE
          #traplot(nimbleModels, params) #traceplot
          if ({rhat(mcmcChain[,'psi.mean']) <= 1.04} & {rhat(mcmcChain[,'p11.mean']) <= 1.04}) {
            conv <- conv + 1
          }
          
          #95% BAYESIAN CREDIBLE INTERVAL COVER AND WIDTH
          #psi
          if ({true_psi >= ci(estim_psi, method = "HDI")$CI_low} & {true_psi <= ci(estim_psi, method = "HDI")$CI_high}) {
            BCI_cover.psi <- BCI_cover.psi + 1
          }
          BCI_width.psi <- c(BCI_width.psi, ci(estim_psi, method = "HDI")$CI_high - ci(estim_psi, method = "HDI")$CI_low)
          #p11
          if ({true_p11 >= ci(estim_p11, method = "HDI")$CI_low} & {true_p11 <= ci(estim_p11, method = "HDI")$CI_high}) {
            BCI_cover.p11 <- BCI_cover.p11 + 1
          }
          BCI_width.p11 <- c(BCI_width.p11, ci(estim_p11, method = "HDI")$CI_high - ci(estim_p11, method = "HDI")$CI_low)
          
          
          # TIME AND PROGRESS of SIMULATION
          tot.reps <- (length(true_psi.list)*length(true_p11.list))*reps
          avc <- (ct-1)*reps + r
          if( round(avc/tot.reps, 3)-(avc/tot.reps) ==  0 ){ 
            print( paste((avc/reps)*100, "% completed", sep = "") )
          }else{}
          flush.console()
          Sys.sleep(0)
          
        } #end of loop across 100 simulations
        
        ### ESTIMATE PERFORMANCE
        #psi
        biais.psi <- mean(estim_psi_biais)
        error.psi <- mean(estim_psi_sd) #sd of estim_psi
        BCI_cover.psi <- BCI_cover.psi*100/conv #percentage of iterations that have psi in BCI among models that have converged 
        BCI_width.psi <- mean(BCI_width.psi)
        mse.psi <- mean( (estim_psi - true_psi)^2)
        rmse.psi <- sqrt(mse.psi)
        #p11
        biais.p11 <- mean(estim_p11_biais)
        error.p11 <- mean(estim_p11_sd) #sd of estim_p11
        BCI_cover.p11 <- BCI_cover.p11*100/conv #percentage of iterations that have psi in BCI among models that have converged 
        BCI_width.p11 <- mean(BCI_width.p11)
        mse.p11 <- mean( (estim_p11 - true_p11)^2)
        rmse.p11 <- sqrt(mse.p11)
        
        
        TABLE[ct,] <- cbind(true_psi, true_p11, true_p10, rel_psi, rel_p11, rel_p10, conv*100/reps,
                            round(mean(estim_lpsi_sd),3),
                            round(mean(estim_lp11_sd),3),
                            round(biais.psi, 3), round(error.psi, 3), round(BCI_cover.psi, 3), round(BCI_width.psi, 3), round(mse.psi, 3), round(rmse.psi, 3),
                            round(biais.p11, 3), round(error.p11, 3), round(BCI_cover.p11, 3), round(BCI_width.p11, 3), round(mse.p11, 3), round(rmse.p11, 3))
        
        LIST_ESTIM_PSI[ct,] <- c(true_psi, true_p11, true_p10, c(mean_estim_psi))
        
      }
    } 
  }
} #end of loops across all parameters


### GET THE DURATION THAT IT TAKES FOR SIMULATIONS TO COMPLETE
time.end <- Sys.time()
duration <- time.end - time.start 
duration


### SAVE RESULTS
#write.table(TABLE, "simulation_outputs/ResultSimul_het_psi3_p1110.2.csv", row.names = F, sep = ";")
#write.table(LIST_ESTIM_PSI, "simulation_outputs/ResultSimul_het_psi3_p1110_list.2.csv", row.names = F, sep = ";")

### RELOAD TABLE
TABLE <- read_delim("simulation_outputs/ResultSimul_het_psi3_p1110.2.csv", 
                    delim = ";", escape_double = FALSE, trim_ws = TRUE)
LIST_ESTIM_PSI <- read_delim("simulation_outputs/ResultSimul_het_psi3_p1110_list.2.csv", 
                    delim = ";", escape_double = FALSE, trim_ws = TRUE)

### PLOT BIAS

TABLE %>%
  as.data.frame() %>%
  dplyr::select(-c(rel_psi,rel_p11,rel_p10, converged, estim_sd_lpsi, estim_sd_lp11,error.psi, biais.p11:rmse.p11)) %>%
  mutate(biais_psi_prop = round(100*(biais.psi)/true_psi,1)) %>%
  relocate(biais_psi_prop, .after = biais.psi) %>%
  rename('psi' = true_psi,
         'p1|1' = true_p11,
         'p1|0' = true_p10,
         'psi bias' = biais.psi,
         'psi bias (%)' = biais_psi_prop,
         'BCI cover' = BCI_cover.psi,
         'BCI width' = BCI_width.psi,
         'MSE' = mse.psi,
         'RMSE' = rmse.psi) %>%
  kbl(format = 'html') %>%
  kable_classic(full_width = F, html_font = 'LM Roman')


LIST_ESTIM_PSI %>% 
  as.data.frame() %>%
  pivot_longer(cols = !starts_with('true'), names_to = 'run', values_to = 'estim_psi') %>%
  mutate(bias_psi = estim_psi - true_psi) %>%
  mutate(bias_psi_prop = 100*(bias_psi)/true_psi) %>%
  mutate(true_psi = as.factor(true_psi),
         true_p11 = as.factor(true_p11),
         true_p10 = as.factor(true_p10)) %>%
  rename('psi' = true_psi,
         'p1|1' = true_p11,
         'p1|0' = true_p10,
         'psi bias (%)' = bias_psi_prop) %>%
  ggplot(aes(y=`psi bias (%)`, x=`p1|0`))+
  geom_boxplot(notch = TRUE) +
  facet_grid(~psi*`p1|1`, labeller = purrr::partial(label_both, sep = " = ")) +
  geom_hline(yintercept=0, size=0.3) +
  theme_bw()
