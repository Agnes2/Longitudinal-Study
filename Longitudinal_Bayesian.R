#Author: Agnes Perez
#Date last modification: 22/06/2020
#Longitudinal study with a Bayesian approach

library(rstan)
library(ggplot2)
library(dplyr)


##################################  STEP 0: DATA #################################
# NOT RUN ALL THE OPTIONS, SELECT ONLY ONE

# ORIGINAL MODEL
adnidata<- read.csv("adnidata.csv", stringsAsFactors = F)

# SIM N MCI:
adnidata<- read.csv("adnidata301.csv", stringsAsFactors = F) 

# SIM N ALL:
adnidata<- read.csv("adnidata199.csv", stringsAsFactors = F) 

# SIM NA:
adnidata<- read.csv("adnidataNAextrem.csv", stringsAsFactors = F)
#The subjects must be sorted
adnidata <- adnidata[with(adnidata, order(adnidata$RID)), ] 


# WITH THE DATA SELECTED WE CAN CARRY ON WITH THE OTHERS STEPS, 
# WHICH ARE EQUAL FOR ALL THE DOCUMENTS


#Rename subjects IDs to have 1,2,3, ...
RIDold<-adnidata$RID
adnidata<-cbind(adnidata,RIDold)
A<-unique(adnidata$RID)
for (i in 1:length(A)){
  adnidata$RID[adnidata$RID == A[i]] <- i
}

#Select just the variables we will use
adnidatab <- select(adnidata, RID, AGE, ICV2,APOE4, PTGENDER, sMCI, cMCI, sHC, cHC, AD,
                    Hippocampus2,time) 

#################### STEP 1: OBTAIN THE PARAMETERS NAMES ############################

#Create a linear model to obtain the parameters names
lmod <- lm(Hippocampus2 ~ cHC*time + sMCI*time + cMCI*time + AD*time + APOE4*time 
           + PTGENDER + AGE + ICV2, adnidatab)
#Parameters names
x <- model.matrix(lmod)

################# STEP 2: CREATE THE DATA FOR THE BAYESIAN LME MODEL #################

adniddat <- list(Nobs = nrow(adnidatab),
                 Npreds = ncol(x),
                 Ngroups = length(unique(adnidatab$RID)),
                 y = adnidatab$Hippocampus2,
                 x = x,
                 timevar = adnidatab$time,
                 group = adnidatab$RID)

############################# STEP 3: STAN FILE ######################################

#Read the stan document, previously created
stanc("longitudinal.stan")
stan_model1 <- "longitudinal.stan"

################################ STEP 4: FIT THE MODEL ##############################

fit=stan(file = stan_model1, data = adniddat, warmup = 4250, iter = 8500, chains = 4, 
         cores = 2, thin = 1)
saveRDS(fit, file = "fit.rds") #save the model

################################ STEP 5: CONVERGENCE ################################

#Summary of the model 
print(fit,par=c("beta","sigmaint","sigmaslope","sigmaeps"))

#Traceplot
traceplot(fit, par=c("beta"))

########################### STEP 6: POSTERIOR DISTRIBUTION ###########################

#CrI plot
plot(fit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, 
     fill_color = "black",par=c("beta[2]","beta[3]","beta[4]","beta[5]",
                                "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]","beta[11]",
                                "beta[12]","beta[13]","beta[14]","beta[15]")) 
+ ggtitle('Credible Intervals') + xlim(-2,1)
