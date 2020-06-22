#Author: Agnes Perez
#Date last modification: 22/06/2020
#Longitudinal study with a Frequentist approach


library(lme4)
library(ggplot2)
library(dplyr)
library(lmerTest)
library(sjPlot)

##################################  STEP 0: DATA #################################

# Read the data prepared for the model
adnidata<- read.csv("adnidata.csv", stringsAsFactors = F)

#############################  STEP 1: RANDOM EFFECTS ############################

#RANDOM INTERCEPT MODEL
model1<-lmer(formula = Hippocampus2 ~ cHC*time + sMCI*time + cMCI*time + AD*time + 
               APOE4*time + PTGENDER + AGE + ICV2 + (1 | RID), data = adnidata, REML=FALSE)
summary(model1)

#RANDOM INTERCEPT + SLOPE MODEL
model2<-lmer(formula = Hippocampus2 ~ cHC*time + sMCI*time + cMCI*time + AD*time +  
               APOE4*time + PTGENDER + AGE + ICV2 + (time | RID), data = adnidata, REML = FALSE)
summary(model2)

#ANOVA
anova(model2,model1) 

#############################  STEP 2: THE MODEL ############################

#RANDOM INTERCEPT + SLOPE MODEL
model<-lmer(formula = Hippocampus2 ~ cHC*time + sMCI*time + cMCI*time + AD*time + 
              APOE4*time + PTGENDER + AGE + ICV2 + (time | RID), data = adnidata, REML = FALSE)
summary(model)


#Plot residuals vs fitted
plot(fitted(model), residuals(model), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2, col='red')
lines(smooth.spline(fitted(model), residuals(model)))


###########################  STEP 3: HYPOTHESIS TESTING ###########################

#Hypothesis testing
contest(model, c(0, 0, 0, 0,0,0,0,0,0,0,0,0,0,1,0)) #AD vs sHC
contest(model, c(0, 0, 0,0,0,0,0,0,0,0,0,1,0,0,0)) #sMCI vs sHC
contest(model, c(0, 0, 0,0,0,0,0,0,0,0,0,0,-1,1,0)) #AD vs cMCI
contest(model, c(0, 0, 0,0,0,0,0,0,0,0,0,-1,1,0,0)) #sMCI vs cMCI
contest(model, c(0, 0, 0,0,0,0,0,0,0,0,0,0,1,0,0)) #cMCI vs sHC
contest(model, c(0, 0, 0,0,0,0,0,0,0,0,1,0,0,0,0)) #cHC vs sHC
contest(model, diag(15)[11:14, ]) #sHC vs all

#### Plot CI ####

confint(model) #CI

my_theme <- theme_bw(base_size = 12)  +
  theme(
    plot.title = element_text(hjust = 0.5),# Center title position and size
  )

#plot
plot_model(model, type = "est",color = "bw", title="Confidence Intervals",
           ci_level = 0.5, outer_level = 0.95, dot.size=1.5, axis.lim=c(-2,1)) + my_theme 

#########################  STEP 4: SIMULATION N MCI #########################

#### Simulation ####

N<-c() #N tot
N1<-c() # N sHC
N2<-c() # N cHC
N3<-c() # N sMCI
N4<-c() # N cMCI
N5<-c() # N AD
pvalue<-c()

# 60 iterations
for (j in 0:60){
  a='F'
  cont=0
  adnidataN<-adnidata
  
  #Remove RID until p>0.05  
  while(a!='T'){
    
    RIDout<-sample(adnidataN$RID, 1) #Random RID
    adnidataN<-adnidataN[adnidataN$RID != RIDout, ] #Remove RID selected
    
    #Model
    modelsim<-lmer(formula = Hippocampus2 ~ cHC*time+ sMCI*time + cMCI*time + 
                     AD*time + APOE4*time + PTGENDER + AGE + ICV2 + (time | RID), 
                   data = adnidataN, REML = FALSE)
    
    #Hypothesis testing
    contrast<-contest(modelsim, c(0, 0, 0,0,0,0,0,0,0,0,0,-1,1,0,0))[6] #sMCI vs cMCI
    
    #Stop if we have p>0.05
    if(contrast>0.05){
      pvalueaux<-contrast
      pvalue<-rbind(pvalue,pvalueaux)
      a='T'
    }
    
    cont=cont+1
    
  }
  
  # N values
  
  #N total
  adninor <- adnidataN[!duplicated(adnidataN$RID),]
  Naux<-dim(adninor)[1]
  N<-append(N,Naux)
  
  # N each group
  Naux1<-table(adninor$DX_final)[1]  
  N1<-append(N1,Naux1)
  
  Naux2<-table(adninor$DX_final)[2]  
  N2<-append(N2,Naux2)
  
  Naux3<-table(adninor$DX_final)[3]   
  N3<-append(N3,Naux3)
  
  Naux4<-table(adninor$DX_final)[4]  
  N4<-append(N4,Naux4)
  
  Naux5<-table(adninor$DX_final)[5]   
  N5<-append(N5,Naux5)
  
}

# Dataframe with the results
Ndata<-cbind.data.frame(N,N1)
Ndata<-cbind(Ndata,N2)
Ndata<-cbind(Ndata,N3)
Ndata<-cbind(Ndata,N4)
Ndata<-cbind(Ndata,N5)
Ndata<-cbind(Ndata,pvalue)

#save results
write.csv(Ndata, file="simdataMCI.csv")

#### Plots ####

#read simulation results data
simdataMCI<- read.csv("simdataMCI.csv", stringsAsFactors = F, sep=';')

simdataMCI$Nnew<-simdataMCI$N+1 #minimum N

# Histogram
hist(simdataMCI$Nnew, main='Histogram of the simulation sMCI vs cMCI', 
     xlab='N', col='skyblue1')

# Outlier
IQR<-summary(simdataMCI$N)[5]-summary(simdataMCI$N)[2]
Q3<-summary(simdataMCI$N)[5]
1.5*IQR+Q3

#Remove Outlier
simdataMCI<-subset(simdataMCI,simdataMCI$N<500)
simdataMCI$Nnew<-simdataMCI$N+1 #minimum N

# Histogram
hist(simdataMCI$Nnew, main='Histogram of the simulation sMCI vs cMCI', 
     xlab='N', col='skyblue1')

#Results
round(median(simdataMCI$Nnew),0)
round(summary(simdataMCI$Nnew)[2],0) #Q1
round(summary(simdataMCI$Nnew)[5],0) #Q3


#########################  STEP 5: SIMULATION N ALL #########################

#### Simulation ####

N<-c() #N tot
N1<-c() # N sHC
N2<-c() # N cHC
N3<-c() # N sMCI
N4<-c() # N cMCI
N5<-c() # N AD
pvalue<-c()

# 60 iterations
for (j in 0:60){
  a='F'
  cont=0
  adnidataN<-adnidata
  
  #Remove RID until p>0.05  
  while(a!='T'){
    
    RIDout<-sample(adnidataN$RID, 1) #Random RID
    adnidataN<-adnidataN[adnidataN$RID != RIDout, ] #Remove RID selected
    
    #Model
    modelsim<-lmer(formula = Hippocampus2 ~ cHC*time+ sMCI*time + cMCI*time + 
                     AD*time + APOE4*time + PTGENDER + AGE + ICV2 + (time | RID), 
                   data = adnidataN, REML = FALSE)
    
    #Hypothesis testing
    contrast<-contest(modelsim, diag(15)[11:14, ])[6] #sHC vs all
    
    #Stop if we have p>0.05
    if(contrast>0.05){
      pvalueaux<-contrast
      pvalue<-rbind(pvalue,pvalueaux)
      a='T'
    }
    
    cont=cont+1
    
  }
  
  # N values
  
  #N total
  adninor <- adnidataN[!duplicated(adnidataN$RID),]
  Naux<-dim(adninor)[1]
  N<-append(N,Naux)
  
  # N each group
  Naux1<-table(adninor$DX_final)[1]  
  N1<-append(N1,Naux1)
  
  Naux2<-table(adninor$DX_final)[2]  
  N2<-append(N2,Naux2)
  
  Naux3<-table(adninor$DX_final)[3]   
  N3<-append(N3,Naux3)
  
  Naux4<-table(adninor$DX_final)[4]  
  N4<-append(N4,Naux4)
  
  Naux5<-table(adninor$DX_final)[5]   
  N5<-append(N5,Naux5)
  
}

# Dataframe with the results
Ndata<-cbind.data.frame(N,N1)
Ndata<-cbind(Ndata,N2)
Ndata<-cbind(Ndata,N3)
Ndata<-cbind(Ndata,N4)
Ndata<-cbind(Ndata,N5)
Ndata<-cbind(Ndata,pvalue)

#save results
write.csv(Ndata, file="simdataall.csv")

#### Plots ####

#read simulation results data
simdataALL<- read.csv("simreddataall.csv", stringsAsFactors = F, sep=';')

simdataALL$Nnew<-simdataALL$N+1 #minimum N

#Histogram
hist(simdataALL$Nnew, main='Histogram of the simulation of all groups', 
     xlab='N', col='salmon')

#Results
round(median(simdataALL$Nnew),0)
round(summary(simdataALL$Nnew)[2],0) #Q1
round(summary(simdataALL$Nnew)[5],0) #Q3


###########################  STEP 6: REDUCED DATA ###########################

#We only keep the subjects that have data for the 4 time points

adnidata1<-subset(adnidata,adnidata$time==0.0)
adnidata2<-subset(adnidata,adnidata$time==0.5)
adnidata3<-subset(adnidata,adnidata$time==1.0)
adnidata4<-subset(adnidata,adnidata$time==2.0)


#Time point 1 + timepoint 2
data1<-adnidata1
data2<-adnidata2

data1$comp <- ifelse(adnidata1$RID %in% adnidata2$RID, 1, 0)
data2$comp <- ifelse(adnidata2$RID %in% adnidata1$RID, 1, 0)

dataaux<-rbind(data1,data2)

dataaux$comp<-as.factor(dataaux$comp)
data<-subset(dataaux,dataaux$comp==1) #T1 + T2
data<-select(data,-comp)

#Time point 1 + timepoint 2 +  timepoint 3
data3<-adnidata3
dataaux<-data

dataaux$comp <- ifelse(dataaux$RID %in% data3$RID, 1, 0)
data3$comp <- ifelse(data3$RID %in% dataaux$RID, 1, 0)

dataaux<-rbind(dataaux,data3)

dataaux$comp<-as.factor(dataaux$comp)
data<-subset(dataaux,dataaux$comp==1) #T1 + T2 +T3
data<-select(data,-comp)


#Time point 1 + timepoint 2 +  timepoint 3 + timepoint 4
data4<-adnidata4
dataaux<-data

dataaux$comp <- ifelse(dataaux$RID %in% data4$RID, 1, 0)
data4$comp <- ifelse(data4$RID %in% dataaux$RID, 1, 0)

dataaux<-rbind(dataaux,data4)

dataaux$comp<-as.factor(dataaux$comp)
data<-subset(dataaux,dataaux$comp==1) #T1 + T2 + T3 + T4
data<-select(data,-comp)

#save data
write.csv(data, file="adinidatared.csv")

###########################  STEP 7: SIMULATION NA MCI ###########################

#### Simulation ####

N<-c() #N tot
N1<-c() # N sHC
N2<-c() # N cHC
N3<-c() # N sMCI
N4<-c() # N cMCI
N5<-c() # N AD
pvalue<-c()

# 60 iterations
for (j in 0:60){
  print(j)
  a='F'
  cont=0
  adnidataN<-data
  
  #Remove Time point until p>0.05  
  while(a!='T'){
    
    RIDout<-sample(adnidataN$X, 1) #Random time point
    adnidataN<-adnidataN[adnidataN$X != RIDout, ] #Remove time point
    
    #Break due to limitation
    if(dim(adnidataN)[1]==length(unique(adnidataN$RID))*2 | 
       (dim(adnidataN)[1]<length(unique(adnidataN$RID))*2)){a<-'T'}
    
    if(dim(adnidataN)[1]>length(unique(adnidataN$RID))*2){
      #Model
      modelsim<-lmer(formula = Hippocampus2 ~ cHC*time+ sMCI*time + cMCI*time + 
                       AD*time + APOE4*time + PTGENDER + AGE + ICV2 + (time | RID), 
                     data = adnidataN, REML = FALSE)
      
      #Hypothesis testing
      contrast<-contest(modelsim, c(0, 0, 0,0,0,0,0,0,0,0,0,-1,1,0,0))[6] #sMCI vs cMCI
      
      #Stop if we have p>0.05
      if(contrast>0.05){
        print('break')
        break
      }
      
    }
    
    cont=cont+1
  }
  
  
  #N total
  adninor <- adnidataN[!duplicated(adnidataN$RID),]
  Naux<-dim(adninor)[1]
  N<-append(N,Naux)
  
  # N each group
  Naux1<-table(adninor$DX_final)[1] 
  N1<-append(N1,Naux1)
  
  Naux2<-table(adninor$DX_final)[2] 
  N2<-append(N2,Naux2)
  
  Naux3<-table(adninor$DX_final)[3]  
  N3<-append(N3,Naux3)
  
  Naux4<-table(adninor$DX_final)[4]   
  N4<-append(N4,Naux4)
  
  Naux5<-table(adninor$DX_final)[5]  
  N5<-append(N5,Naux5)
  
  pvalue<-contrast
  
}


# Dataframe with the results
Ndata<-cbind.data.frame(N,N1)
Ndata<-cbind(Ndata,N2)
Ndata<-cbind(Ndata,N3)
Ndata<-cbind(Ndata,N4)
Ndata<-cbind(Ndata,N5)
Ndata<-cbind(Ndata,pvalue)

#save results
write.csv(Ndata, file="simreddataMCI.csv")

#### Plots ####

#read simulation results data
simdataMCI<- read.csv("simreddataMCI.csv", stringsAsFactors = F, sep=';')

#Histogram
hist(simdataMCI$N, main='Histogram of the simulation sMCI vs cMCI', 
     xlab='N', col='skyblue1')


###########################  STEP 8: SIMULATION NA ALL ###########################

N<-c() #N tot
N1<-c() # N sHC
N2<-c() # N cHC
N3<-c() # N sMCI
N4<-c() # N cMCI
N5<-c() # N AD
pvalue<-c()

# 60 iterations
for (j in 0:60){
  print(j)
  a='F'
  cont=0
  adnidataN<-data
  
  #Remove Time point until p>0.05  
  while(a!='T'){
    
    RIDout<-sample(adnidataN$X, 1) #Random time point
    adnidataN<-adnidataN[adnidataN$X != RIDout, ] #Remove time point
    
    #Break due to limitation
    if(dim(adnidataN)[1]==length(unique(adnidataN$RID))*2 | 
       (dim(adnidataN)[1]<length(unique(adnidataN$RID))*2)){a<-'T'}
    
    if(dim(adnidataN)[1]>length(unique(adnidataN$RID))*2){
      #Model
      modelsim<-lmer(formula = Hippocampus2 ~ cHC*time+ sMCI*time + cMCI*time + 
                       AD*time + APOE4*time + PTGENDER + AGE + ICV2 + (time | RID), 
                     data = adnidataN, REML = FALSE)
      
      #Hypothesis testing
      contrast<-contest(modelsim, diag(15)[11:14, ])[6] #sHC vs all
      
      #Stop if we have p>0.05
      if(contrast>0.05){
        print('break')
        break
      }
      
    }
    
    cont=cont+1
  }
  
  
  #N total
  adninor <- adnidataN[!duplicated(adnidataN$RID),]
  Naux<-dim(adninor)[1]
  N<-append(N,Naux)
  
  # N each group
  Naux1<-table(adninor$DX_final)[1] 
  N1<-append(N1,Naux1)
  
  Naux2<-table(adninor$DX_final)[2] 
  N2<-append(N2,Naux2)
  
  Naux3<-table(adninor$DX_final)[3]  
  N3<-append(N3,Naux3)
  
  Naux4<-table(adninor$DX_final)[4]   
  N4<-append(N4,Naux4)
  
  Naux5<-table(adninor$DX_final)[5]  
  N5<-append(N5,Naux5)
  
  pvalue<-contrast
  
}


# Dataframe with the results
Ndata<-cbind.data.frame(N,N1)
Ndata<-cbind(Ndata,N2)
Ndata<-cbind(Ndata,N3)
Ndata<-cbind(Ndata,N4)
Ndata<-cbind(Ndata,N5)
Ndata<-cbind(Ndata,pvalue)

#save results
write.csv(Ndata, file="simreddataall.csv")

#### Plots ####

#read simulation results data
simdataALL<- read.csv("simreddataall.csv", stringsAsFactors = F, sep=';')

#Histogram
hist(simdataALL$N, main='Histogram of the simulation of all groups', 
     xlab='N', col='salmon')

###########################  STEP 9: NA BAYESIAN ###########################

#Creation of a dataset where a frequentist approach is not possible

a='F'
adnidataN<-data

while(a!='T'){
  
  RIDout<-sample(adnidataN$X, 1) #Random time point
  adnidataN<-adnidataN[adnidataN$X != RIDout, ] #Remove time point
  
  #stop when we have a dataset of the limitation
  if(dim(adnidataN)[1]==length(unique(adnidataN$RID))*2 |
     (dim(adnidataN)[1]<length(unique(adnidataN$RID))*2)){a<-'T'}
}

#To check that the frequentist approximation does not work
modelsim<-lmer(formula = Hippocampus2 ~ cHC*time+ sMCI*time + cMCI*time + 
                 AD*time + APOE4*time + PTGENDER + AGE + ICV2 + (time | RID), 
               data = adnidataN, REML = FALSE)

#save dataset of the limitation
write.csv(adnidataN, file="adnidataNAextrem.csv")

######################  STEP 10: N MCI FREQ vs BAYESIAN ######################

# We want to compare the same dataset for the two approximations
# We change the iterations to obtain different datasets
# These lines is for obtain the dataset of 301 subjects
# For example to obtain 305 it would be: for(i in 0:944)

adnidataN<-adnidata

for (i in 0:948){
  RIDout<-sample(adnidataN$RID, 1) #Random RID
  adnidataN<-adnidataN[adnidataN$RID != RIDout, ] #Remove RID
}

#Model for the frequentist results
modelsim<-lmer(formula = Hippocampus2 ~ cHC*time+ sMCI*time + cMCI*time + AD*time 
               + APOE4*time + PTGENDER + AGE + ICV2 + (time | RID), 
               data = adnidataN, REML = FALSE)

#Hypothesis testing
contest(modelsim, c(0, 0, 0,0,0,0,0,0,0,0,0,-1,1,0,0)) #sMCI vs cMCI
#contest(modelsim, diag(15)[11:14, ])[6] #sHC vs all

#N
adninor <- adnidataN[!duplicated(adnidataN$RID),]#Baseline RID
table(adninor$DX_final) # N for each group

#save dataset
write.csv(adnidataN, file="adnidata301.csv")

######################  STEP 11: N ALL FREQ vs BAYESIAN ######################

# We want to compare the same dataset for the two approximations
# These lines is for obtain the dataset of 199 subjects

adnidataN<-adnidata

for (i in 0:1050){
  RIDout<-sample(adnidataN$RID, 1) #Random RID
  adnidataN<-adnidataN[adnidataN$RID != RIDout, ] #Remove RID
}

#Model for the frequentist results
modelsim<-lmer(formula = Hippocampus2 ~ cHC*time+ sMCI*time + cMCI*time + AD*time 
               + APOE4*time + PTGENDER + AGE + ICV2 + (time | RID), 
               data = adnidataN, REML = FALSE)

#Hypothesis testing
contest(modelsim, diag(15)[11:14, ])[6] #sHC vs all

#N
adninor <- adnidataN[!duplicated(adnidataN$RID),]#Baseline RID
table(adninor$DX_final) # N for each group

#save dataset
write.csv(adnidataN, file="adnidata199.csv")
