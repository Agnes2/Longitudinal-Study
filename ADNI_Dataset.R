#Author: Agnes Perez
#Date last modification: 22/06/2020
#ADNI Dataset analysis


library(ggplot2)
library(dplyr)

#Function to delete NA
delete.na <- function(df, n=0) {
  df[rowSums(is.na(df)) <= n,]
}

##################################  STEP 0: DATA #################################

#Read the data download from ADNI
adnioriginal<- read.csv("dadeslogvariablesRM.csv", stringsAsFactors = F)

#Ral AGE subjects
RealAge<-adnioriginal$AGE+adnioriginal$Years_bl
adnioriginal<-cbind(adnioriginal,RealAge)

#Select the most important variables for this study
adni <- select(adnioriginal, RID, DX_bl, DX, RealAge,Years_bl, AGE, ICV,APOE4, 
               PTGENDER,sMCI, cMCI, sHC, cHC, AD,Ventricles,Hippocampus,
               WholeBrain,Entorhinal,Fusiform,MidTemp,Ventricles_bl,
               Hippocampus_bl,WholeBrain_bl,Entorhinal_bl,Fusiform_bl,MidTemp_bl,
               ICV_bl,Month_bl, Month, M) 
head(adni)


################################  STEP 1: VARIABLES ###############################

#Diagnostic labels as numeric variables
#Legend 0:CN, 1:AD i 2:MCI
adni$DX[adni$DX == "CN"] <- 0
adni$DX[adni$DX == "AD"] <- 1
adni$DX[adni$DX == "MCI"] <- 2

adni$DX_bl[adni$DX_bl == "CN"] <- 0
adni$DX_bl[adni$DX_bl == "AD"] <- 1
adni$DX_bl[adni$DX_bl == "MCI"] <- 2

adni$DX<-as.numeric(adni$DX)
adni$DX_bl<-as.numeric(adni$DX_bl)


#Sex as numeric variable
#Legend 0:M, 1:F
adni$PTGENDER[adni$PTGENDER == "Male"] <- 0
adni$PTGENDER[adni$PTGENDER == "Female"] <- 1

adni$PTGENDER<-as.numeric(adni$PTGENDER)

#APOE-E4 variable
#Legend carriers: 1+2 no carriers:0 --> carriers: 1 no carriers:0
adni$APOE4[adni$APOE4 == 2] <- 1


#New diagnostic labels with the five clinical groups
#Legend sHC: 1, cHC: 2, sMCI: 3, cMCI: 4, AD: 5
DX_final<-(adni$sHC)+2*(adni$cHC)+3*(adni$sMCI)+4*(adni$cMCI)+5*(adni$AD)
adni<-cbind(DX_final,adni)


##############################  STEP 2: HV ##############################

#Select variables needed
adnidata <- select(adni, RID, DX_final, DX_bl, DX, RealAge,Years_bl, AGE, ICV,APOE4, 
                   PTGENDER,sMCI, cMCI, sHC, cHC, AD,Hippocampus,Hippocampus_bl,
                   ICV_bl,Month_bl, Month, M) 
head(adnidata)

#To temove subjects that present NA
adnidata<-delete.na(adnidata) 

#To remove subjects that are no longer longitudinal due to the elimination of NA
x <- adnidata$RID[duplicated(adnidata$RID)]
unique(x)

long <- adnidata$RID %in% c(x) #RID TRUE-->duplicated, longitudinals

adnidata<-cbind(long,adnidata)
adnidata$long<-as.factor(adnidata$long)
adnidata<-subset(adnidata,adnidata$long==TRUE) #keep TRUE


############################  STEP 3: ALL DEMOGRAPHICS ############################

adninor <- adnidata[!duplicated(adnidata$RID),]#Baseline RID

adninor1<-subset(adninor,adninor$DX_final==1) #Baseline sHC
adninor2<-subset(adninor,adninor$DX_final==2) #Baseline cHC
adninor3<-subset(adninor,adninor$DX_final==3) #Baseline sMC
adninor4<-subset(adninor,adninor$DX_final==4) #Baseline cMCI
adninor5<-subset(adninor,adninor$DX_final==5) #Baseline AD

#N each group:
table(adninor$DX_final) 

estudi<-adninor4$AGE #change (choose study)

#Normality study
qqnorm(estudi, main = "HC")
qqline(estudi)
shapiro.test(estudi) #shapiro test


#### p-value ####

#APOE4: fisher test
apoe4table<-t(matrix(c(table(adninor1$APOE4),table(adninor2$APOE4),
                       table(adninor3$APOE4),table(adninor4$APOE4),
                       table(adninor5$APOE4)),2,5))

fisher.test(apoe4table,simulate.p.value=TRUE)

#SEX: fisher test
sextable<-t(matrix(c(table(adninor1$PTGENDER),table(adninor2$PTGENDER),
                     table(adninor3$PTGENDER),table(adninor4$PTGENDER),
                     table(adninor5$PTGENDER)),2,5))

fisher.test(sextable,simulate.p.value=TRUE)

#AGE:
anova<-aov(AGE ~ DX_final, data=adninor)
summary(anova)


#### plots ####

#Plot points REAL AGE: 3 grups
ggplot(data  = adnidata,
       aes(x = RealAge,
           y = Hippocampus,
           col = DX))+
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter")+
  theme_minimal()+
  scale_color_gradientn(colours = c('red','blue','green'))+
  labs(title = "Hippocampus for each group", x= 'Age (years)', 
       y= 'Hippocampus (mm3)')+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  theme(legend.position='none')


#Plot points REAL AGE: 5 grups
ggplot(data  = adnidata,
       aes(x = RealAge,
           y = Hippocampus,
           col = DX_final))+
  geom_point(size = 1.2,
             alpha = .8,
             position = "jitter")+
  theme_minimal()+
  scale_color_gradientn(colours = c('red','turquoise1','blue','yellow','green'))+
  labs(title = "Hippocampus for each group", x= 'Age (years)', 
       y= 'Hippocampus (mm3)')+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))+
  theme(legend.position='none')


#### all time points ####

table(adnidata$M) #Time point N

#time points
adnit0<-subset(adnidata,adnidata$M==0) #time point 0 (t0)
adnit1<-subset(adnidata,adnidata$M==6) #t1
adnit2<-subset(adnidata,adnidata$M==12) #t2
adnit3<-subset(adnidata,adnidata$M==18) #t3
adnit4<-subset(adnidata,adnidata$M==24) #t4
adnit5<-subset(adnidata,adnidata$M==36) #t5
adnit6<-subset(adnidata,adnidata$M==48) #t6
adnit7<-subset(adnidata,adnidata$M==60) #t7
adnit8<-subset(adnidata,adnidata$M==72) #t8
adnit9<-subset(adnidata,adnidata$M==84) #t9
adnit10<-subset(adnidata,adnidata$M==96) #t10
adnit11<-subset(adnidata,adnidata$M==108) #t11
adnit12<-subset(adnidata,adnidata$M==120) #t12

estudi<-adnit0$DX_final #change (choose study)
table(estudi) #N time point according to group


estudi<-adnit12$Years_bl #change (choose study)
mean(estudi) # mean time point
sd(estudi) # std time point


############################  STEP 4: 4 TIME POINTS ############################

#Select the 4 time points
adnidata<-subset(adnidata,adnidata$M==0 | adnidata$M==6 | 
                   adnidata$M==12 | adnidata$M==24) 

#To remove subjects that are no longer longitudinal 
#due to the elimination of time points
x <- adnidata$RID[duplicated(adnidata$RID)]
unique(x)

long1 <- adnidata$RID %in% c(x) #RID TRUE-->duplicated, longitudinals

adnidata<-cbind(long1,adnidata)
adnidata$long1<-as.factor(adnidata$long1)
adnidata<-subset(adnidata,adnidata$long1==TRUE) #keep TRUE


##########################  STEP 5: STUDY DEMOGRAPHICS ##########################

adninor <- adnidata[!duplicated(adnidata$RID),]#Baseline RID

adninor1<-subset(adninor,adninor$DX_final==1) #Baseline sHC
adninor2<-subset(adninor,adninor$DX_final==2) #Baseline cHC
adninor3<-subset(adninor,adninor$DX_final==3) #Baseline sMC
adninor4<-subset(adninor,adninor$DX_final==4) #Baseline cMCI
adninor5<-subset(adninor,adninor$DX_final==5) #Baseline AD

#N each group:
table(adninor$DX_final) 

estudi<-adninor4$AGE #change (choose study)

#Normality study
qqnorm(estudi, main = "HC")
qqline(estudi)
shapiro.test(estudi) #shapiro test

mean(estudi) # mean AGE
sd(estudi) # std AGE

#### p-value ####

#APOE4: fisher test
apoe4table<-t(matrix(c(table(adninor1$APOE4),table(adninor2$APOE4),
                       table(adninor3$APOE4),table(adninor4$APOE4),
                       table(adninor5$APOE4)),2,5))

fisher.test(apoe4table,simulate.p.value=TRUE)

#SEX: fisher test
sextable<-t(matrix(c(table(adninor1$PTGENDER),table(adninor2$PTGENDER),
                     table(adninor3$PTGENDER),table(adninor4$PTGENDER),
                     table(adninor5$PTGENDER)),2,5))

fisher.test(sextable,simulate.p.value=TRUE)

#AGE:
anova<-aov(AGE ~ DX_final, data=adninor)
summary(anova)


#### N Time points according to clinical group ####

table(adnidata$M) #N de cada un dels timepoints

adnit0<-subset(adnidata,adnidata$M==0) #time point 0 (t0)
adnit1<-subset(adnidata,adnidata$M==6) #t1
adnit2<-subset(adnidata,adnidata$M==12) #t2
adnit3<-subset(adnidata,adnidata$M==24) #t3

estudi<-adnit0$DX_final #change (choose study)
table(estudi) #N time point according to DX


#### Time points with Years to Baseline ####

#New variable time: Years to Baseline for each time point
adnidata<-cbind(adnidata,adnidata$M) 
noms=names(adnidata)
noms[24]='time'
names(adnidata) = noms

adnidata$time[adnidata$time == 0] <- 0.0
adnidata$time[adnidata$time == 6] <- 0.5
adnidata$time[adnidata$time == 12] <- 1.0
adnidata$time[adnidata$time == 24] <- 2.0



#### Boxplot + post-hoc ####

HVt1<-subset(adnidata,adnidata$time==0.0)
HVt2<-subset(adnidata,adnidata$time==0.5)
HVt3<-subset(adnidata,adnidata$time==1.0)
HVt4<-subset(adnidata,adnidata$time==2.0)

#Boxplot
boxplot(Hippocampus ~ DX_final, data = HVt1, xlab='Group', 
        ylab='Hippocampus (mm3)', main='Year 0',
        col=(c('red','turquoise1','blue','yellow','green')))

#post-hoc
pairwise.t.test(HVt3$Hippocampus, HVt3$DX_final, p.adj = "BH")


##########################  STEP 6: DATA FOR THE MODEL ##########################

#Scale ICV and HV
adnidata$ICV2<-scale(adnidata$ICV) 
adnidata$Hippocampus2<-scale(adnidata$Hippocampus) 

#save data
write.csv(adnidata, file="adnidata.csv")