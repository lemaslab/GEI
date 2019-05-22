# This code is intended to analyze the CANHR data with respect to the 
# PUFA SNPs analyses

# Plot Interactions
# rs11XX80->chol p=0.006
# rs59XX61-->WHR p=0.006
# rs7849-->WHR p=0.007
# rs22XX80-->HOMAIR p=0.002

# **************************************************************************** #
# ***************                Directory Variables           *************** #
# **************************************************************************** #

# Set working directory
# Work Computer
setwd("C:/Users/lemasd/Dropbox/Analysis/CANHR/R/PUFA_SNPs/Data")
# Laptop
#setwd("C:/Users/Dominick/Dropbox/Analysis/CANHR/R/PUFA_SNPs/Data")

# Clear the slate #
rm(list=ls())

# **************************************************************************** #
# ***************                Import Data                   *************** #
# **************************************************************************** #

data= read.table("AssocPed_Full3_10Oct13.dat", header=T,na="",sep="\t",as.is=TRUE)


# **************************************************************************** #
# ***************                  Library                     *************** #
# **************************************************************************** #

library("coxme")
library("plyr")

# **************************************************************************** #
# ***************     Subset and recode data.frames            *************** #
# **************************************************************************** #

# Subset the canhr data.frame
# List variables in canhr
names(data)

# Create canhr1 and include only pedigree and phenotypes variables
# Pedigree---> 1=Fam, 2=Ind, 3=FA, 4=MO, 5=Sex, 464=age, 467:469=A1-A2
# n-3 PUFA---> 465=d15n, 470:472=DN1-DN2
# Phenotype--> 376:463
# VillageGroup > 476
# canhr-id > 655
# Genotype--> 132:133= rs599961 alleles & 142:143 rs7849 allels
# Genotype--> 120:121 rs1502593 alleles & 126:127 rs3071 alleles

canhr=data
names(canhr)

# TG_HDL Ratio
tg_hdl=canhr$lab_triglyceride/canhr$lab_dhdl
hist(tg_hdl)
mean(tg_hdl, na.rm=T)
sd(tg_hdl, na.rm=T)
which(tg_hdl>10)
tg_hdl[370]=NA
tg_hdl[1978]=NA
plot(canhr$hba1c,tg_hdl)


# Leptin_Adiponectin Ratio
hist(canhr$lep_adipo)
hist(canhr$lep_adipo_log)
hist(canhr$lep_adipo_sqrt)


#### Check WHR
# Waist to Hip Ratio untransformed
canhr$WHR.c=canhr$waistcirc_avg/canhr$thighcirc_avg
#plot(canhr$WHR,canhr$WHR.c) # General Agreement with outlier: c-1182
# Waist to Hip Ratio transformed
#plot(canhr$WHR,canhr$WHRp) # Spot on.

# Concatenate the SNPs into a single variable
rs599961=paste(canhr$rs599961_A1, canhr$rs599961_A2, sep="");table(rs599961) #T>G
# Condense the variable down into a single catagorical variable
canhr$rs599961[rs599961==00]=NA;canhr$rs599961[rs599961==44]=1
canhr$rs599961[rs599961==43]=2;canhr$rs599961[rs599961==33]=3;table(canhr$rs599961)

# Concatenate the SNPs into a single variable
rs1502593=paste(canhr$rs1502593_A1,canhr$rs1502593_A2,sep="");table(rs1502593) # G>A
# Condense the variable down into a single catagorical variable
canhr$rs1502593[rs1502593==00]=NA; canhr$rs1502593[rs1502593==33]=1
canhr$rs1502593[rs1502593==13]=2;canhr$rs1502593[rs1502593==11]=3;table(canhr$rs1502593)

# Concatenate the SNPs into a single variable
rs3071=paste(canhr$rs3071_A1,canhr$rs3071_A2,sep="");table(rs3071) #A>C
# Condense the variable down into a single catagorical variable
canhr$rs3071[rs3071==00]=NA;canhr$rs3071[rs3071==11]=1
canhr$rs3071[rs3071==21]=2;canhr$rs3071[rs3071==22]=3;table(canhr$rs3071)

# Concatenate the SNPs into a single variable
rs7849=paste(canhr$rs7849_A1, canhr$rs7849_A2, sep="");table(rs7849) #T>C
# Condense the variable down into a single catagorical variable
canhr$rs7849[rs7849==00]=NA;canhr$rs7849[rs7849==44]=1
canhr$rs7849[rs7849==42]=2;canhr$rs7849[rs7849==22]=3;table(canhr$rs7849)

# Concatenate the SNPs into a single variable
rs11190480=paste(canhr$rs11190480_A1, canhr$rs11190480_A2, sep="");table(rs11190480) #A>G
# Condense the variable down into a single catagorical variable
canhr$rs11190480[rs11190480==00]=NA;canhr$rs11190480[rs11190480==11]=1
canhr$rs11190480[rs11190480==31]=2;canhr$rs11190480[rs11190480==33]=3;table(canhr$rs11190480)

# Concatenate the SNPs into a single variable
rs2167444=paste(canhr$rs2167444_A1, canhr$rs2167444_A2, sep="");table(rs2167444) #T>A
# Condense the variable down into a single catagorical variable
canhr$rs2167444[rs2167444==00]=NA;canhr$rs2167444[rs2167444==44]=1
canhr$rs2167444[rs2167444==14]=2;canhr$rs2167444[rs2167444==11]=3;table(canhr$rs2167444)

# Concatenate the SNPs into a single variable
rs3829160=paste(canhr$rs3829160_A1, canhr$rs3829160_A2, sep="");table(rs3829160) #G>A
# Condense the variable down into a single catagorical variable
canhr$rs3829160[rs3829160==00]=NA;canhr$rs3829160[rs3829160==33]=1
canhr$rs3829160[rs3829160==13]=2;canhr$rs3829160[rs3829160==11]=3;table(canhr$rs3829160)

# Concatenate the SNPs into a single variable
rs2234970=paste(canhr$rs2234970_A1, canhr$rs2234970_A2, sep="");table(rs2234970) #A>C
# Condense the variable down into a single catagorical variable
canhr$rs2234970[rs2234970==00]=NA;canhr$rs2234970[rs2234970==11]=1
canhr$rs2234970[rs2234970==21]=2;canhr$rs2234970[rs2234970==22]=3;table(canhr$rs2234970)

# Concatenate the SNPs into a single variable
rs11557927=paste(canhr$rs11557927_A1, canhr$rs11557927_A2, sep="");table(rs11557927) #T>G
# Condense the variable down into a single catagorical variable
canhr$rs11557927[rs11557927==00]=NA;canhr$rs11557927[rs11557927==44]=1
canhr$rs11557927[rs11557927==34]=2;canhr$rs11557927[rs11557927==33]=3;table(canhr$rs11557927)

# Concatenate the SNPs into a single variable
rs3978768=paste(canhr$rs3978768_A1, canhr$rs3978768_A2, sep="");table(rs3978768) #A>G
# Condense the variable down into a single catagorical variable
canhr$rs3978768[rs3978768==00]=NA;canhr$rs3978768[rs3978768==11]=1
canhr$rs3978768[rs3978768==31]=2;canhr$rs3978768[rs3978768==33]=3;table(canhr$rs3978768)

# Concatenate the SNPs into a single variable
rs5435=paste(canhr$rs5435_A1, canhr$rs5435_A2, sep="");table(rs5435) #T>C
# Condense the variable down into a single catagorical variable
canhr$rs5435[rs5435==00]=NA;canhr$rs5435[rs5435==44]=1
canhr$rs5435[rs5435==24]=2;canhr$rs5435[rs5435==22]=3;table(canhr$rs5435)

# Concatenate the SNPs into a single variable
rs5415=paste(canhr$rs5415_A1, canhr$rs5415_A2, sep="");table(rs5415) #T>C
# Condense the variable down into a single catagorical variable
canhr$rs5415[rs5415==00]=NA;canhr$rs5415[rs5415==44]=1
canhr$rs5415[rs5415==24]=2;canhr$rs5415[rs5415==22]=3;table(canhr$rs5415)

# Concatenate the SNPs into a single variable
rs4925114=paste(canhr$rs4925114_A1, canhr$rs4925114_A2, sep="");table(rs4925114) #G>A
# Condense the variable down into a single catagorical variable
canhr$rs4925114[rs4925114==00]=NA;canhr$rs4925114[rs4925114==33]=1
canhr$rs4925114[rs4925114==13]=2;canhr$rs4925114[rs4925114==11]=3;table(canhr$rs4925114)

# Concatenate the SNPs into a single variable
rs2282180=paste(canhr$rs2282180_A1, canhr$rs2282180_A2, sep="");table(rs2282180) #G>A
# Condense the variable down into a single catagorical variable
canhr$rs2282180[rs2282180==00]=NA;canhr$rs2282180[rs2282180==33]=1
canhr$rs2282180[rs2282180==13]=2;canhr$rs2282180[rs2282180==11]=2;table(canhr$rs2282180)

# Rename data.frame 
data1=canhr;names(data1)

# Make SNPs numeric variables
data1$rs599961=as.numeric(data1$rs599961); table(data1$rs599961)
data1$rs1502593=as.numeric(data1$rs1502593);table(data1$rs1502593)
data1$rs3071=as.numeric(data1$rs3071);table(data1$rs3071)
data1$rs7849=as.numeric(data1$rs7849);table(data1$rs7849)
data1$rs11190480=as.numeric(data1$rs11190480);table(data1$rs11190480)
data1$rs2167444=as.numeric(data1$rs2167444);table(data1$rs2167444)
data1$rs3829160=as.numeric(data1$rs3829160);table(data1$rs3829160)
data1$rs2234970=as.numeric(data1$rs2234970);table(data1$rs2234970)
data1$rs11557927=as.numeric(data1$rs11557927);table(data1$rs11557927)
data1$rs3978768=as.numeric(data1$rs3978768);table(data1$rs3978768)
data1$rs5435=as.numeric(data1$rs5435);table(data1$rs5435)
data1$rs5415=as.numeric(data1$rs5415);table(data1$rs5415)
data1$rs4925114=as.numeric(data1$rs4925114);table(data1$rs4925114)
data1$rs2282180=as.numeric(data1$rs2282180);table(data1$rs2282180)

# Create new variable that is catagorical age
data1$Age_cat=data1$A1+data1$A2+data1$A3
data1$Age_cat=as.factor(data1$Age_cat)
data1$Age_cat=factor(data1$Age_cat, labels=1:4)

# Create new variable that is catagorical d15n
data1$D15N_cat=data1$DN1+data1$DN2+data1$DN3
as.factor(data1$D15N_cat)
data1$D15N_cat=factor(data1$D15N_cat, labels=1:4)

# VillageGroup as factor
data1$VillageGroup=as.factor(data1$VillageGroup)

# uniqueid as factor
data1$uniqueid=paste(data1$Fam, data1$Ind, sep=".")
length(unique(data1$uniqueid))
data1$d15n2=data1$d15n*data1$d15n

# **************************************************************************** #
# ***************               Create Genotyped variable      *************** #
# **************************************************************************** # 

# PUFA SNPs Gene
geno_sum= data1$rs1502593_A1+data1$rs1502593_A2+        #rs1502593  # SCD
  data1$rs522951_A1+data1$rs522951_A2+          #rs522951
  data1$rs11190480_A1+data1$rs11190480_A2+      #rs11190480
  data1$rs3071_A1+data1$rs3071_A2+              #rs3071
  data1$rs3829160_A1+data1$rs3829160_A2+        #rs3829160
  data1$rs2234970_A1+data1$rs2234970_A2+        #rs2234970
  data1$rs599961_A1+data1$rs599961_A2+          #rs599961
  data1$rs41290540_A1+data1$rs41290540_A2+      #rs41290540
  data1$rs3978768_A1+data1$rs3978768_A2+        #rs3978768
  data1$rs11557927_A1+data1$rs11557927_A2+      #rs11557927
  data1$rs7849_A1+data1$rs7849_A2+              #rs7849
  data1$rs2167444_A1+data1$rs2167444_A2+        #rs2167444
  data1$rs2654185_A1+data1$rs2654185_A2+        #rs2654185  # SLC2A4
  data1$rs5415_A1+data1$rs5415_A2+              #rs5415
  data1$rs5417_A1+data1$rs5417_A2+              #rs5417
  data1$rs16956647_A1+data1$rs16956647_A2+      #rs16956647
  data1$rs5435_A1+data1$rs5435_A2+              #rs5435
  data1$rs3744405_A1+data1$rs3744405_A2+        #rs3744405
  data1$rs2297508_A1+data1$rs2297508_A2+        #rs2297508  # SREBF1
  data1$rs2282180_A1+data1$rs2282180_A2+        #rs16956647
  data1$rs9899634_A1+data1$rs9899634_A2+        #rs9899634
  data1$rs8066560_A1+data1$rs8066560_A2+        #rs8066560
  data1$rs9902941_A1+data1$rs9902941_A2        #rs9902941

# **************************************************************************** #
# ***************               Look at Covariates             *************** #
# **************************************************************************** # 

# who was genotyped?
data1$geno=ifelse(geno_sum>0,1,0)
table(data1$geno) #1138

# who has canhr_id? 
data1$canhr_id
data1$geno.c=ifelse(is.na(data1$canhr_id),0,1)
table(data1$geno.c)  #1259

# who has canhr_id and is genotyped?
geno_canhr=ifelse(data1$geno.c==1 & data1$geno==1,1,0)
table(geno_canhr)  #1138
names(data1)

# who has d15n measure?
data1$d15n
data1$PUFA=ifelse(is.na(data1$d15n),0,1)
table(data1$PUFA) #1141

# who has d15n and genotyped
geno_PUFA=ifelse(data1$geno==1 & data1$PUFA==1,1,0)
table(geno_PUFA) #1135
data1$geno_PUFA=geno_PUFA

# who has village variable?
village=ifelse(is.na(data1$VillageGroup),0,1)
table(village) #1144
#who has village and PUFA and Geno
geno_PUFA_village=ifelse(data1$geno==1 & data1$PUFA==1 & village==1,1,0)
table(geno_PUFA_village) #1135

# Subset to only include individuals in the analysis
data_sub=subset(data1,geno_PUFA_village==1)
  # How many participants have medications?
    table(data_sub$meds_chol)
  # Whata re the CANHR-ID for these participants?
    #write.table(data_sub$canhr_id,file="CANHR_ID_PUFA_Genes_21Oct13.csv",sep="",row.names=F)

# subset to include only individuals with Geno & PUFA
#---------------------------------------------------
  newdata <- subset(data1, data1$geno_PUFA==1) 
    # How many participants have lipid medications?        
      table(newdata$meds_chol) #55  
    # How many participants have T2D drugs?
      table(newdata$meds_diab) #13
    # what is gender break out?
      table(newdata$Sex) # 538 men/597 women
    # what is d15n (mean/sd) by PUFA quartile?
      newdata$D15N_cat=as.numeric(newdata$D15N_cat)
      ddply(newdata,.(D15N_cat), summarize, mean = mean(d15n),sd = sd(d15n),len=length(d15n),max=max(d15n,na.rm=T),min=min(d15n,na.rm=T) )
      newdata$D15N_cat=as.numeric(newdata$D15N_cat)
      #D15N_cat      mean        sd       len
        #1            7.327860 0.3244182 272
        #2            8.170209 0.2360636 278
        #3            9.127514 0.3291831 288
        #4            11.030758 1.0540485 297
    # what is total mean/sd?
      mean(newdata$d15n) #8.95
      sd(newdata$d15n) # 1.5
      range(newdata$d15n, na.rm=T) #6.2-15.2
    # what is d15n by sex
      ddply(newdata,.(Sex), summarize, mean = mean(d15n),sd = sd(d15n),len=length(d15n),max=max(d15n,na.rm=T),min=min(d15n,na.rm=T) )
        #Sex     mean       sd len
        #8.773188 1.497302 538 (men)
        #9.127940 1.501683 597 (women)
# **************************************************************************** #
# ***************                          Analysis            *************** #
# **************************************************************************** #

# uniqueid as factor
data1$uniqueid=paste(data1$Fam, data1$Ind, sep=".");length(unique(data1$uniqueid))

# WHR=rs599961*d15n
  lmekin=lmekin(WHR ~ age + Sex + VillageGroup + d15n + rs599961 + rs599961*d15n +(1|uniqueid), data=data1, na.action=na.exclude)
  lmekin;summary(lmekin);
  # Interaction model is significant: time to plot!
    table(data$rs599961) # Resulting data.sets should be same as genotype breakdown
    rs599961_0 <- subset(data1, rs599961=="1");rs599961_1 <- subset(data1, rs599961=="2")
    rs599961_2 <- subset(data1, rs599961=="3")
  # construct model 1 from subset "0" data.set
    reg0 <- lmekin(WHR ~ age + Sex + VillageGroup + d15n + (1|uniqueid), data=rs599961_0, na.action=na.exclude); reg0
    # extract intercept and slope from model 1
      int_0=reg0$coefficients$fixed[1]; slope_0=reg0$coefficients$fixed[5]
  # construct model 2 from subset "1" data.set
      reg1 <- lmekin(WHR ~ age + Sex + VillageGroup + d15n + (1|uniqueid), data=rs599961_1, na.action=na.exclude); reg1
  # extract intercept and slope from model 2
      int_1=reg1$coefficients$fixed[1]; slope_1=reg1$coefficients$fixed[5]
  # construct model 3 from subset "2" data.set
      reg2 <- lmekin(WHR ~ age + Sex + VillageGroup + d15n + (1|uniqueid), data=rs599961_2, na.action=na.exclude); reg2
  # extract intercept and slope from model 2
      int_2=reg2$coefficients$fixed[1]; slope_2=reg2$coefficients$fixed[5]
  # plot regression lines
    range(data1$d15n,na.rm=T); range(data1$WHR,na.rm=T) # c(1.4,1.55) (WHR)
    plot(WHR ~ age + Sex + VillageGroup + rs599961 + d15n,
     xlim=range(data1$d15n,na.rm=T),ylim=c(1.4,1.55), xlab="",ylab="",
     data=data1, type='n'); 
    mtext('Waist-to-Hip Ratio', 2, line=2.3, cex = 1.2)
    mtext(expression(paste('n-3 PUFA intake (',{delta}^15,'N)')), 1, line=2.6, cex = 1.2)
    abline(int_0,slope_0,lty=1,col="red");
    abline(int_1,slope_1,lty=2,col="black");
    abline(int_2,slope_2,lty=5,col="blue")
    legend("center", c("T/T (n=394)","T/G (n=512)","G/G (n=166)"),title ="Genotype [T>G]", lty=c(1,2,6),col=c("red","black","blue"),cex=0.9)

############

# WHR=rs7849*PUFA
  lmekin=lmekin(WHR ~ age + Sex + VillageGroup + d15n + rs7849 + rs7849*d15n + (1|uniqueid), data=data1, na.action=na.exclude)
    lmekin
  # Interaction model is significant: time to plot!
    table(data1$rs7849);rs7849_0 <- subset(data1, rs7849=="1");rs7849_1 <- subset(data1, rs7849=="2")
    rs7849_2 <- subset(data1, rs7849=="3")
  # construct model 1 from subset "0" data.set
    reg0 <- lmekin(WHR ~ age + Sex + VillageGroup + d15n + (1|uniqueid), data=rs7849_0, na.action=na.exclude); reg0
  # extract intercept and slope from model 1
    int_0=reg0$coefficients$fixed[1]; slope_0=reg0$coefficients$fixed[5]
  # construct model 2 from subset "1" data.set
    reg1 <- lmekin(WHR ~ age + Sex + VillageGroup + d15n + (1|uniqueid), data=rs7849_1, na.action=na.exclude); reg1
  # extract intercept and slope from model 2
    int_1=reg1$coefficients$fixed[1]; slope_1=reg1$coefficients$fixed[5]
  # construct model 3 from subset "2" data.set
    reg2 <- lmekin(WHR ~ age + Sex + VillageGroup + d15n + (1|uniqueid), data=rs7849_2, na.action=na.exclude); reg2
  # extract intercept and slope from model 2
    int_2=reg2$coefficients$fixed[1]; slope_2=reg2$coefficients$fixed[5]
  # plot regression lines
    range(data1$d15n,na.rm=T);range(data1$WHR,na.rm=T) # c(1.4,1.55) (WHR)
    plot(WHR ~ age + Sex + VillageGroup + rs7849 + d15n,
     xlim=range(data1$d15n,na.rm=T),ylim=c(1.4,1.6), xlab="",ylab="",
     data=data1, type='n')
    mtext('Waist-to-Hip Ratio', 2, line=2.3, cex = 1.2)
    mtext(expression(paste('n-3 PUFA intake (',{delta}^15,'N)')), 1, line=2.6, cex = 1.2)

    abline(int_0,slope_0,lty=1,col="red");abline(int_1,slope_1,lty=2,col="black")
    abline(int_2,slope_2,lty=5,col="blue")
    legend("center", c("T/T (n=279)","T/C (n=557)","C/C (n=238)"),title ="Genotype [T>C]", lty=c(1,2,6),col=c("red","black","blue"))

############
data1$d15n2=data1$d15n*data1$d15n

# HOMA-IR=rs2282180*PUFA
  lmekin=lmekin(HOMAIR ~ age + Sex + VillageGroup + d15n + rs2282180 + rs2282180*d15n + (1|uniqueid), data=data1, na.action=na.exclude)
    lmekin
  # Interaction model is significant: time to plot!
    table(data1$rs2282180);rs2282180_0 <- subset(data1, rs2282180=="1");rs2282180_1 <- subset(data1, rs2282180=="2")
    #rs2282180_2 <- subset(data1, rs2282180=="3")
  # construct model 1 from subset "0" data.set
    reg0 <- lmekin(HOMAIR ~ age + Sex + VillageGroup + d15n + (1|uniqueid), data=rs2282180_0, na.action=na.exclude); reg0
  # extract intercept and slope from model 1
    int_0=reg0$coefficients$fixed[1]; slope_0=reg0$coefficients$fixed[5]
  # construct model 2 from subset "1" data.set
    reg1 <- lmekin(HOMAIR ~ age + Sex + VillageGroup + d15n + (1|uniqueid), data=rs2282180_1, na.action=na.exclude); reg1
  # extract intercept and slope from model 2
    int_1=reg1$coefficients$fixed[1]; slope_1=reg1$coefficients$fixed[5]

    range(data1$d15n,na.rm=T);range(data1$HOMAIR,na.rm=T) # c(1.4,1.55) (WHR)
  plot(HOMAIR ~ age + Sex + VillageGroup + rs2282180 + d15n,
     xlim=range(data1$d15n,na.rm=T),ylim=c(0,10), xlab="",ylab="",
     data=data1, type='n')
    mtext('HOMA-IR', 2, line=2.3, cex = 1.2)
    mtext(expression(paste('n-3 PUFA intake (',{delta}^15,'N)')), 1, line=2.6, cex = 1.2)
  abline(int_0,slope_0,lty=1,col="red");abline(int_1,slope_1,lty=2,col="black")
  
  legend("center", c("G/G (n=963)", "G/A (n=109) \n A/A"),title ="Genotype [G>A]", lty=c(1,2),col=c("red","black"),cex=0.9)

#######

  # CHOL=rs11190480*PUFA
    lmekin=lmekin(cholesterol ~ age + Sex + VillageGroup + d15n + rs11190480 + rs11190480*d15n + (1|uniqueid), data=data1, na.action=na.exclude)
    lmekin
  # Interaction model is significant: time to plot!
    table(data1$rs11190480);rs11190480_0 <- subset(data1, rs11190480=="1");rs11190480_1 <- subset(data1, rs11190480=="2")
    rs11190480_2 <- subset(data1, rs11190480=="3")
  # construct model 1 from subset "0" data.set
    reg0 <- lmekin(cholesterol ~ age + Sex + VillageGroup + d15n + (1|uniqueid), data=rs11190480_0, na.action=na.exclude); reg0
  # extract intercept and slope from model 1
    int_0=reg0$coefficients$fixed[1]; slope_0=reg0$coefficients$fixed[5]
  # construct model 2 from subset "1" data.set
    reg1 <- lmekin(cholesterol ~ age + Sex + VillageGroup + d15n + (1|uniqueid), data=rs11190480_1, na.action=na.exclude); reg1
  # extract intercept and slope from model 2
    int_1=reg1$coefficients$fixed[1]; slope_1=reg1$coefficients$fixed[5]
  # construct model 3 from subset "2" data.set
    reg2 <- lmekin(cholesterol ~ age + Sex + VillageGroup + d15n + (1|uniqueid), data=rs11190480_2, na.action=na.exclude); reg2
  # extract intercept and slope from model 2
    int_2=reg2$coefficients$fixed[1]; slope_2=reg2$coefficients$fixed[5]
  # plot regression lines
    range(data1$d15n,na.rm=T);range(data1$cholesterol,na.rm=T) # c(1.4,1.55) (WHR)
    plot(Cholp ~ age + Sex + VillageGroup + rs11190480 + d15n,
     xlim=range(data1$d15n,na.rm=T),ylim=c(100,200), xlab="",ylab="",
     data=data1, type='n')
    mtext('Cholesterol (mg/dL)', 2, line=2.3, cex = 1.2)
    mtext(expression(paste('n-3 PUFA intake (',{delta}^15,'N)')), 1, line=2.6, cex = 1.2)
    abline(int_0,slope_0,lty=1,col="red");abline(int_1,slope_1,lty=2,col="black")
    abline(int_2,slope_2,lty=5,col="blue")

    # output figure and legend seperately.

    # need to figure out how to get plot to right.
    legend("center", c("A/A (n=871)","A/G (n=241)","G/G (n=16)"),title ="Genotype [A>G]", lty=c(1,2,6),col=c("red","black","blue"),cex=0.9)
