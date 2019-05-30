
#' ---
#' title: "Data management for ADIPOQ SNPS"
#' author: "Dominick Lemas"
#' date: "May 22, 2019"
#' ---

# **************************************************************************** #
# ***************                Library                       *************** #
# **************************************************************************** #

library(tidyr)
library(tidyverse)
library(dplyr)
library(snpReady)
library(naniar)
library(kinship2)
library(coxme)
library(nlme)
library(stargazer)
library(effects)
library(ggplot2)

library(sjPlot)
library(sjmisc)
library(ggplot2)
data(efc)

# **************************************************************************** #
# ***************                Directory Variables           *************** #
# **************************************************************************** #

work.dir=paste0(Sys.getenv("USERPROFILE"),"\\Dropbox (UFL)\\02_Projects\\CANHR\\data\\");work.dir
data.dir=paste0(Sys.getenv("USERPROFILE"),"\\Dropbox (UFL)\\02_Projects\\CANHR\\data\\");data.dir
result.dir=paste0(Sys.getenv("USERPROFILE"),"\\Dropbox (UFL)\\02_Projects\\CANHR\\results\\");result.dir

# Set Working Directory
setwd(work.dir)
list.files()

# **************************************************************************** #
# ***************         AssocPed_Full3_10Oct13.dat                                              
# **************************************************************************** # 

# Read Data
data.file.name="AssocPed_Full3_10Oct13.dat";data.file.name
data.file.path=paste0(data.dir,data.file.name);data.file.path
canhr<- read_tsv(data.file.path);

# look at data
dat=canhr
head(dat); str(dat); names(dat)

# **************************************************************************** #
# ***************  Select ADIPOQ SNPs                                             
# **************************************************************************** #

names(dat)

df=dat%>%
  select(Fam,Ind,FA,MO,Sex,d15n,Genotyped,    # covariates
         VillageGroup,d15n,age,
         Adipp,ApoA1p,BMIp,BodFp,HCAp,HDLp,   # outcomes
         LCholp,LDLp,LTrigp,ThCAp,VLDLp,WCAp,
         rs10865710_A1, rs10865710_A2,        # PPARG
         rs12497191_A1, rs12497191_A2,
         rs1801282_A1, rs1801282_A2,
         rs3856806_A1, rs3856806_A2,
         rs10937273_A1, rs10937273_A2,        # ADIPOQ
         rs822387_A1, rs822387_A2,
         rs822387_A1, rs822387_A2,
         APM11426_A1, APM11426_A2,
         rs17300539_A1, rs17300539_A2,
         rs266729_A1, rs266729_A2,
         rs182052_A1, rs182052_A2,
         rs822393_A1, rs822393_A2,
         rs822394_A1, rs822394_A2,
         rs822395_A1, rs822395_A2,
         rs822396_A1, rs822396_A2,
         rs17846866_A1, rs17846866_A2,
         Pasker_SNP2_A1, Pasker_SNP2_A2,
         rs2241766_A1, rs2241766_A2,
         rs1501299_A1, rs1501299_A2,
         rs2241767_A1, rs2241767_A2,
         rs3774261_A1, rs3774261_A2,
         rs3774262_A1, rs3774262_A2,
         rs35554619_A1, rs35554619_A2,
         rs8192678_A1, rs8192678_A2,           # PPARGC1A
         rs17574213_A1, rs17574213_A2,
         rs2970847_A1, rs2970847_A2,
         rs135549_A1, rs135549_A2,             # PPARA
         rs135539_A1, rs135539_A2, 
         rs1800206_A1, rs1800206_A2,
         rs4253778_A1, rs4253778_A2)

# check
names(df)

# merge SNP alleles
df1=df%>%
  mutate(unique_id=as.character(paste0(Fam,".",Ind)),
         unique_fa=as.character(paste0(Fam,".",FA)),
         unique_mo=as.character(paste0(Fam,".",MO)),
         snp1.rs10865710=paste0(rs10865710_A1,"_",rs10865710_A2),
         snp2.rs12497191=paste0(rs12497191_A1,"_",rs12497191_A2),
         snp3.rs1801282=paste0(rs1801282_A1,"_",rs1801282_A2),
         snp4.rs3856806=paste0(rs3856806_A1,"_",rs3856806_A2),
         snp5.rs10937273=paste0(rs10937273_A1,"_",rs10937273_A2),
         snp6.rs822387=paste0(rs822387_A1,"_",rs822387_A2),
         snp7.APM11426=paste0(APM11426_A1,"_",APM11426_A2),
         snp8.rs17300539=paste0(rs17300539_A1,"_",rs17300539_A2),
         snp9.rs266729=paste0(rs266729_A1,"_",rs266729_A2),
         snp10.rs182052=paste0(rs182052_A1,"_",rs182052_A2),
         snp11.rs822393=paste0(rs822393_A1,"_",rs822393_A2),
         snp12.rs822394=paste0(rs822394_A1,"_",rs822394_A2),
         snp13.rs822395=paste0(rs822395_A1,"_",rs822395_A2),
         snp14.rs822396=paste0(rs822396_A1,"_",rs822396_A2),
         snp15.rs17846866=paste0(rs17846866_A1,"_",rs17846866_A2),
         snp16.SNP2=paste0(Pasker_SNP2_A1,"_",Pasker_SNP2_A2),
         snp17.rs2241766=paste0(rs2241766_A1,"_",rs2241766_A2),
         snp18.rs1501299=paste0(rs1501299_A1,"_",rs1501299_A2),
         snp19.rs2241767=paste0(rs2241767_A1,"_",rs2241767_A2),
         snp20.rs3774261=paste0(rs3774261_A1,"_",rs3774261_A2),
         snp21.rs3774262=paste0(rs3774262_A1,"_",rs3774262_A2),
         snp22.rs35554619=paste0(rs35554619_A1,"_",rs35554619_A2),
         snp23.rs8192678=paste0(rs8192678_A1,"_",rs8192678_A2),
         snp24.rs17574213=paste0(rs17574213_A1,"_",rs17574213_A2),
         snp25.rs2970847=paste0(rs2970847_A1,"_",rs2970847_A2),
         snp26.rs135549=paste0(rs135549_A1,"_",rs135549_A2),
         snp27.rs135539=paste0(rs135539_A1,"_",rs135539_A2),
         snp28.rs1800206=paste0(rs1800206_A1,"_",rs1800206_A2),
         snp29.rs4253778=paste0(rs4253778_A1,"_",rs4253778_A2))%>%
  select(unique_id,unique_fa,unique_mo,Fam,Ind,FA,MO,Sex,d15n,Genotyped,VillageGroup,age,
         Adipp,ApoA1p,BMIp,BodFp,HCAp,HDLp,   # outcomes
         LCholp,LDLp,LTrigp,ThCAp,VLDLp,WCAp, # outcomes
         snp1.rs10865710, snp2.rs12497191, snp3.rs1801282, 
         snp4.rs3856806, snp5.rs10937273, snp6.rs822387,
         snp7.APM11426, snp8.rs17300539,  snp9.rs266729,
         snp10.rs182052, snp11.rs822393, snp12.rs822394,
         snp13.rs822395, snp14.rs822396, snp15.rs17846866, 
         snp16.SNP2, snp17.rs2241766, snp18.rs1501299,
         snp19.rs2241767, snp20.rs3774261, snp21.rs3774262,
         snp22.rs35554619, snp23.rs8192678, snp24.rs17574213,
         snp25.rs2970847, snp26.rs135549, snp27.rs135539,
         snp28.rs1800206, snp29.rs4253778)

df.geno=df1%>%
  gather(snp, allele, snp1.rs10865710:snp29.rs4253778)%>%
  separate(allele, c("A1", "A2"), sep="_")%>%
  mutate(A1=recode(A1, # A=1, C=2, G=3, T=4
                      "1"="A","2"="C", 
                      "3"="G","4"="T",
                      "0"="NA"),
         A2=recode(A2, # A=1, C=2, G=3, T=4
                   "1"="A","2"="C", 
                   "3"="G","4"="T",
                   "0"="NA"))%>%na_if("NA")%>%
  mutate(snp=factor(snp, levels = c("snp1.rs10865710","snp2.rs12497191","snp3.rs1801282", 
                                    "snp4.rs3856806","snp5.rs10937273","snp6.rs822387",
                                    "snp7.APM11426","snp8.rs17300539","snp9.rs266729",
                                    "snp10.rs182052","snp11.rs822393","snp12.rs822394",
                                    "snp13.rs822395","snp14.rs822396","snp15.rs17846866", 
                                    "snp16.SNP2","snp17.rs2241766","snp18.rs1501299",
                                    "snp19.rs2241767","snp20.rs3774261","snp21.rs3774262",
                                    "snp22.rs35554619","snp23.rs8192678","snp24.rs17574213",
                                    "snp25.rs2970847","snp26.rs135549","snp27.rs135539",
                                    "snp28.rs1800206","snp29.rs4253778")))
levels(df.geno$snp)
names(df.geno)

# **************************************************************************** #
# ***************                  SNP ready                                             
# **************************************************************************** #


# SNP ready
# https://cran.r-project.org/web/packages/snpReady/vignettes/snpReady-vignette.html

names(df.geno)
geno=df.geno%>%
  select(unique_id,snp,A1,A2)%>%
  as.matrix()

geno.ready <- raw.data(data =geno, frame = "long", 
                       base = TRUE, sweep.sample = 0.5, 
                       call.rate = 0.95, maf = 0.05, 
                       imput = FALSE)

Mwrth <- data.frame(geno.ready$M.clean)
Mwrth[1:10,1:5]

# $report$maf$r
# [1] "7 Markers removed by MAF = 0.05"
# 
# $report$maf$whichID
# [1] "snp16.SNP2"       "snp22.rs35554619" "snp25.rs2970847"  "snp28.rs1800206" 
# [5] "snp29.rs4253778"  "snp6.rs822387"    "snp8.rs17300539" 
# 

# $report$cr$r
# [1] "2 Markers removed by Call Rate = 0.95"
# 
# $report$cr$whichID
# [1] "snp22.rs35554619" "snp5.rs10937273" 

# **************************************************************************** #
# ***************                  popgen                                             
# **************************************************************************** #

pop.gen <- popgen(M=Mwrth)
head(pop.gen$whole$Markers)
qc_dat=pop.gen$whole$Markers
analysis_snps=row.names(qc_dat)
str(Mwrth)
summary(Mwrth)
names(Mwrth)

# convert row.names to column
dat=rownames_to_column(Mwrth,var="unique_id")
names(dat)
dim(dat)

genotypes=dat%>%
  as_tibble()%>%
  gather(snp, allele, snp1.rs10865710:snp9.rs266729)%>%
  group_by(snp,allele)%>%
  count()%>%
  ungroup() %>%
  spread(allele, n, fill=0)%>%
  rename("BB"=`0`,"AB"=`1`,"AA"=`2`,"Missing"=`<NA>`)%>%
  mutate(snp=factor(snp, levels = c("snp1.rs10865710","snp2.rs12497191","snp3.rs1801282", 
                                    "snp4.rs3856806","snp5.rs10937273","snp6.rs822387",
                                    "snp7.APM11426","snp8.rs17300539","snp9.rs266729",
                                    "snp10.rs182052","snp11.rs822393","snp12.rs822394",
                                    "snp13.rs822395","snp14.rs822396","snp15.rs17846866", 
                                    "snp16.SNP2","snp17.rs2241766","snp18.rs1501299",
                                    "snp19.rs2241767","snp20.rs3774261","snp21.rs3774262",
                                    "snp22.rs35554619","snp23.rs8192678","snp24.rs17574213",
                                    "snp25.rs2970847","snp26.rs135549","snp27.rs135539",
                                    "snp28.rs1800206","snp29.rs4253778")))%>%
  arrange(snp)%>%
  select(snp, AA, AB, BB, Missing)%>%
  write_csv(path =paste0(result.dir,"adipoq_genotypes.csv",na=""))

# merge genotype data back with outcomes data
# genotype data
dim(dat)
names(dat)
length(unique(dat$unique_id)) # 967

# outcomes data
dim(df1)
names(df1)
length(unique(df1$unique_id)) # 2231

df.m=df1%>%
  select(unique_id,unique_fa,unique_mo, Fam,Ind,FA,MO,Sex,d15n,Genotyped,VillageGroup,age,
         Adipp,ApoA1p,BMIp,BodFp,HCAp,HDLp,   # outcomes
         LCholp,LDLp,LTrigp,ThCAp,VLDLp,WCAp)

# merge
df.analysis=left_join(df.m, dat, by="unique_id")
names(df.analysis)
range(df.analysis$Fam)
df.analysis$unique_id

# extract function
extract_coxme_table <- function (mod){
  beta <- mod$coefficients$fixed
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- signif(1 - pchisq((beta/se)^2, 1), 2)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}

str(df.analysis)
length(unique(df.analysis$unique_id))

# outcome
outcome.index=c("Adipp","ApoA1p","BMIp","BodFp","HCAp","HDLp","LCholp",
                "LDLp","LTrigp","ThCAp","VLDLp","WCAp")
out.ind=length(outcome.index)

# index
snp.index = c("snp1.rs10865710","snp2.rs12497191","snp3.rs1801282", 
           "snp4.rs3856806", "snp7.APM11426","snp9.rs266729",
           "snp10.rs182052","snp11.rs822393","snp12.rs822394",
           "snp13.rs822395","snp14.rs822396","snp15.rs17846866", 
           "snp17.rs2241766","snp18.rs1501299",
           "snp19.rs2241767","snp21.rs3774262",
             "snp23.rs8192678",
           "snp26.rs135549","snp27.rs135539")

TABLE.1<-data.frame(outcome=character(),
                    snp=character(),
                     model=character(),
                     beta=numeric(),
                     se=numeric(),
                     z=numeric(),
                     p=numeric(),
                     stringsAsFactors=FALSE);TABLE.1

# Create index for loops
index=snp.index;index; 
myIndex<-length(index) 

# create table index
n=myIndex*out.ind

# create model index 
formula.snp=lapply(index, function(x){paste0("~",x,"+Sex+VillageGroup+age+d15n+(1|unique_id)")})
formula.outcome=lapply(outcome.index, function(x){paste0(x,formula.snp)})
formula.list=unlist(formula.outcome, recursive = TRUE)
formula.index=length(formula.list)


# Start the Loop
for (i in 1:(formula.index))
{
    # select model: kinda ugly
    formula.final=as.formula(formula.list[i])
    outcome=str_split(formula.final,"~")[[2]]
    snp.tmp=str_split(formula.final,"~")[[3]]
    snp=str_split_fixed(snp.tmp," ", n=2)[1]
    
    # run model
    fit <- lmekin(formula.final, data=df.analysis, method="ML")

    # extract params
    model.param=fit%>%extract_coxme_table()
    model=as_tibble(rownames_to_column(model.param,var="param"))

      # model output
      output=model%>%gather(key,value,beta:p)%>%
              mutate(outcome=outcome,
                     snp=snp,
                     model="main_effect")%>%
              select(outcome,snp,model, param, key, value)%>%
              filter(param==snp)%>%
              spread(key, value)%>%
              select(outcome,snp,model,beta,se,z,p)
  
      # create table index
      TABLE.1[i,1]=outcome
      TABLE.1[i,2]=output$snp
      TABLE.1[i,3]=output$model
      TABLE.1[i,4]=output$beta
      TABLE.1[i,5]=output$se
      TABLE.1[i,6]=output$z
      TABLE.1[i,7]=output$p
} # END

# modify table1 output
names(TABLE.1)
head(TABLE.1)
str(TABLE.1)

# round
table=TABLE.1%>%
  mutate(beta=round(beta, 1),
         se=round(se, 1),
         p=round(p, 3),
         beta_se=paste0("(β= ",beta," ± S.E.= ",se,")"),
         p_final=paste0(p," ",beta_se))%>%
  write_csv(path=paste0(result.dir,"adipoq_snps_main_effect.csv"),na="")

# interactions
TABLE.2<-data.frame(outcome=character(),
                    snp=character(),
                    model=character(),
                    beta=numeric(),
                    se=numeric(),
                    z=numeric(),
                    p=numeric(),
                    stringsAsFactors=FALSE);TABLE.2


# Create index for loops
index=snp.index;index; 
myIndex<-length(index) 

# create table index
n=myIndex*out.ind

# create model index 
formula.snp=lapply(index, function(x){paste0("~",x,"*d15n+Sex+VillageGroup+age+(1|unique_id)")})
formula.outcome=lapply(outcome.index, function(x){paste0(x,formula.snp)})
formula.list=unlist(formula.outcome, recursive = TRUE)
formula.index=length(formula.list)


# Start the Loop
for (i in 1:(formula.index))
{
  # select model: kinda ugly
  formula.final=as.formula(formula.list[i])
  outcome=str_split(formula.final,"~")[[2]]
  snp.tmp=str_split(formula.final,"~")[[3]]
  snp=str_split_fixed(snp.tmp," ", n=2)[1]
  
  # run model
  fit <- lmekin(formula.final, data=df.analysis, method="ML")
  
  # extract params
  model.param=fit%>%extract_coxme_table()
  model=as_tibble(rownames_to_column(model.param,var="param"))
  
  # model output
  output=model%>%gather(key,value,beta:p)%>%
    mutate(outcome=outcome,
           snp=snp,
           model="interaction")%>%
    select(outcome,snp,model, param, key, value)%>%
    filter(param==paste0(snp,":d15n"))%>%
    spread(key, value)%>%
    select(outcome,snp,model,beta,se,z,p)
  
  # create table index
  TABLE.2[i,1]=outcome
  TABLE.2[i,2]=output$snp
  TABLE.2[i,3]=output$model
  TABLE.2[i,4]=output$beta
  TABLE.2[i,5]=output$se
  TABLE.2[i,6]=output$z
  TABLE.2[i,7]=output$p
} # END

# modify table1 output
names(TABLE.2)
head(TABLE.2)
str(TABLE.2)
range(TABLE.2$p, na.rm=T)

# round
table=TABLE.2%>%
  mutate(beta=round(beta, 1),
         se=round(se, 1),
         p=round(p, 3),
         beta_se=paste0("(β= ",beta," ± S.E.= ",se,")"),
         p_final=paste0(p," ",beta_se))%>%
  write_csv(path=paste0(result.dir,"adipoq_snps_interaction.csv"),na="")

# explore interactions 
# https://ademos.people.uic.edu/Chapter13.html#3_continuous_x_categorical_regression
GPA.2.Model.1 <- lm(LDLp~snp11.rs822393+d15n+Sex+VillageGroup+age, df.analysis)
GPA.2.Model.2 <- lm(LDLp~snp11.rs822393*d15n+Sex+VillageGroup+age, df.analysis)

stargazer(GPA.2.Model.1, GPA.2.Model.2,type="html", 
          column.labels = c("Main Effects", "Interaction"), 
          intercept.bottom = FALSE, 
          single.row=FALSE,     
          notes.append = FALSE, 
          header=FALSE,
          out=paste0(result.dir,"LDLp_snp11.rs822393_d15n.html")) 


Inter.GPA.2 <- effect('snp11.rs822393*d15n', GPA.2.Model.2,
                      xlevels=list(genotype = c(-1, 0, 1)),
                      se=TRUE, confidence.level=.95, typical=mean)
Inter.GPA.2<-as.data.frame(Inter.GPA.2)

# interaction plot
# https://cran.r-project.org/web/packages/sjPlot/vignettes/plot_interactions.html

data(efc)
theme_set(theme_sjplot())

# make categorical
efc$c161sex <- to_factor(efc$c161sex)
efc$barthtot

# fit model with interaction
fit <- lm(neg_c_7 ~ c12hour + barthtot * c161sex, data = efc)

plot_model(fit, type = "pred", terms = c("barthtot", "c161sex"))


GPA.2.Model.2 <- lm(LDLp~snp11.rs822393*d15n+Sex+VillageGroup+age, df.analysis)
plot_model(GPA.2.Model.2, type = "pred", terms = c("snp11.rs822393", "d15n"))
plot_model(GPA.2.Model.2, type = "pred", terms = c("d15n", "snp11.rs822393"))



