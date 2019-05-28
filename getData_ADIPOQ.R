
#' ---
#' title: "Data management for ADIPOQ SNPS"
#' author: "Dominick Lemas"
#' date: "May 22, 2019"
#' ---


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

# expect monomorphic among (ADIPOQ:rs822387, ADIPOQ:rs822395, ADIPOQ:rs3774262, 
#                           ADIPOQ:rs35554619, PPARGC1A:rs8192678)

df=dat%>%
  select(Fam,Ind,FA,MO,Sex,d15n,Genotyped,
         WHRp,BMIp,BodFp,VillageGroup,d15n,age,Adipp, 
         rs10865710_A1, rs10865710_A2,        # PPARG
         rs12497191_A1, rs12497191_A2,
         rs1801282_A1, rs1801282_A2,
         rs3856806_A1, rs3856806_A2,
         rs10937273_A1, rs10937273_A2,   # ADIPOQ
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
         rs8192678_A1, rs8192678_A2,     # PPARGC1A
         rs17574213_A1, rs17574213_A2,
         rs2970847_A1, rs2970847_A2,
         rs135549_A1, rs135549_A2,       # PPARA
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
  select(unique_id,unique_fa,unique_mo,Fam,Ind,FA,MO,Sex,d15n,Genotyped,
         WHRp,BMIp,BodFp,VillageGroup,d15n,age,Adipp,
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
  select(unique_id,unique_fa,unique_mo, Fam,Ind,FA,MO,Sex,d15n,Genotyped,WHRp,
        BMIp,BodFp,VillageGroup,age,Adipp)

# merge
df.analysis=left_join(df.m, dat, by="unique_id")
names(df.analysis)
range(df.analysis$Fam)
df.analysis$unique_id

# test analysis 
data(minnbreast)
str(minnbreast)
minnbreast$id

str(df.analysis)
unique(df.analysis$unique_id)
which(duplicated(df.analysis$unique_id))
which(duplicated(df.analysis$Ind))
df.analysis$Sex=as.character(df.analysis$Sex)
df.analysis$unique_id=as.numeric(df.analysis$unique_id)

gped <- with(df.analysis, pedigree(unique_id, unique_fa, unique_mo, Sex, Fam))

# http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.592.1810&rep=rep1&type=pdf



fit2 <- lmekin(Adipp~snp7.APM11426+Sex+VillageGroup+age+d15n+(1|unique_id), data=df.analysis, method="ML")
summary(fit2)
names(df.analysis)
tidy(fit2)
