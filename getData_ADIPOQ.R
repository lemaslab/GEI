
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
  mutate(unique_id=paste0(Fam,".",Ind),
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
  select(unique_id,Fam,Ind,FA,MO,Sex,d15n,Genotyped,
         snp1.rs10865710, snp2.rs12497191, snp3.rs1801282, 
         snp4.rs3856806, snp5.rs10937273, snp6.rs822387,
         snp7.APM11426, snp8.rs17300539,  snp9.rs266729,
         snp10.rs182052, snp11.rs822393, snp12.rs822394,
         snp13.rs822395, snp14.rs822396, snp15.rs17846866, 
         snp16.SNP2, snp17.rs2241766, snp18.rs1501299,
         snp19.rs2241767, snp20.rs3774261, snp21.rs3774262,
         snp22.rs35554619, snp23.rs8192678, snp24.rs17574213,
         snp25.rs2970847, snp26.rs135549, snp27.rs135539,
         snp28.rs1800206, snp29.rs4253778)%>%
  gather(snp, allele, snp1.rs10865710:snp29.rs4253778)%>%
  separate(allele, c("A1", "A2"), sep="_")%>%
  mutate(A1=recode(A1, # A=1, C=2, G=3, T=4
                      "1"="A",
                      "2"="C", 
                      "3"="G",
                      "4"="T",
                      "0"="NA"),
         A2=recode(A2, # A=1, C=2, G=3, T=4
                   "1"="A",
                   "2"="C", 
                   "3"="G",
                   "4"="T",
                   "0"="NA"))%>%na_if("NA")  %>%
  # spread(snp, allele)

names(df1)

geno=df1%>%
  select(unique_id,snp1.rs10865710,snp2.rs12497191,snp3.rs1801282)%>%
  as.matrix()

geno.ready <- raw.data(data =geno, frame = "wide") 
                       base = TRUE, sweep.sample = 0.5, 
                       call.rate = 0.95, maf = 0.10, 
                       imput = FALSE)



names(df1)
head(df1)
str(df1)
dim(df1)
df1$snp1.rs10865710
length(unique(df1$unique_id))
unique(df1$allele)
table(df1$allele)

# ready to recode
df2=df1%>%
  select(unique_id,snp1.rs10865710, snp2.rs12497191)
  
  geno.readySTR <- raw.data(data = as.matrix(geno), frame = "long", base = TRUE, sweep.sample = 0.5, call.rate = 0.95, maf = 0.10, imput = FALSE, outfile = "structure")
Mstr <- geno.readySTR$M.clean
  
  

# what is the ordering of redcap_events
levels(dat$clinic_visit)
unique(as.character(dat$Aliquot.Type))
unique(as.character(dat$tube.type))
table(dat$tube.type)
table(dat$Freezer.Section)

# look for "2_weeks" vs "2_week"
which(dat$clinic_visit=="2_weeks") 
change=dat%>%
  filter(clinic_visit=="2_weeks")
  print(change$Participant_ID)

# which "3rd Trimester" vs "3rd_trimester"  # BLS059A needs to be changed
  change=dat%>%
    filter(clinic_visit=="3rd Trimester")
  print(change$Participant_ID)
  
# set the order of redcap_events
df <- dat %>% 
  mutate(clinic_visit = factor(clinic_visit, 
                                    levels = c("3rd_trimester", 
                                               "2_week", 
                                               "2_months",
                                               "6_months",
                                               "12_months")))
# Factored Variables: 
# clinic visit
levels(df$clinic_visit)
table(df$clinic_visit)

# tube type
levels(df$tube.type)
table(df$tube.type)

# freezer section
levels(dat$Freezer.Section)
table(dat$Freezer.Section)

# change dates
df$Clinic.visit.date=as.Date(df$Clinic.visit.date, "%m/%d/%Y")
dim(df) # 2021
length(unique(df$Participant_ID)) #85

names(df)
df$Mom_Baby

# create paired mom-baby
df$part_id_link=gsub("A","",df$Participant_ID)
df$part_id_link=gsub("B$","",df$part_id_link)
 
# drop NA observations
dat.s=df %>%
  group_by(Participant_ID, clinic_visit) %>%
  arrange(Clinic.visit.date) 
  dim(dat.s) # 2021
  length(unique(dat.s$Participant_ID)) #85
  names(dat.s)
  
# how many visits
  dat.s%>%
    group_by(clinic_visit)%>%
    summarize(count=n_distinct(part_id_link))
  table(dat.s$clinic_visit)

# how many tubes per participant?
part_count=dat.s %>%
  group_by(part_id_link) %>%
  summarize(count=n_distinct(crc_specimen_barcode))
  mean(part_count$count) # 43.9 tubes
            
# how many sample types
dat.s %>%
  group_by(Aliquot.Type) %>%
  summarize(count=n_distinct(crc_specimen_barcode))

# how many tubes per sample type
dat.s %>%
  group_by(Freezer.Section) %>%
  summarize(count=n_distinct(crc_specimen_barcode))

# what about those that have completed 12-month visit
table(dat.s$clinic_visit)  #43 sample
year.complete=dat.s %>%
  filter(any(clinic_visit %in% "12_months"))
dim(year.complete)
length(unique(year.complete$Participant_ID)) # 9 participants (mom-baby)
length(unique(year.complete$part_id_link)) # 5 participants (mom-baby)
# why is there a missing mom-baby?
year.all=dat.s%>%
filter(part_id_link%in%c("BLS001","BLS002","BLS003","BLS011","BLS016"))

# how many tubes per sample type- AMONG - people completed 12-month visit
d2=year.all %>%
  group_by(part_id_link,Freezer.Section)%>%
  summarize(count=n_distinct(crc_specimen_barcode))

d2%>%
  group_by(Freezer.Section)%>%
  summarize(mean=mean(count),
            min=min(count),
            max=max(count))

# output data for import to redcap
redcap=dat.s %>%
  select(Participant_ID,clinic_visit,Clinic.visit.date,Mom_Baby,Aliquot.Type,crc_specimen_barcode,crc_specimen_number,
         tube.type,Aliquot.Number) %>%
  rename(test_id=Participant_ID, 
         biosample_study_visit=clinic_visit,
         biosample_collection_date=Clinic.visit.date,
         biosample_mom_baby=Mom_Baby,
         biosample_aliquot_type=Aliquot.Type,
         crc_specimen_barcode=crc_specimen_barcode,
         crc_specimen_number=crc_specimen_number,
         biosample_tube_type=tube.type,
         biosample_aliquot_numb=Aliquot.Number)%>%
  mutate(redcap_event_name=NA,
         redcap_repeat_instrument="biological_specimen_collection")%>%
  arrange(test_id,biosample_study_visit,biosample_collection_date,biosample_aliquot_type)%>%
  #group_by(test_id,biosample_study_visit) %>% 
  mutate(redcap_repeat_instance=row_number())%>%
  select(test_id,redcap_event_name,redcap_repeat_instrument,redcap_repeat_instance,everything())%>%
  mutate(redcap_event_name=case_when(biosample_study_visit=="3rd_trimester" ~ "third_trimester_arm_1",
                                     biosample_study_visit=="2_week" ~ "two_week_arm_1",
                                     biosample_study_visit=="2_months" ~ "two_month_arm_1",
                                     biosample_study_visit=="6_months" ~ "six_month_arm_1",
                                     biosample_study_visit=="12_months" ~ "twelve_month_arm_1"))%>%
  mutate(redcap_event_name = factor(redcap_event_name, 
                               levels = c("third_trimester_arm_1", 
                                          "two_week_arm_1", 
                                          "two_month_arm_1",
                                          "six_month_arm_1",
                                          "twelve_month_arm_1")))%>%
  mutate(biosample_collection_date=format(biosample_collection_date, "%m/%d/%Y"))%>%
  ungroup()%>%
  mutate(biosample_study_visit=recode(biosample_study_visit, 
                     "3rd_trimester"="1", 
                     "2_week"="2",
                     "2_months"="3",
                     "6_months"="4",
                     "12_months"="5"),
         biosample_mom_baby=recode(biosample_mom_baby,
                      "mom"="0",
                      "baby"="1"),
         biosample_aliquot_type=recode(biosample_aliquot_type,
                      "plasma"="1",
                      "urine"="2",
                      "saliva"="3",
                      "milk- skim"="4",
                      "milk- whole"="5",
                      "milk-lipid"="6",
                      "stool"="7",
                      "vaginal"="8",
                      "blood"="9",
                      "formula"="10"),
         biosample_tube_type=recode(biosample_tube_type,
                      "2ml"="1",
                      "ez sample"="2",
                      "vaginal vial"="3",
                      "5ml"="4",
                      "tiny"="5",
                      "blood card"="6",
                      "other"="7",
                      "15ml"="8",
                      "saliva tube"="9",
                      "50ml"="10"))
# check data
names(redcap)
dim(redcap)

# replace NA with blanks
df <- sapply(redcap, as.character)
df[is.na(df)] <- " "
df1=as.data.frame(df)

# checks
unique(df1$redcap_event_name)
table(df1$redcap_event_name)
table(df1$biosample_study_visit)

# export test data: BLS001A
redcap.bls001=df1%>%
  filter(test_id=="BLS001A")%>%
write_csv(path =paste0(work.dir,"redcap.bls001.csv",na = ""))

# 

# need to create a report with data that needs to be followed up.
# output as html.

# redcap_event_name
# crc_specimen_number
# crc_specimen_barcode
# biosample_collection_date


 


