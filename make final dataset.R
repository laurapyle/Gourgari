# make dataset for final analysis

library(lubridate)
library(dplyr)
library(tableone)
library(knitr)
library(data.table)
library(tidyr)
library(Hmisc)


####  PROTEOMICS DATA ######
filename <- "H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Proteomics data/experimentcodes.csv"
codes <- read.csv(filename)
filename <-  "H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Proteomics data//proteinquant.csv"
proteins <- read.csv(filename,header = FALSE)
# Transpose and format protein concentration data
proteins <- t(proteins)
proteins[1,1] <- "subjectid"
colnames(proteins) <- proteins[1,]
proteins <- proteins[-1,]
rownames(proteins) <- 1:nrow(proteins)
proteins <- as.data.frame(proteins)
# Get protein subjcet ids
proteins$subjectid <- substr(proteins$subjectid,7,8)
proteins[,2:ncol(proteins)] <- sapply(proteins[,2:ncol(proteins)],function(x) as.numeric(as.character(x)))
# Add subject id, date, etc.
groups <- codes[,c("MS.ID","ANALYTICID","COLLECTIONDT","SAMPLEGROUP")]
colnames(groups) <- c("subjectid","analyticid","date","group")
proteins <- merge(proteins,groups,by = "subjectid")
proteins$group <- ifelse(proteins$group == "Group 3","Metformin","Placebo")
proteins$date <- mdy(proteins$date)
# DROP GROUP - WILL GET THIS FROM THE RANDOMIZATION FILE
proteins <- subset(proteins,select=-c(group))

####  CYTOKINE DATA ######
all <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Preliminary CVD cytokine data/Jenny manifest v4 GRP2 9th Jan 2019.csv",na.strings = c("","NA"," ","error"))
# delete records with missing data
cytokine <- all[!is.na(all$IFN.g),]
cytokine$date <- cytokine$COLLECTIONDT
cytokine$analyticid <- cytokine$ANALYTICID
# get rid of unneeded variables
cytokine <- select(cytokine,-c("RN","PARENT_ID","SAMPLE_ID","STORAGETYPE","SAMPLEGROUP","BOX.ID","NO.","ROW","COL","pg.ml","X",
                               "COLLECTIONDT","ANALYTICID"))
cytokine$date <- mdy(cytokine$date)

# merge cytokine data and protein data and check if there are any records in one but not the other
cytokine_protein <- merge(cytokine, proteins, by=c("analyticid","date"),all.x = TRUE,all.y = TRUE)
test <- cytokine_protein[,c("analyticid","date","IFN.g","sp|P01024|CO3_HUMAN")]
write.csv(test,"H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/checking cytokine protein dates.csv")

####  RANDOMIZATION DATA ######
rand <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/metformin vs placebo.csv")

####  CLINICAL DATA ######
####  vars needed: height, weight, sex, age in months, HbA1c, visceral fat, LDL, HDL, waist, triglycerides, adipo, DBP, insulin dose #####
####  plus demographics: race/ethnicity, Tanner #####
anthro <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/031319_BOC031_LabValuesVisits_anthro.csv")
a1c <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/031319_BOC031_LabValuesVisits_a1c.csv")
labs <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/031319_BOC031_LabValuesVisits_Labs.csv")
adipo <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/030719_BOC031_Analytes_adipokines.csv")
jan <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/011518_BOC031_Data Pull.csv")
oct <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/BOC031Data Pull 10_11_18.csv")

# need link between visit number and dates in order to merge all the data
# here are the variables and the files they are found in
# height - anthro
# weight - anthro
# sex    - 011518 pull
# age    - 011518 pull
# A1c    - a1c
# % fat  - anthro
# LDL    - labs
# HDL    - labs
# trig   - labs
# waist  - 101118 pull
# adip   - adipokines
# DBP    - 101118 pull
# ins    - 101118 pull
# Tanner - 101118 pull


####  EFFLUX DATA ######