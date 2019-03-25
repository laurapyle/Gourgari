# make dataset for final analysis

library(lubridate)
library(dplyr)
library(tableone)
library(knitr)
library(data.table)
library(tidyr)
library(Hmisc)
library(childsds)

####  EFFLUX DATA ######
efflux1 <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Efflux data/Jenny proteome efflux for me.csv")
efflux1 <- efflux1[!is.na(efflux1$efflux.value),]
efflux1$analyticid <- efflux1$ANALYTICID
efflux1$date <- efflux1$COLLECTIONDT
efflux1 <- select(efflux1,c("analyticid","date","efflux.value"))
efflux2 <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Efflux data/Copy of Jenny manifest v4 GRP2 march 2019.csv")
efflux2 <- efflux2[!is.na(efflux2$cholesterol.efflux),]
efflux2$analyticid <- efflux2$ANALYTICID
efflux2$date <- efflux2$COLLECTIONDT
efflux2$efflux.value <- efflux2$cholesterol.efflux
efflux2 <- select(efflux2,c("analyticid","date","efflux.value"))
alldata <- rbind(efflux1,efflux2)
alldata <- alldata[order(alldata$analyticid, alldata$date),]
# emailed Jenny - there are two different sets of results for the same patient in these two files
# after she responds, decide which records to keep
# then go through the rest of this code and make sure it merges everything correctly

# figure out who has duplicate efflux results
check <- alldata %>% add_count(analyticid) %>% filter(n==4)
write.csv(check,"H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Efflux data/duplicate results.csv")


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
# add visit number to protein dataset - first visit is baseline, second is 6 mo
proteins <- proteins[order(proteins$analyticid, proteins$date),]
# only keep paired observations, otherwise we can't tell what visit is which
paired <- proteins$analyticid[duplicated(proteins$analyticid)]
proteins <- proteins[which(proteins$analyticid %in% paired),]
visits <- rep(c("Baseline","6 month"),length(paired))
proteins$visit <- visits
alldata <- merge(alldata,proteins,by=c("analyticid","date"),all.x=TRUE, all.y=FALSE)
alldata <- alldata[order(alldata$analyticid, alldata$date),]

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
alldata <- merge(alldata,cytokine,by=c("analyticid","visit"),all.x=TRUE, all.y=FALSE)
alldata <- alldata[order(alldata$analyticid, alldata$date),]

####  RANDOMIZATION DATA ######
rand <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/metformin vs placebo.csv")
rand$analyticid <- rand$ANALYTICID
rand <- select(rand,-"ANALYTICID")
alldata <- merge(alldata,rand,by="analyticid",all.x = TRUE,all.y=TRUE)
# get rid of missing visit numbers - these were in cytokine file but not proteins and we won't know what visit they are
alldata <- alldata[!is.na(alldata$visit),]

####  CLINICAL DATA ######
####  vars needed: height, weight, sex, age in months, HbA1c, visceral fat, LDL, HDL, waist, triglycerides, adipo, DBP, insulin dose #####
####  plus demographics: race/ethnicity, Tanner #####

# anthro file - only visit needed is 6 months
anthro <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/031319_BOC031_LabValuesVisits_anthro.csv")
anthro <- anthro[anthro$Visit %in% c("26 week","26 Week"),]
anthro$Visit <- rep("6 month",nrow(anthro))
anthro$analyticid <- anthro$ID
anthro$visit <- anthro$Visit
anthro <- select(anthro,-c("ID","Visit"))
alldata <- merge(alldata,anthro,by=c("analyticid","visit"),all.x=TRUE, all.y=FALSE)
alldata <- alldata[order(alldata$analyticid, alldata$date),]

# a1c data
a1c <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/031319_BOC031_LabValuesVisits_a1c.csv")
a1c <- a1c[a1c$Visit %in% c("26 week","Screening","26 Week"),]
a1c$visit[a1c$Visit %in% c("26 week","26 Week")] <- "6 month"
a1c$visit[a1c$Visit=="Screening"] <- "Baseline"
a1c$analyticid <- a1c$ID
a1c <- select(a1c,-c("ID","Visit","Method"))
alldata <- merge(alldata,a1c,by=c("analyticid","visit"),all.x=TRUE, all.y=FALSE)
alldata <- alldata[order(alldata$analyticid, alldata$date),]

# labs
labs <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/031319_BOC031_LabValuesVisits_Labs.csv")
# adipo data
adipo <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/030719_BOC031_Analytes_adipokines.csv")
colnames(adipo) <- c("ID","Visit","Analyte","Value","Units")
# deduplicate the adipokines
adipo <- unique(adipo)
adipo <- adipo[order(adipo$ID,adipo$Visit,adipo$Analyte),]
labs <- rbind(labs,adipo)
# make lab dataset one record per visit
labs <- labs[labs$Visit %in% c("Randomization","26 week","26 Week"),]
labs$visit[labs$Visit=="Randomization"] <- "Baseline"
labs$visit[labs$Visit %in% c("26 week","26 Week")] <- "6 month"
labs$analyticid <- labs$ID
labs <- select(labs,-c("ID","Visit","Units"))
labs <- labs[order(labs$analyticid, labs$visit),]
# subset data by patient and reshape
by_patient_sort<-function(ID,visit,data){
  #for each unique ID in the dataset,
  temp<-lapply(unique(ID), function(x){
    #create a dat.temp that subsets the data by id. 
    dat.temp <- subset(data, ID == x )
    ##test on single subject
    #dat.temp<-subset(nont1d,nont1d$Random_ID==146184)
    ###new code
    dat.temp <- reshape(dat.temp,idvar = c("analyticid","visit"),timevar = "Analyte",v.names = "Value",direction="wide")
  })
  #this binds together all of the mini patient datasets
  dat<-do.call(rbind,temp)
}
#use this to call the function
labs_wide <- by_patient_sort(labs$analyticid,labs$visit,labs) 
alldata <- merge(alldata,labs_wide,by=c("analyticid","visit"),all.x=TRUE, all.y=FALSE)
alldata <- alldata[order(alldata$analyticid, alldata$date),]

# first clinical data pull
jan <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/011518_BOC031_Data Pull.csv")
jan$analyticid <- jan$ANALYTICID
jan <- select(jan,c("analyticid","Age","Gender","Race","Ethnicity"))
alldata <- merge(alldata,jan,by="analyticid",all.x=TRUE, all.y=FALSE)
alldata <- alldata[order(alldata$analyticid, alldata$date),]

# second clinical data pull
oct <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/BOC031Data Pull 10_11_18.csv")
oct <- oct[oct$Visit=="A)baseline",]
oct$analyticid <- oct$Analytic.ID
oct$visit <- rep("Baseline",nrow(oct))
oct$weight <- as.numeric(as.character(oct$Weight))
oct$height <- as.numeric(as.character(oct$Height))
oct$DEXAPercFat <- as.numeric(as.character(oct$DEXAPercFat))
oct <- select(oct,-c("Analytic.ID","Visit","Weight","Height","UnitsInsTotalPumpDload"))
alldata <- merge(alldata,oct_keep, by=c("analyticid","visit"),all.x=TRUE,all.y=FALSE)
alldata$weight.x[is.na(alldata$weight.x)] <- alldata$weight.y[is.na(alldata$weight.x)] 
alldata$height.x[is.na(alldata$height.x)] <- alldata$height.y[is.na(alldata$height.x)]
alldata$weight <- alldata$weight.x
alldata$height <- alldata$height.x
alldata <- select(alldata,-c("weight.x","weight.y","height.x","height.y"))

# calculate BMI
alldata$bmi <- alldata$weight/((alldata$height/100)^2)
alldata$bmi_perc <- sds(alldata$bmi,
                       age = alldata$Age,
                       sex = alldata$Gender, male = "M", female = "F",
                       ref = cdc.ref,
                       item = "bmi",
                       type = "perc")*100

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

