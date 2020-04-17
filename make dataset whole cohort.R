# make dataset for final analysis

library(lubridate)
library(dplyr)
library(tableone)
library(knitr)
library(data.table)
library(tidyr)
library(Hmisc)
library(childsds)
library(zoo)


####  RANDOMIZATION DATA ######
rand <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/metformin vs placebo.csv")
rand$analyticid <- rand$ANALYTICID
rand <- select(rand,-"ANALYTICID")

####  CLINICAL DATA ######
####  vars needed: height, weight, sex, age in months, HbA1c, visceral fat, LDL, HDL, waist, triglycerides, adipo, DBP, insulin dose #####
####  plus demographics: race/ethnicity, Tanner #####

# a1c data
a1c <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/031319_BOC031_LabValuesVisits_a1c.csv")
a1c <- a1c[a1c$Visit %in% c("26 week","Screening","26 Week"),]
a1c$visit[a1c$Visit %in% c("26 week","26 Week")] <- "6 month"
a1c$visit[a1c$Visit=="Screening"] <- "Baseline"
a1c$analyticid <- a1c$ID
a1c <- select(a1c,-c("ID","Visit","Method"))
alldata <- merge(a1c,rand,by=c("analyticid"),all.x=TRUE, all.y=TRUE)

# anthro file - only visit needed is 6 months
anthro <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/031319_BOC031_LabValuesVisits_anthro.csv")
anthro <- anthro[anthro$Visit %in% c("26 week","26 Week"),]
levels(anthro$Visit) <- c("26 week","26 Week","6 month")
anthro$visit[anthro$Visit %in% c("26 week","26 Week")] <- "6 month"
anthro$analyticid <- anthro$ID
#anthro$visit <- anthro$Visit
anthro <- select(anthro,-c("ID","Visit"))
alldata <- merge(alldata,anthro,by=c("analyticid","visit"),all.x=TRUE, all.y=FALSE)



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


# first clinical data pull
jan <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/011518_BOC031_Data Pull.csv")
jan$analyticid <- jan$ANALYTICID
jan <- select(jan,c("analyticid","Age","Gender","Race","Ethnicity"))
alldata <- merge(alldata,jan,by="analyticid",all.x=TRUE, all.y=FALSE)

# second clinical data pull
oct <- read.csv("H:/Endocrinology/Nadeau/T1D Exchange metformin and lipids/Data/Clinical data/BOC031Data Pull 10_11_18.csv")
oct <- oct[oct$Visit %in% c("A)baseline","I)26 week"),]
oct$analyticid <- oct$Analytic.ID
oct$visit[oct$Visit=="A)baseline"] <- "Baseline"
oct$visit[oct$Visit=="I)26 week"] <- "6 month"
oct$weight <- as.numeric(as.character(oct$Weight))
oct$height <- as.numeric(as.character(oct$Height))
oct$DEXAPercFat <- as.numeric(as.character(oct$DEXAPercFat))
oct <- select(oct,-c("Analytic.ID","Visit","Weight","Height","UnitsInsTotalPumpDload"))
alldata <- merge(alldata,oct, by=c("analyticid","visit"),all.x=TRUE,all.y=FALSE)
alldata$weight.x[is.na(alldata$weight.x)] <- alldata$weight.y[is.na(alldata$weight.x)] 
alldata$height.x[is.na(alldata$height.x)] <- alldata$height.y[is.na(alldata$height.x)]
alldata$weight <- alldata$weight.x
alldata$height <- alldata$height.x
alldata$DEXAPercFat.x[is.na(alldata$DEXAPercFat.x)] <- alldata$DEXAPercFat.y[is.na(alldata$DEXAPercFat.x)]
alldata$DEXAPercFat <- alldata$DEXAPercFat.x
alldata <- select(alldata,-c("weight.x","weight.y","height.x","height.y","DEXAPercFat.x","DEXAPercFat.y"))

# calculate BMI
alldata$bmi <- alldata$weight/((alldata$height/100)^2)
alldata$bmi_perc <- sds(alldata$bmi,
                       age = alldata$Age,
                       sex = alldata$Gender, male = "M", female = "F",
                       ref = cdc.ref,
                       item = "bmi",
                       type = "perc")*100

# calculate eIS
#exp (4.06154 ??? 0.01317 * waist [cm] ??? 1.09615 * insulin dose [daily units per kg] 
#     + 0.02027 * adiponectin [??g/mL] ??? 0.27168 * triglycerides [mmol/L (???0.00307 for mg/dL)] ??? 0.00733 * DBP [mm Hg])
alldata$units_per_kg <- as.numeric(as.character(alldata$UnitsInsTotal))/alldata$weight
alldata$adipo_ugml <- as.numeric(as.character(alldata$Value.ADIP))/1000
alldata$eis[!is.na(alldata$adipo_ugml)] <- exp(4.06154 - 0.01317 * as.numeric(as.character(alldata$WaistCircum[!is.na(alldata$adipo_ugml)])) 
                    - 1.09615 * alldata$units_per_kg[!is.na(alldata$adipo_ugml)]
                    + 0.02027 * alldata$adipo_ugml[!is.na(alldata$adipo_ugml)] 
                    -0.00307 * as.numeric(as.character(alldata$`Value.TG-NET`[!is.na(alldata$adipo_ugml)])) 
                    - 0.00733 * as.numeric(as.character(alldata$BldPrDia[!is.na(alldata$adipo_ugml)])))
# for visits missing adiponectin
alldata$eis[is.na(alldata$adipo_ugml)] <- exp(4.1075 - 0.01299 * as.numeric(as.character(alldata$WaistCircum[is.na(alldata$adipo_ugml)])) 
                                              - 1.05819 * alldata$units_per_kg[is.na(alldata$adipo_ugml)]
                                              -0.00354 * as.numeric(as.character(alldata$`Value.TG-NET`[is.na(alldata$adipo_ugml)])) 
                                              - 0.00802 * as.numeric(as.character(alldata$BldPrDia[is.na(alldata$adipo_ugml)])))

# fill in baseline characteristics across visits
vars <- alldata[,c("analyticid","visit","Age","Gender","Race","Ethnicity","TannerPubicH","TannerBreGen")]
vars <- vars[vars$visit=="Baseline",]
vars <- select(vars,-"visit")
alldata <- select(alldata,-c("Age","Gender","Race","Ethnicity","TannerPubicH","TannerBreGen"))
alldata <- merge(alldata,vars,by="analyticid",all.x = TRUE,all.y = TRUE)
alldata <- alldata[order(alldata$analyticid, alldata$visit),]

# fix non-numeric variables
alldata$bmiperc <- alldata$bmi_perc
alldata <- select(alldata,-bmi_perc)
nonnum <- c("analyticid","visit","date","Treatment.Group","UnitsInsTotalPumpOrLog","InsDeliveryMethod","Gender",
            "Race","Ethnicity","TannerPubicH","TannerBreGen")
temp <- colnames(alldata)
num <- temp[!(temp %in% nonnum)]
alldata[,num] <-  sapply(alldata[,num],function(x) as.numeric(as.character(x)))

# need to combine race/ethnicity 
alldata$raceeth[alldata$Ethnicity=="Hispanic or Latino"] <- "Hispanic"
alldata$raceeth[alldata$Ethnicity !="Hispanic or Latino" & alldata$Race=="White"] <- "Non-hispanic White"
alldata$raceeth[is.na(alldata$raceeth)] <- "Other"

# labels
var.labels <- colnames(alldata)
for(i in seq_along(alldata)){
  Hmisc::label(alldata[, i]) <- var.labels[i]
}
label(alldata$raceeth) = "Race/ethnicity"
label(alldata$TannerBreGen) = "Tanner by breasts/genitals"
label(alldata$TannerPubicH) = "Tanner by pubic hair"
label(alldata$Age)="Age"
label(alldata$Gender)="Gender"

# write final dataset
write.csv(alldata,file="H:\\Endocrinology\\Nadeau\\T1D Exchange metformin and lipids\\Data\\clinical data whole cohort.csv")


