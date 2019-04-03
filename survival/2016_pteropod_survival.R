#load libraries
library(readxl)
library(tidyr)
library(ggplot2)
library(FSA)

d <- read_xlsx("~/Documents/GitHub/pteropod_pHxDO_metabolomics/survival/pteropod_pHxDO2016_masterdatasheet.xlsx", sheet = "sample IDs living", skip = 2)
weird_date <- as.numeric(colnames(d)[3:11], quote = FALSE)
good_dates <- as.Date(weird_date, origin = '1899-12-30')
d <- data.frame(d)
colnames(d)[3:11] <- as.character(good_dates)
colnames(d)[1:2] <- c("MOATS", "Jar")
#remove rows containing repetitive date info
d <- d[-grep("^42",d$`2016-11-22`),]
#remove rows containing NAs
d <- d[which(!is.na(d$`2016-11-22`)),]
#add column for first day
d$`2016-11-21` <- "L"
#reorder columns so dates are all next to each other
d <- d[,c(1,2,16,3:15)]

#remove MOATS 2 and 9
d <- d[-grep("M2|M9", d$MOATS),]

#reshape data in long format
STACKED_d <- gather(d, "Date", "state", 3:12)

#how many states are there?
states <- unique(STACKED_d$state)
length(states)

dsummary <- data.frame()
for(i in 1:nrow(d)){
  dtemp <- data.frame(table(t(d[i,3:12])))
  dtemp$ID <- paste(d$MOATS[i],d$Jar[i], sep = "_")
  dsummary <- rbind(dsummary, dtemp)
}

dsummary <- spread(dsummary, "Var1", "Freq")

#exclude animals that went missing and live but lost since we can't determine what happened with them
#and also 2L because we don't know where they came from

dsummary <- dsummary[which(is.na(dsummary$M) & is.na(dsummary$LL) & is.na(dsummary$`2L`)),]


#create MOATS and jar columns
dsummary$MOATS <- gsub("_.*","",dsummary$ID)
dsummary$MOATS <- gsub("M","",dsummary$MOATS)
dsummary$Jar <- gsub(".*_","",dsummary$ID)

#plot 
ggplot(dsummary, aes(factor(MOATS, levels = c(3,7,13,6,12,5,8,10,1,4,11)),L)) + geom_violin(aes(fill = MOATS)) + geom_boxplot(aes(fill = MOATS), width = 0.1) + theme_bw() + scale_fill_manual(values = c("mediumpurple1","turquoise3","mediumpurple1","forestgreen","coral1","coral1","mediumpurple1","turquoise3","forestgreen","coral1","turquoise3")) + xlab("MOATS ordered by treatment") + ylab("Survival (days)") + ggtitle("Treatment effect on Pteropod survival")

#read in treatment info
treatments <- read.csv("~/Documents/GitHub/Seawater-Chemistry-Analysis/2016-17_PteropodExp_WaterChem/PteropodWaterChem/Treatments.csv", stringsAsFactors = FALSE)
#merge treatment info with dsummary
dsummary <- merge(dsummary, treatments, by = "MOATS")

#plot by treatment
ggplot(dsummary, aes(Target_Treatment,L)) + geom_violin(aes(fill = Target_Treatment)) + geom_boxplot(aes(fill = Target_Treatment), width = 0.1) + theme_bw() + theme(axis.text.x = element_blank(), strip.text.x = element_text(size = 7))
ggplot(dsummary, aes(Treatment_abbv,L)) + geom_violin(aes(fill = Treatment_abbv)) + geom_boxplot(aes(fill = Treatment_abbv), width = 0.1) + theme_bw() + theme(axis.text.x = element_blank(), strip.text.x = element_text(size = 7)) + xlab("Treatment") + ylab("Survival (days)") + ggtitle("Treatment effect on Pteropod survival")

#plot density to check normality
ggplot(dsummary) + geom_density(aes(L, color = Treatment_abbv)) + theme(axis.text.x = element_blank(), strip.text.x = element_text(size = 7)) + xlab("Treatment") + ylab("Survival (days)") + ggtitle("Treatment effect on Pteropod survival")


#data not normal so kruskal test
#kruskal test
kruskal.test(L ~ as.factor(Treatment_abbv), data = dsummary)
#Kruskal-Wallis rank sum test
#data:  L by as.factor(Treatment_abbv)
#Kruskal-Wallis chi-squared = 11.983, df = 3, p-value = 0.007441

#test significant so dunn test 

#dunn test
#https://rcompanion.org/rcompanion/d_06.html
dunnTest(L ~ as.factor(Treatment_abbv), data = dsummary, method = "bh")

#Comparison         Z     P.unadj      P.adj
#1    HH - HL -1.719614 0.085502718 0.17100544
#2    HH - LH -1.131275 0.257939271 0.30952712
#3    HL - LH  0.736434 0.461466642 0.46146664
#4    HH - LL  1.594922 0.110729738 0.16609461
#5    HL - LL  3.093264 0.001979681 0.01187809
#6    LH - LL  2.714395 0.006639695 0.01991909

#

#alternative pairwise mann-whitney u test
#wilcox rank sum to determine which is different, need FDR
#HH:HL
pairwise.wilcox.test(dsummary[grep("HH|HL",dsummary$Treatment_abbv),"L"], as.factor(dsummary[grep("HH|HL",dsummary$Treatment_abbv),"Treatment_abbv"]), p.adjust.method = "BH")
#Pairwise comparisons using Wilcoxon rank sum test 
#data:  dsummary[grep("HH|HL", dsummary$Treatment_abbv), "L"] and as.factor(dsummary[grep("HH|HL", dsummary$Treatment_abbv), "Treatment_abbv"]) 
#HH:HL 0.068
#P value adjustment method: BH

wilcox.test(L ~ as.factor(Treatment_abbv), data = dsummary[grep("HH|HL",dsummary$Treatment_abbv),])
#Wilcoxon rank sum test with continuity correction
#data:  L by as.factor(Treatment_abbv)
#W = 1661.5, p-value = 0.06828
#alternative hypothesis: true location shift is not equal to 0

#HH:LH
wilcox.test(L ~ as.factor(Treatment_abbv), data = dsummary[grep("HH|LH",dsummary$Treatment_abbv),])
#Wilcoxon rank sum test with continuity correction
#data:  L by as.factor(Treatment_abbv)
#W = 2999.5, p-value = 0.2367
#alternative hypothesis: true location shift is not equal to 0

#HH:LL
wilcox.test(L ~ as.factor(Treatment_abbv), data = dsummary[grep("HH|LL",dsummary$Treatment_abbv),])
#Wilcoxon rank sum test with continuity correction
#data:  L by as.factor(Treatment_abbv)
#W = 3742.5, p-value = 0.1022
#alternative hypothesis: true location shift is not equal to 0


#HL:LH
wilcox.test(L ~ as.factor(Treatment_abbv), data = dsummary[grep("HL|LH",dsummary$Treatment_abbv),])
#Wilcoxon rank sum test with continuity correction
#data:  L by as.factor(Treatment_abbv)
#W = 2100, p-value = 0.5174
#alternative hypothesis: true location shift is not equal to 0

#HL:LL
wilcox.test(L ~ as.factor(Treatment_abbv), data = dsummary[grep("HL|LL",dsummary$Treatment_abbv),])
#Wilcoxon rank sum test with continuity correction
#data:  L by as.factor(Treatment_abbv)
#W = 2549, p-value = 0.00215
#alternative hypothesis: true location shift is not equal to 0




#LH:LL
wilcox.test(L ~ as.factor(Treatment_abbv), data = dsummary[grep("LH|LL",dsummary$Treatment_abbv),])
#Wilcoxon rank sum test with continuity correction
#data:  L by as.factor(Treatment_abbv)
#W = 3939, p-value = 0.01038
#alternative hypothesis: true location shift is not equal to 0



for(i in 1:length(STACKED_d$state)){
  if(STACKED_d$state[i] == "D"){
    STACKED_d$state_num[i] <- 0
  }
  if(STACKED_d$state[i] == "M"){
    STACKED_d$state_num[i] <- 0.5
  }
  if(STACKED_d$state[i] == "L"){
    STACKED_d$state_num[i] <- 1
  }
  if(STACKED_d$state[i] == "2L"){
    STACKED_d$state_num[i] <- 2
  }
  if(STACKED_d$state[i]== "LL"){
    STACKED_d$state_num[i] <- 0.75
  }
  else{
    STACKED_d$state_num[i] <- -1
  }
}

STACKED_d$Date <- as.POSIXct(STACKED_d$Date)
ggplot(STACKED_d, aes(Date, state_num)) + geom_line() + facet_wrap(~MOATS)

