#load libraries
library(readxl)
library(tidyr)
library(ggplot2)

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
#reshape data in long format
STACKED_d <- gather(d, "Date", "state", 3:11)

#how many states are there?
unique(STACKED_d$state)

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
}

STACKED_d$Date <- as.POSIXct(STACKED_d$Date)
ggplot(STACKED_d, aes(Date, state_num)) + geom_line() + facet_wrap(~MOATS)

