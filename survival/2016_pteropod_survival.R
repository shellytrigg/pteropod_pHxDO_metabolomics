#load libraries
library(readxl)
library(tidyr)

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

