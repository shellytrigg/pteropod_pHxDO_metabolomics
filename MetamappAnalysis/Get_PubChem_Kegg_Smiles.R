###Script for getting pubchem, inchi key, and KEGG names for all compounds

#STEP 1: Retrieve known compounds from processed data files from WCMC
#Copy and paste names, inchi key, pubchem, and KEGG names for all known compounds in WCMC processed data files 
#saved this as All_known_compounds.csv in the data folder in the repo
#Changed Gal-Gal-Cer(d18:1/16:0) or Lactosylceramide(d18:1/16:0) [M+Na]+ to Gal-Gal-Cer(d18:1/16:0) and took only the first inchikey listed


#STEP 2: Retrieve up-to-date metabolite and lipid data to try to find more IDs
#download .xml file for ALl metabolites from HMDB (http://www.hmdb.ca/downloads) and .sdf file from LipidMaps (https://www.lipidmaps.org/data/structure/download.php)
#Downloaded metabolite data from HMDB (Jan 2019 release) and lipid data from LipidMaps (June 2019 release)
#saved in github repo under data

####get inchi key and pubchem compound id for all hmdb metabolites

#run this code in bash
awk '{if($0~/<creation_date>/)print "@@@"$0;else print $0}' /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/hmdb_metabolites.xml | \
grep -E '<.{0,1}creation_date.{0,1}>|<.{0,1}inchikey.{0,1}>|<.{0,1}pubchem_compound_id.{0,1}>' | tr -d ' ' | tr '\n' '\t' | tr '@@@' '\n' | \
awk '{if($1~/creation_date/)print $0}' |\
sed 's/<//g' | sed 's/>//g' | sed 's/\///g' | sed 's/pubchem_compound_id//g'| sed 's/inchikey//g' |sed 's/creation_date//g' > /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/hmdb_entry_inchikey_pubchemID.tsv

###QC for HMDB table
awk -F"\t" '{print length($1), length($2), length($3)}' /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/hmdb_entry_inchikey_pubchemID.tsv | sort | uniq -c 
94048 21 27 0
45 21 27 1
29 21 27 2
319 21 27 3
1573 21 27 4
2486 21 27 5
2929 21 27 6
3409 21 27 7
8761 21 27 8
501 21 27 9

###investigating why 45 lines have 1 character...
awk -F"\t" '{if(length($3)==1)print $0}' /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/hmdb_entry_inchikey_pubchemID.tsv | less
### These are all zeros. So now searching if <pubchem_compound_id>0<pubchem_compound_id/> exists
#it does
grep "<pubchem_compound_id>0</pubchem_compound_id>" /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/hmdb_metabolites.xml | less

<pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  <pubchem_compound_id>0</pubchem_compound_id>
  
####get inchi key and pubchem compound id for all lipid maps metabolites
  
awk '{if($0~/> <LM_ID>/)print "@@@> <LM_ID>";else print $0}' /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSD_2019-06-06.sdf |\
grep -E -A1 'LM_ID|PUBCHEM_CID|INCHI_KEY' | \
tr '\n' '\t' | tr '@@@' '\n' | tr -d '\t' | sed 's/--//g' | tr -d '>' | sed 's/<PUBCHEM_CID//g' | sed 's/<INCHI_KEY//g' \
| awk '{if($1~/</)print $2"\t"$3}' > /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSDF_pubchem_inchi.tsv

#QC

2 1 27
8 2 27
313 27 0
42 3 27
230 4 27
888 5 27
1653 6 27
8342 7 27
27400 8 27
3978 9 27

#What has a zero inchikey?
grep -B50 -A100 PTMZTTPJFDLIOR-HCCCIJMNSA-N /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSD_2019-06-06.sdf | less
##This shows that these entries don't have pubchem IDs in the LMDB
#So I need to enter 0 for these pubchem IDs
awk '{if($0~/> <LM_ID>/)print "@@@> <LM_ID>";else print $0}' /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSD_2019-06-06.sdf |\
grep -E -A1 'LM_ID|PUBCHEM_CID|INCHI_KEY' | \
tr '\n' '\t' | tr '@@@' '\n' | tr -d '\t' | sed 's/--//g' | tr -d '>' | sed 's/<PUBCHEM_CID//g' | sed 's/<INCHI_KEY//g' \
| awk '{if($1~/</)print $2"\t"$3}' |\
awk -F"\t" '{if(length($2)==0)print 0FS$1;else print $0}' > /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSDF_pubchem_inchi.tsv

#QC again
awk -F"\t" '{print length($1), length($2)}' /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSDF_pubchem_inchi.tsv | sort | uniq -c 

315 1 27
8 2 27
42 3 27
230 4 27
888 5 27
1653 6 27
8342 7 27
27400 8 27
3978 9 27

#GOOD!
#Print all inchikeys with pubchem ID == 0 and look up pubchem ID in Pubchem ID exchange
awk -F"\t" '{if($1==0)print $2}' /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSDF_pubchem_inchi.tsv

#Go to identifier exchange and convert inchi keys to CIDs and download
#move downloaded file to data folder in repo
#replace entries with pubchem ID of 0 to the new retrieved IDS

#QC to make sure it works
zcat < /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/3561635504380663618.txt.gz | awk -F"\t" '{print $2FS$1}' | cat /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSDF_pubchem_inchi.tsv - | awk -F"\t" '{if($1!=0)print $0}' | awk -F"\t" '{print length($1), length($2)}'| sort | uniq -c
2 1 27
8 2 27
42 3 27
230 4 27
889 5 27
1654 6 27
8345 7 27
27445 8 27
4241 9 27

#it does so save new file
zcat < /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/3561635504380663618.txt.gz | awk -F"\t" '{print $2FS$1}' | cat /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSDF_pubchem_inchi.tsv - | awk -F"\t" '{if($1!=0)print $0}' > LMSDF_pubchem_inchi_nozeros.tsv

#rename file

mv LMSDF_pubchem_inchi_nozeros.tsv /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSDF_pubchem_inchi.tsv



#STEP 3: read in data
All_names <- read.csv("~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/All_known_compounds.csv", stringsAsFactors = FALSE)
HMDB <- read.csv("~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/hmdb_entry_inchikey_pubchemID.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
LMDB <- read.csv("~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSDF_pubchem_inchi.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
crab <- read.csv("~/Desktop/UCSD-SALK-UW/NOAA/metabolomics/all_cmpds_analyzed/2groups_compared_FC/name_files/All_known_444cmpds_pub_inch_kegg_binbase.csv", header = TRUE, stringsAsFactors = TRUE)

#STEP 4: format data
#add headers
colnames(HMDB) <- c("CreationDate","InChI.Key","PubChem", "NA")
colnames(LMDB) <- c("PubChem", "InChI.Key")

#remove creation date and extra column in HMDB
HMDB <- HMDB[,2:3]

#STEP 5: merge HMDB and LMDB files with All_names

All_names <- merge(All_names, HMDB,by = "InChI.Key", all.x = TRUE)
All_names <- merge(All_names, LMDB, by = "InChI.Key", all.x = TRUE)
All_names <- merge(All_names, crab[,c(4:6)], by = "InChI.Key", all.x = TRUE)

#rename columns

colnames(All_names) <- c("InChI.Key", "BinBase.name", "PubChem", "KEGG", "HMDBpc", "LMDBpc", "CrabPC", "CrabKegg")


#remove duplicate lines with same pubchem IDs
All_names_STACKED <-tidyr::gather(All_names, "source", "ID", c(3,5,6,7))
All_names_STACKED <- All_names_STACKED[,-5]
All_names_STACKED <- unique(All_names_STACKED)

#rename column back to pubchem
colnames(All_names_STACKED)[5] <- "PubChem"
#convert back to wide format excluding KEGG ids and pubchem IDs = NA
All_names_uniquePubChem <- All_names_STACKED[which(!(is.na(All_names_STACKED$PubChem))),]



##compounds are sometimes detected in both lipidomics and metabolomics and don't have the same binbase name, but do have the same pubchem ID
#these should probably be considered only once...


#consolidate KEGG column
All_names_uniquePubChem$CrabKegg <- gsub("not found","", All_names_uniquePubChem$CrabKegg)

All_names_uniquePubChem_STACKED <- tidyr::gather(All_names_uniquePubChem, "KEGG","ID", 3:4)

All_names_uniquePubChem_STACKED <- All_names_uniquePubChem_STACKED[,-4]
All_names_uniquePubChem_STACKED <- unique(All_names_uniquePubChem_STACKED)
colnames(All_names_uniquePubChem_STACKED)[4] <- "KEGG"

### remove binbase names since this isn't important right now and is causing duplicates
All_names_uniquePubChem_STACKED <- All_names_uniquePubChem_STACKED[,-2]
All_names_uniquePubChem_STACKED <- unique(All_names_uniquePubChem_STACKED)


All_names_uniquePubChem_KEGG <- All_names_uniquePubChem_STACKED[grep("C", All_names_uniquePubChem_STACKED$KEGG),]

All_names_uniquePubChem_noCKEGG<- All_names_uniquePubChem_STACKED[-grep("C", All_names_uniquePubChem_STACKED$KEGG),]

All_names_uniquePubChem_noCKEGGnotin <- All_names_uniquePubChem_noCKEGG[which(!(All_names_uniquePubChem_noCKEGG$InChI.Key %in% All_names_uniquePubChem_KEGG$InChI.Key)),]
#add NAs to KEGG column where they were blank
All_names_uniquePubChem_noCKEGGnotin[All_names_uniquePubChem_noCKEGGnotin==''] <- NA
#remove duplicate lines 
All_names_uniquePubChem_noCKEGGnotin <- unique(All_names_uniquePubChem_noCKEGGnotin)

All_names_uniquePubChem_KEGG <- rbind(All_names_uniquePubChem_KEGG, All_names_uniquePubChem_noCKEGGnotin)

View(All_names_uniquePubChem_KEGG[which(duplicated(All_names_uniquePubChem_KEGG$InChI.Key)),])

#only keep inchikeys with the lowest pubchem ID
df.agg <- aggregate(PubChem ~ InChI.Key, All_names_uniquePubChem_KEGG, min)

All_names_uniquePubChem_KEGG <- merge(All_names_uniquePubChem_KEGG, df.agg)

#now we have only unique inchikeys
duplicated(All_names_uniquePubChem_KEGG[,1])

#remove lines without KEGG IDs if they are duplicated
dup_cpds <- All_names_uniquePubChem[duplicated(All_names_uniquePubChem[,1]),1]


#get smiles :
write.csv(All_names_uniquePubChem_KEGG,"~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/unique_inchi_pubchem_kegg.csv", quote = FALSE, row.names = FALSE)
#copy pubchem IDs and go to exchange to get smiles
#download file and save
#zcat < /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/1087426400801374154.txt.gz | tr '\t' ',' | awk -F, 'NR==FNR{a[$1]=$2;next}{print $0","a[$2]}' - /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/unique_inchi_pubchem_kegg.csv > unique_inchi_pubchem_kegg_smiles.csv



###Read in univariate data

library(readxl)
metab_uni<- read_xlsx("~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/pteropod_genmetabANDlipids_univariate_RESULTSsupplement.xlsx", sheet =1)
lipid_uni<- read_xlsx("~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/pteropod_genmetabANDlipids_univariate_RESULTSsupplement.xlsx", sheet =2)

colnames(metab_uni)[1] <- "BinBasename_R"
colnames(lipid_uni)[1] <- "Identifier_R"


#read in ID keys
metab_ID <- read.table("../MetamappAnalysis/metab_ID_Key.txt", sep ="\t", stringsAsFactors = FALSE,header = TRUE, quote ="")
  
lipid_ID <- read.table("../MetamappAnalysis/lipid_ID_Key.txt", sep ="\t", stringsAsFactors = FALSE,header = TRUE, quote ="",comment.char = "")




#get names with decimals instead of spaces
#names_df <- All_names[,1:2]
#names_df$alt <- trimws(names_df$BinBase.name)
#names_df$alt <- gsub(" ",".",names_df$alt)
#names_df$alt <- gsub("\\(",".",names_df$alt)
#names_df$alt <- gsub("\\)",".",names_df$alt)
#names_df$alt <- gsub("-",".",names_df$alt)
#names_df$alt <- gsub(":",".",names_df$alt)
#names_df$alt <- gsub(",",".",names_df$alt)
#names_df$alt <- gsub("\\/",".",names_df$alt)
#names_df$alt <- gsub(";",".",names_df$alt)

#names_df <- unique(names_df)
#colnames(names_df)[3] <- "analyte"

#add cholesterol inchikey to cholesterol.nist compound
#names_df$InChI.Key[2] <- "HVYWMOMLDIMFJA-DPAQBDIFSA-N"

#add methionine inchikey to methionine minor compound
#names_df$InChI.Key[1] <- "FFEARJCKVFRZRR-BYPYZUCNSA-N"


metab_uni_names <- merge(metab_uni[,c(1,4,8,14,20)],metab_ID[,-c(3:4)], by = "BinBasename_R")
#reorder columns
metab_uni_names <- metab_uni_names[,c(1:5,7,6,8)]
#rename columns to match lipid data colnames
colnames(metab_uni_names)[1] <- "Identifier_R"
colnames(metab_uni_names)[6] <- "Identifier"
colnames(metab_uni_names)[7] <- "Annotation"


lipid_uni_names <- merge(lipid_uni[,c(1,4,8,14,20)],lipid_ID, by = "Identifier_R")


all_analytes <- rbind(metab_uni_names,lipid_uni_names)
###Combine data with other data
pubch_kegg_smiles <- read.csv("~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/unique_inchi_pubchem_kegg_smiles.csv", stringsAsFactors = FALSE)
colnames(pubch_kegg_smiles)[4] <- "SMILES"

#pubchem table still contains 1 duplicate where one line has a kegg and the other doesn't. It is pubchem ID 24699
#remove the line missing the kegg ID
pubch_kegg_smiles <- pubch_kegg_smiles[-which(pubch_kegg_smiles$PubChem=="24699" & is.na(pubch_kegg_smiles$KEGG)),]


#merge pubchem, kegg, and smiles IDs with all analytes by their InChI.Key
###NEED TO VALIDATE THIS
df <- merge(all_analytes,pubch_kegg_smiles, by = "InChI.Key")


#read in RSD, FC, and Cohen D info
RSD_FC_cohen <- read.table("../FoldChange_CohenD/all_cmpds_FC_RSD_cohen.tsv", sep = "\t",header = TRUE,stringsAsFactors = FALSE, quote = "", comment.char = "")


#merge RSD data with anova data
RSD_FC_cohen_pub <- merge(RSD_FC_cohen,df, by = "Identifier")


#df <- merge(df, RSD_FC_cohen, by = "InChI.Key")

#df <- unique(df)

#df <- df[,c(8,9,10,7,3)]

#colnames(df) <- c("PubChem_ID", "KEGG_ID", "SMILES", "Compound_Name", "Chi_overall_pvalue")
#write.csv(df, "~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/metamapp_input.csv", quote = FALSE, row.names = FALSE)


###METAMAPP
#copied and pasted arbitrary FC values and uploaded to metamapp

##need to resolve duplicated compound names (see email from Dinesh)

#print duplicated compound names
#dups <- df[which(duplicated(df$Compound_Name)),"Compound_Name"]
dups <- df[which(duplicated(df$PubChem)),"PubChem"]

View(df[which(df$PubChem %in% dups),])

#removed duplicates manually and tried to keep the ones that had KEGG IDs

#went through and manually removed duplicated pubchem IDs (most of these were from compounds detected in both GC/MS and LC/MS with different compound names but same pubchem IDs)
#saved this file as metamapp_input_nodups.csv

#copied and pasted header from Metamapp example data into metamapp_input_nodups.csv
#Uploaded to MetaMapp and downloaded sif and attributes file

node_attr <- read.table("data/node_attributes_chemsim_krp_07-1.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)


