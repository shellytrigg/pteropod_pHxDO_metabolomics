###Script for getting pubchem, inchi key, and KEGG names for all compounds

########################################################
####################### STEP 1 #########################
########################################################

#Retrieve known compounds from processed data files from WCMC
#Copy and paste names, inchi key, pubchem, and KEGG names for all known compounds in WCMC processed data files 
#saved this as All_known_compounds.csv in the data folder in the repo
#Changed Gal-Gal-Cer(d18:1/16:0) or Lactosylceramide(d18:1/16:0) [M+Na]+ to Gal-Gal-Cer(d18:1/16:0) and took only the first inchikey listed

########################################################
####################### STEP 2 #########################
########################################################
# Retrieve up-to-date metabolite and lipid data to try to find more IDs
#download .xml file for ALl metabolites from HMDB (http://www.hmdb.ca/downloads) and .sdf file from LipidMaps (https://www.lipidmaps.org/data/structure/download.php)
#Downloaded metabolite data from HMDB (Jan 2019 release) and lipid data from LipidMaps (June 2019 release)
#saved in github repo under data

####get inchi key and pubchem compound id for all hmdb metabolites

#run this code in bash

```{bash, eval = FALSE}
awk '{if($0~/<creation_date>/)print "@@@"$0;else print $0}' /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/hmdb_metabolites.xml | \
grep -E '<.{0,1}creation_date.{0,1}>|<.{0,1}inchikey.{0,1}>|<.{0,1}pubchem_compound_id.{0,1}>' | tr -d ' ' | tr '\n' '\t' | tr '@@@' '\n' | \
awk '{if($1~/creation_date/)print $0}' |\
sed 's/<//g' | sed 's/>//g' | sed 's/\///g' | sed 's/pubchem_compound_id//g'| sed 's/inchikey//g' |sed 's/creation_date//g' > /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/hmdb_entry_inchikey_pubchemID.tsv
```
###QC for HMDB table
```{bash, eval = FALSE}
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
```


 investigating why 45 lines have 1 character...
```{bash, eval = FALSE}
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
```
# get inchi key and pubchem compound id for all lipid maps metabolites
```{bash, eval = FALSE} 
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
```

#What has a zero inchikey?
```{bash, eval = FALSE}
grep -B50 -A100 PTMZTTPJFDLIOR-HCCCIJMNSA-N /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSD_2019-06-06.sdf | less
##This shows that these entries don't have pubchem IDs in the LMDB
#So I need to enter 0 for these pubchem IDs
awk '{if($0~/> <LM_ID>/)print "@@@> <LM_ID>";else print $0}' /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSD_2019-06-06.sdf |\
grep -E -A1 'LM_ID|PUBCHEM_CID|INCHI_KEY' | \
tr '\n' '\t' | tr '@@@' '\n' | tr -d '\t' | sed 's/--//g' | tr -d '>' | sed 's/<PUBCHEM_CID//g' | sed 's/<INCHI_KEY//g' \
| awk '{if($1~/</)print $2"\t"$3}' |\
awk -F"\t" '{if(length($2)==0)print 0FS$1;else print $0}' > /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSDF_pubchem_inchi.tsv
```

#QC again
```{bash, eval = FALSE}
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
```

#GOOD! Print all inchikeys with pubchem ID == 0 and look up pubchem ID in Pubchem ID exchange
```{bash, eval = FALSE}
awk -F"\t" '{if($1==0)print $2}' /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSDF_pubchem_inchi.tsv
```
#Go to identifier exchange and convert inchi keys to CIDs and download
#move downloaded file to data folder in repo
#replace entries with pubchem ID of 0 to the new retrieved IDS

#QC to make sure it works
```{bash, eval = FALSE}
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
```
#it does so save new file
```{bash, eval = FALSE}
zcat < /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/3561635504380663618.txt.gz | awk -F"\t" '{print $2FS$1}' | cat /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSDF_pubchem_inchi.tsv - | awk -F"\t" '{if($1!=0)print $0}' > LMSDF_pubchem_inchi_nozeros.tsv
```
#rename file
```{bash, eval = FALSE}
mv LMSDF_pubchem_inchi_nozeros.tsv /Users/Shelly/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSDF_pubchem_inchi.tsv
```

########################################################
####################### STEP 3 #########################
########################################################

# Create table with all chemical identifiers needed for MetaMapp for all known compounds
# This table will later be merged wtih metabolite stats data


# Read in data
All_names <- read.csv("~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/All_known_compounds.csv", stringsAsFactors = FALSE)
#HMDB <- read.csv("~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/hmdb_entry_inchikey_pubchemID.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#LMDB <- read.csv("~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/LMSDF_pubchem_inchi.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
HMDB <- read.table("data/hmdb_metabolites_classes.txt", header = F, sep = "\t", stringsAsFactors = F,quote = "")

LMDB <- read.table("data/LMSDF_lipid_classes.txt", header = F, sep = "\t", stringsAsFactors = F)
crab <- read.csv("~/Documents/UCSD-SALK-UW/NOAA/metabolomics/all_cmpds_analyzed/2groups_compared_FC/name_files/All_known_444cmpds_pub_inch_kegg_binbase.csv", header = TRUE, stringsAsFactors = FALSE)

CTS <- read.csv("data/cts-20201107221856.csv", header = T, stringsAsFactors = F)

#STEP 4: format data
#add headers
#colnames(HMDB) <- c("CreationDate","InChI.Key","PubChem", "NA")
#colnames(LMDB) <- c("PubChem", "InChI.Key")

colnames(HMDB) <- c("Creation.Date", "InChI.Key", "Super.Class", "Main.Class", "Sub.Class", "Molecular.Framework", "KEGG.ID","PubChem.ID")

colnames(LMDB) <- c("LM.ID", "Super.Class", "Main.Class", "Sub.Class", "InChI.Key", "PubChem", "KEGG")

colnames(crab)[4] <- "InChI.Key" 


colnames(CTS) <- c("InChI.Key", "KEGG_CTS", "PubChem.CID")



#replace the space preceding  FIJFPUAJUDAZEY-MNDXXDKYSA-N with nothing
All_names$InChI.Key <- gsub(" ","", All_names$InChI.Key) 
#remove ? at the end of inchi-key
All_names$InChI.Key <- gsub("\\?","", All_names$InChI.Key) 
#fix typo of -M to -N
All_names$InChI.Key <- gsub("-M","-N", All_names$InChI.Key) 
#fix typo of - to -N
All_names$InChI.Key <- gsub("-$","-N", All_names$InChI.Key) 



#try to find PubChem IDs for compounds that don't have them listed by searching 
# the identifier exchange service (https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi)
# using their InkChiKeys

All_names[which(All_names$InChI.Key=="AHWDQDMGFXRVFB-UHFFFAOYSA-N"),"PubChem"] <- "13225"
All_names[which(All_names$InChI.Key=="BTVWZWFKMIUSGS-UHFFFAOYSA-N"),"PubChem"] <- "68410"
All_names[which(All_names$InChI.Key=="DLIKSSGEMUFQOK-SFTVRKLSSA-N"),"PubChem"] <- "92794"
All_names[which(All_names$InChI.Key=="DSCFFEYYQKSRSV-FEPQRWDDSA-N"),"PubChem"] <- "164619"
All_names[which(All_names$InChI.Key=="LKDRXBCSQODPBY-OEXCPVAWSA-N"),"PubChem"] <- "439312"
All_names[which(All_names$InChI.Key=="OOWQBDFWEXAXPB-UHFFFAOYSA-N"),"PubChem"] <- "72733"
All_names[which(All_names$InChI.Key=="SRBFZHDQGSBBOR-SOOFDHNKSA-N"),"PubChem"] <- "10975657"
All_names[which(All_names$InChI.Key=="WSHCLKDGQAEXNQ-MRGADFHSSA-N"),"PubChem"] <- "52923862"






#remove creation date and extra column in HMDB
#HMDB <- HMDB[,2:3]

# merge HMDB and LMDB files with All_names
All_names_hmdb <- merge(HMDB,All_names,by = "InChI.Key", all.y = TRUE)
#there are three compounds where the pubchem IDs don't match; so need to correct these
# I manually looked them up on pubchem
View(All_names_hmdb[which(All_names_hmdb$PubChem != All_names_hmdb$PubChem.ID),])
#keep 1110 for succinic acid
All_names_hmdb[which(All_names_hmdb$BinBase.name == "succinic acid"),"PubChem.ID"] <- 1110
#keep 1118 for sulfuric acid
All_names_hmdb[which(All_names_hmdb$BinBase.name == "sulfuric acid"),"PubChem.ID"] <- 1118
#keep 6262 for ornithine
All_names_hmdb[which(All_names_hmdb$BinBase.name == "ornithine"),"PubChem"] <- 6262
#check to make sure conversions worked
View(All_names_hmdb[which(All_names_hmdb$PubChem != All_names_hmdb$PubChem.ID),])
#it did!
#replace PubChem NAs with LMDB pubchem IDs
All_names_hmdb_LMDB <- merge(All_names_hmdb, LMDB, by = "InChI.Key", all.x = TRUE)

for (i in 1:nrow(All_names_hmdb_LMDB)){
  if(is.na(All_names_hmdb_LMDB$PubChem.ID[i]) & !(is.na(All_names_hmdb_LMDB$PubChem.y[i]))){
    All_names_hmdb_LMDB$PubChem.ID[i] <- All_names_hmdb_LMDB$PubChem.y[i]
  }
}

#All_names <- merge(All_names, HMDB,by = "InChI.Key", all.x = TRUE)

#remove pubchem columns from HMDB and LMDB now that the info is all saved in the PubChem.ID column
All_names_hmdb_LMDB <- All_names_hmdb_LMDB[,!names(All_names_hmdb_LMDB) %in% c("PubChem.x", "PubChem.y")]
All_names_hmdb_LMDB_crab <- merge(unique(All_names_hmdb_LMDB), unique(crab[,c(4:6)]), by = "InChI.Key", all.x = TRUE)

#remove the compounds where the Crab pubchem ID does not match the HMDB/WCMC/LMDB PubChem ID:
#this code keeps all lines where there was no crab pubchem ID and lines where the crab pubchem ID matched the other pubchem ID
#All_names_hmdb_LMDB_crab <- All_names_hmdb_LMDB_crab[which(!(!(is.na(All_names_hmdb_LMDB_crab$PubChem.ID)) & All_names_hmdb_LMDB_crab$PubChem.ID!=All_names_hmdb_LMDB_crab$PubChem) | is.na(All_names_hmdb_LMDB_crab$PubChem)),]

#replace NAs in PubChem.ID column with crab pubchem values
for (i in 1:nrow(All_names_hmdb_LMDB_crab)){
  if(is.na(All_names_hmdb_LMDB_crab$PubChem.ID[i]) & !(is.na(All_names_hmdb_LMDB_crab$PubChem[i]))){
    All_names_hmdb_LMDB_crab$PubChem.ID[i] <- All_names_hmdb_LMDB_crab$PubChem[i]
  }
}

#Save a data frame of compounds with discrepancies between PubChem.ID and PubChem columns
pubchem_diff <- All_names_hmdb_LMDB_crab[which(All_names_hmdb_LMDB_crab$PubChem.ID != All_names_hmdb_LMDB_crab$PubChem & !(is.na(All_names_hmdb_LMDB_crab$PubChem))),c("InChI.Key","PubChem.ID", "BinBase.name", "PubChem")]

#go through each line of the main data frame and substitute the PubChem.ID value for the PubChem value if the inchikey is in the discrepancy df
# And doesn't contain the patterns LPE, PC (36 , and PC (38 in the BinBase.name
for (i in 1:nrow(All_names_hmdb_LMDB_crab)){
  if(All_names_hmdb_LMDB_crab$InChI.Key[i] %in% pubchem_diff[-grep("LPE|PC \\(36|PC \\(38", pubchem_diff$BinBase.name),"InChI.Key"]){
    All_names_hmdb_LMDB_crab$PubChem.ID[i] <- All_names_hmdb_LMDB_crab$PubChem[i]
  }
}


#replace empty KEGG fields with crab keggs

All_names_hmdb_LMDB_crab$KEGG.y <- gsub("not found", NA, All_names_hmdb_LMDB_crab$KEGG.y)

for (i in 1:nrow(All_names_hmdb_LMDB_crab)){
  if(All_names_hmdb_LMDB_crab$KEGG.x[i]==""){
    All_names_hmdb_LMDB_crab$KEGG.x[i] <- as.character(All_names_hmdb_LMDB_crab$KEGG.y[i])
  }
}

#remove extra columns
All_names_hmdb_LMDB_crab <- All_names_hmdb_LMDB_crab[,!names(All_names_hmdb_LMDB_crab) %in% c("PubChem", "KEGG.y", "LM.ID","Creation.Date")]

#change binbase name to annotation

colnames(All_names_hmdb_LMDB_crab)[8] <- "Annotation"

#replace 'not found' with NA
All_names_hmdb_LMDB_crab$KEGG <- gsub("not found",NA, All_names_hmdb_LMDB_crab$KEGG)


#Find compounds with different KEGG IDs
All_names_hmdb_LMDB_crab[which(All_names_hmdb_LMDB_crab$KEGG.ID!=All_names_hmdb_LMDB_crab$KEGG.x),c(6,9,13)]
#there are 5 IDs that come up, one that has no value and 4 that are different:
# C03547 : omega-Hydroxy fatty acid; other columns have C00160 Glycolate; should be glycolate
# C00134 : putrescine; other columns have C00138 ferredoxin; should be putrescine as annotated by WCMC
# C01991: 4-Hydroxyacid; other columns have 4-Hydroxybutanoate; compounds look identical and 4-Hydroxybutanoate has more biological info so keep that one 
# C00539: aromatic acid; other columns have C00180 benzoate; benzoate is the more specific name with bio annotations so use that one
# So only putrescine needs to be swapped
All_names_hmdb_LMDB_crab[which(All_names_hmdb_LMDB_crab$Annotation == "putrescine"),c("KEGG.x","KEGG")] <- "C00134"
#now remove the KEGG.ID column
All_names_hmdb_LMDB_crab <- All_names_hmdb_LMDB_crab[,-grep("KEGG.ID",colnames(All_names_hmdb_LMDB_crab))]

#now check other two KEGG columns
All_names_hmdb_LMDB_crab[which(All_names_hmdb_LMDB_crab$KEGG.x!=All_names_hmdb_LMDB_crab$KEGG),c(8,12)]

#there are no discrepancies so only keep one kegg column
All_names_hmdb_LMDB_crab <- All_names_hmdb_LMDB_crab[,-grep("KEGG.x",colnames(All_names_hmdb_LMDB_crab))]

#fill in empty fields with NAs
All_names_hmdb_LMDB_crab[All_names_hmdb_LMDB_crab ==""] <- NA

#See if CTS data shows different KEGG and CIDs

#change no result to NA
CTS[CTS=="No result"] <- NA

#remove duplicates
CTS <- unique(CTS)

#remove rows with NA in KEGG and PubChem columns
CTS <- CTS[which(!(is.na(CTS$KEGG_CTS)) & substr(CTS$KEGG_CTS,1,1)!="D" | !(is.na(CTS$PubChem.CID))),]

#combine CTS data with all data
All_names_hmdb_LMDB_crab_CTS <- merge(All_names_hmdb_LMDB_crab,CTS, by = "InChI.Key", all.x = T)

#preview missing pubchem IDs and kegg ids that are in the CTS data
View(All_names_hmdb_LMDB_crab_CTS[which(is.na(All_names_hmdb_LMDB_crab_CTS$KEGG)& !(is.na(All_names_hmdb_LMDB_crab_CTS$KEGG_CTS))),c("Annotation", "InChI.Key", "PubChem.ID", "PubChem.CID","KEGG", "KEGG_CTS")])

#loop through and replace NAs in PubChem.ID and KEGG column with CTS data
for (i in 1:nrow(All_names_hmdb_LMDB_crab_CTS)){
  if(is.na(All_names_hmdb_LMDB_crab_CTS$PubChem.ID[i]) & !(is.na(All_names_hmdb_LMDB_crab_CTS$PubChem.CID[i]))){
    All_names_hmdb_LMDB_crab_CTS$PubChem.ID[i] <- All_names_hmdb_LMDB_crab_CTS$PubChem.CID[i]
  }
  if(is.na(All_names_hmdb_LMDB_crab_CTS$KEGG[i]) & !(is.na(All_names_hmdb_LMDB_crab_CTS$KEGG_CTS[i]))){
    All_names_hmdb_LMDB_crab_CTS$KEGG[i] <- All_names_hmdb_LMDB_crab_CTS$KEGG_CTS[i]
  }
}

#confirm this shows nothing now
View(All_names_hmdb_LMDB_crab_CTS[which(is.na(All_names_hmdb_LMDB_crab_CTS$KEGG)& !(is.na(All_names_hmdb_LMDB_crab_CTS$KEGG_CTS))),c("Annotation", "InChI.Key", "PubChem.ID", "PubChem.CID","KEGG", "KEGG_CTS")])
#it does

#remove CTS columns
All_names_hmdb_LMDB_crab_CTS <- All_names_hmdb_LMDB_crab_CTS[,!names(All_names_hmdb_LMDB_crab_CTS) %in% c("KEGG_CTS", "PubChem.CID")]

#need to convert inosine's pubchem ID
All_names_hmdb_LMDB_crab_CTS[which(All_names_hmdb_LMDB_crab_CTS_SMILES$Annotation == "inosine"),"PubChem.ID"] <- "135398641"



#rename columns
#colnames(All_names) <- c("InChI.Key", "BinBase.name", "PubChem", "KEGG", "HMDBpc", "LMDBpc", "CrabPC", "CrabKegg")


#remove duplicate lines with same pubchem IDs
#All_names_STACKED <-tidyr::gather(All_names, "source", "ID", c(3,5,6,7))
#All_names_STACKED <- All_names_STACKED[,-5]
#All_names_STACKED <- unique(All_names_STACKED)

#rename column back to pubchem
#colnames(All_names_STACKED)[5] <- "PubChem"
#convert back to wide format excluding KEGG ids and pubchem IDs = NA
#All_names_uniquePubChem <- All_names_STACKED[which(!(is.na(All_names_STACKED$PubChem))),]



##compounds are sometimes detected in both lipidomics and metabolomics and don't have the same binbase name, but do have the same pubchem ID
#these should probably be considered only once...


dup_pubchem <- All_names_hmdb_LMDB_crab_CTS[duplicated(All_names_hmdb_LMDB_crab_CTS$PubChem.ID),"PubChem.ID"]
View(All_names_hmdb_LMDB_crab_CTS[which(All_names_hmdb_LMDB_crab_CTS$PubChem.ID %in% dup_pubchem),])

#for Fatty Acids, use only the compounds detected in lipidomics
#for PCs, use 36:4A, 36:5B, 38:6B, 38:5A

All_names_hmdb_LMDB_crab_CTS <- unique(All_names_hmdb_LMDB_crab_CTS[-grep("heptadecanoic acid|arachidic acid|pentadecanoic acid|arachidonic acid|palmitoleic acid|oleic acid|stearic acid|cis-gondoic acid|behenic acid|PC \\(38:5\\) B|PC \\(38:6\\)|PC \\(38:6\\) A|	
PC \\(38:6\\) C|PC \\(36:4\\) B|PC \\(36:5\\) A",All_names_hmdb_LMDB_crab_CTS$Annotation),])


#write out table
write.table(All_names_hmdb_LMDB_crab_CTS, "data/All_Known_Compound_Annotations.txt", sep ="\t", quote = F, row.names = F)

#go to PubChem Identifier Exchange and convert pubchem IDs to SMILES
#https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi

#read in table
pubchem_SMILES <- read.table("data/237110758970980098.txt", sep ="\t", stringsAsFactors = FALSE,header = F, quote ="",comment.char = "")

colnames(pubchem_SMILES) <- c("PubChem.ID", "SMILES")

#merge SMILES with data
All_names_hmdb_LMDB_crab_CTS_SMILES <- merge(All_names_hmdb_LMDB_crab_CTS, pubchem_SMILES, by = "PubChem.ID", all.x = T)

#for some reason inosine and guanosine are missing their smiles so I looked them up manually,

#guanosine C1=NC2=C(N1C3C(C(C(O3)CO)O)O)NC(=NC2=O)N
#inosince C1=NC2=C(N1C3C(C(C(O3)CO)O)O)NC(=NC2=O)N

# I updated these in these in the file 237110758970980098.txt

#realized inosine had an incorrect pubchem ID so I updated it in the CTS file 
All_names_hmdb_LMDB_crab_CTS[which(All_names_hmdb_LMDB_crab_CTS$InChI.Key=="UGQMRVRMYYASKQ-KQYNXXCUSA-N"),"PubChem.ID"] <- "135398641"
All_names_hmdb_LMDB_crab_CTS[which(All_names_hmdb_LMDB_crab_CTS$InChI.Key=="RHGKLRLOHDJJDR-BYPYZUCNSA-N"),"PubChem.ID"] <- "9750"
  


#read in table
pubchem_SMILES <- read.table("data/237110758970980098.txt", sep ="\t", stringsAsFactors = FALSE,header = F, quote ="",comment.char = "")

colnames(pubchem_SMILES) <- c("PubChem.ID", "SMILES")

#merge SMILES with data
All_names_hmdb_LMDB_crab_CTS_SMILES <- merge(All_names_hmdb_LMDB_crab_CTS, pubchem_SMILES, by = "PubChem.ID", all.x = T)


View(All_names_hmdb_LMDB_crab_CTS_SMILES[duplicated(All_names_hmdb_LMDB_crab_CTS_SMILES$SMILES),])
View(All_names_hmdb_LMDB_crab_CTS_SMILES[duplicated(All_names_hmdb_LMDB_crab_CTS_SMILES$PubChem.ID),])
View(All_names_hmdb_LMDB_crab_CTS_SMILES[duplicated(All_names_hmdb_LMDB_crab_CTS_SMILES$Annotation),])

#There are duplicate annotations; this is because of adducts
#assigning unique names for adducts
All_names_hmdb_LMDB_crab_CTS_SMILES[which(All_names_hmdb_LMDB_crab_CTS_SMILES$InChI.Key=="IIZPXYDJLKNOIY-JXPKJXOSSA-N"),"Annotation"] <- "PC (36:4) A..M.Hac.H"
All_names_hmdb_LMDB_crab_CTS_SMILES[which(All_names_hmdb_LMDB_crab_CTS_SMILES$InChI.Key=="CNNSEHUKQJCGTE-UPPWDXJYSA-N"),"Annotation"] <- "PC (34:3) M.Hac.H"
All_names_hmdb_LMDB_crab_CTS_SMILES[which(All_names_hmdb_LMDB_crab_CTS_SMILES$InChI.Key=="VRBJSZDBPQLNEM-PFTGJCAASA-N"),"Annotation"] <- "PE (38:6) M-H"
All_names_hmdb_LMDB_crab_CTS_SMILES[which(All_names_hmdb_LMDB_crab_CTS_SMILES$InChI.Key=="YLWBKBDNHWQEFU-YJXJLLHLSA-N"),"Annotation"] <- "PC (38:5) A..M.Hac.H"
All_names_hmdb_LMDB_crab_CTS_SMILES[which(All_names_hmdb_LMDB_crab_CTS_SMILES$InChI.Key=="KCNBSSYOJRUKOM-JALPNNRCSA-N"),"Annotation"] <- "PE (p-38:5) or PE (o-38:6) M-H"
All_names_hmdb_LMDB_crab_CTS_SMILES[which(All_names_hmdb_LMDB_crab_CTS_SMILES$InChI.Key=="ATTCDOPAYPGSLE-LQULQHAGSA-N"),"Annotation"] <- "PC (p-38:5) or PC (o-38:6) M.Hac.H"
All_names_hmdb_LMDB_crab_CTS_SMILES[which(All_names_hmdb_LMDB_crab_CTS_SMILES$InChI.Key=="DZPHTQGGRSWHLG-USOCFQKQSA-N"),"Annotation"] <- "PE (p-40:5) or PE (o-40:6) M-H"
All_names_hmdb_LMDB_crab_CTS_SMILES[which(All_names_hmdb_LMDB_crab_CTS_SMILES$InChI.Key=="XVXISDREVDGQPX-PISDLAQISA-N"),"Annotation"] <- "PE (p-36:2) or PE (o-36:3) M-H"

#check names are unique
View(All_names_hmdb_LMDB_crab_CTS_SMILES[duplicated(All_names_hmdb_LMDB_crab_CTS_SMILES$Annotation),])

#they are except for the ones with no info
#only keep rows with values in PubChem ID and and SMILES

All_names_hmdb_LMDB_crab_CTS_SMILES_metamapp <- All_names_hmdb_LMDB_crab_CTS_SMILES[which(!(is.na(All_names_hmdb_LMDB_crab_CTS_SMILES$PubChem.ID))& !(is.na(All_names_hmdb_LMDB_crab_CTS_SMILES$SMILES))),]
All_names_hmdb_LMDB_crab_CTS_SMILES_metamapp <- All_names_hmdb_LMDB_crab_CTS_SMILES_metamapp[,c("PubChem.ID","KEGG", "SMILES", "Annotation")]









#prepare metamapp input

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

########################################################
####################### STEP 4 #########################
########################################################

#Combine chemical identifiers with stats data

###Read in stats data

library(readxl)
metab_uni<- read_xlsx("~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/pteropod_genmetabANDlipids_univariate_RESULTSsupplement.xlsx", sheet =1)
lipid_uni<- read_xlsx("~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/pteropod_genmetabANDlipids_univariate_RESULTSsupplement.xlsx", sheet =2)

colnames(metab_uni)[1] <- "BinBasename_R"
colnames(lipid_uni)[1] <- "Identifier_R"


#read in ID keys
metab_ID <- read.table("../MetamappAnalysis/metab_ID_Key.txt", sep ="\t", stringsAsFactors = FALSE,header = TRUE, quote ="")
lipid_ID <- read.table("../MetamappAnalysis/lipid_ID_Key.txt", sep ="\t", stringsAsFactors = FALSE,header = TRUE, quote ="",comment.char = "")

#read in chemical identifier info
pubch_kegg_smiles <- read.csv("~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/unique_inchi_pubchem_kegg_smiles.csv", stringsAsFactors = FALSE)
colnames(pubch_kegg_smiles)[4] <- "SMILES"

#pubchem table still contains 1 duplicate where one line has a kegg and the other doesn't. It is pubchem ID 24699
#remove the line missing the kegg ID
pubch_kegg_smiles <- pubch_kegg_smiles[-which(pubch_kegg_smiles$PubChem=="24699" & is.na(pubch_kegg_smiles$KEGG)),]


#read in RSD, FC, and Cohen D info
RSD_FC_cohen <- read.table("../FoldChange_CohenD/all_cmpds_FC_RSD_cohen.tsv", sep = "\t",header = TRUE,stringsAsFactors = FALSE, quote = "", comment.char = "")


#change name of column 2
colnames(metab_ID_table)[2] <- "Identifier"




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

#merge univariate data with ID keys
metab_uni_names <- merge(metab_uni[,c(1,4,8,14,20)],metab_ID[,-c(3:4)], by = "BinBasename_R")
#reorder columns
metab_uni_names <- metab_uni_names[,c(1:5,7,6,8)]
#rename columns to match lipid data colnames
colnames(metab_uni_names)[1] <- "Identifier_R"
colnames(metab_uni_names)[6] <- "Identifier"
colnames(metab_uni_names)[7] <- "Annotation"

#merge univariate data with ID keys
lipid_uni_names <- merge(lipid_uni[,c(1,4,8,14,20)],lipid_ID, by = "Identifier_R")

#bind lipid and metab data into one df
all_analytes <- rbind(metab_uni_names,lipid_uni_names)


###Combine univariate data with chemical identifier data
#merge pubchem, kegg, and smiles IDs with all analytes by their InChI.Key
df <- merge(all_analytes,pubch_kegg_smiles, by = "InChI.Key")

#count the number of rows in df
nrow(df)
#There are 239 analytes that match to InChI.Keys

#Count the number of duplicate InChI.Keys
length(unique(df$InChi.Key))
#There are 23/239 InChi.Keys that are duplicated

## Resolve duplicated InChI.Keys
View(df[which(duplicated(df$InChI.Key)),])


# top_n depricated https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&cad=rja&uact=8&ved=2ahUKEwjx8Nmx2qboAhXLo54KHW0bAvEQFjACegQIAxAB&url=https%3A%2F%2Fgithub.com%2Ftidyverse%2Fdplyr%2Fblob%2Fmaster%2FR%2Ftop-n.R&usg=AOvVaw0OzD_NLLLtsi-F3o7AGlA_
# Using aggregate and max https://nsaunders.wordpress.com/2013/02/13/basic-r-rows-that-contain-the-maximum-value-of-a-variable/

df_agg <- aggregate(`Pr(>Chisq)_overall` ~ InChI.Key, df, min)
df <- merge(df_agg, df)

# Check for unique pubchem IDs

View(df[which(duplicated(df$PubChem)),])

# There is one compound with duplicated pubchem ID "PC..38.5..B" and "PC..38.5..A.1"
# select the compound with the lowest pvalue

df_agg <- aggregate(`Pr(>Chisq)_overall` ~ PubChem, df, min)
df <- merge(df_agg, df)

#Check for unique SMILES
nrow(df[which(duplicated(df$SMILES)),])
View(df[which(duplicated(df$SMILES)),])
# these are compounds that have no SMILES entry so just leave these




#merge RSD data with lm data
RSD_FC_cohen_pub <- merge(RSD_FC_cohen,df, by = "Identifier")

#create table for Metamapp
#metamapp_input <- RSD_FC_cohen_pub[,c("PubChem", "KEGG", "SMILES", "Annotation","Pr(>Chisq)_overall","logFC_HL")]
metamapp_input <- RSD_FC_cohen_pub[,c("PubChem", "KEGG", "SMILES", "Identifier_R","Pr(>Chisq)_overall","logFC_HL")]

colnames(metamapp_input) <- c("PubChem_ID", "KEGG_ID",	"SMILES",	"Compound_Name", "Condition_A_pvalue",	"Condition_A_foldchange")
## Thes column names below didn't work because they need need to match exactly with metamapp example
#colnames(df) <- c("PubChem_ID", "KEGG_ID", "SMILES", "Compound_Name", "Chi_overall_pvalue")

#remove rows without SMILES because they are problematic in uploading the data
write.table(metamapp_input, "~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/metamapp_input.csv", sep = "\t",quote = FALSE, row.names = FALSE)


#format file for chemRICH input
chemrich_input <- RSD_FC_cohen_pub[,c("Identifier_R","InChI.Key.y","PubChem", "SMILES","Pr(>Chisq)_overall","logFC_HL")]
colnames(chemrich_input) <- c("Compound Name", "InChiKeys",	"Pubchem ID",	"SMILES", "pvalue",	"foldchange")
write.table(chemrich_input, "~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/chemrich_input.csv", sep = "\t",quote = FALSE, row.names = FALSE)

#for DO (20200824)
### manually added the SMILES IDs from PUBCHEM database to guanosine, inosine, and X2.hydroxypyrazinyl.2.propenoic.acid.ethyl.ester.NIST

# X2.hydroxypyrazinyl.2.propenoic.acid.ethyl.ester.NIST, CCOC(=O)C(=O)C=C1C=NC=CN1
# guanosine, C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)NC(=NC2=O)N
# inosine, C1=NC(=O)C2=C(N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O

RSD_FC_cohen_pub[which(RSD_FC_cohen_pub$Annotation == "guanosine"),"SMILES"] <- "C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)NC(=NC2=O)N"
RSD_FC_cohen_pub[which(RSD_FC_cohen_pub$Annotation == "inosine"),"SMILES"] <- "C1=NC(=O)C2=C(N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O"
RSD_FC_cohen_pub[which(RSD_FC_cohen_pub$Annotation == "2-hydroxypyrazinyl-2-propenoic acid ethyl ester NIST"),"SMILES"] <- "CCOC(=O)C(=O)C=C1C=NC=CN1"



chemrich_input_DO <- RSD_FC_cohen_pub[,c("Identifier_R","InChI.Key.y","PubChem", "SMILES","DO_p.value","logFC_HL")]
chemrich_input_pH <- RSD_FC_cohen_pub[,c("Identifier_R","InChI.Key.y","PubChem", "SMILES","pH_p.value","logFC_LH")]
chemrich_input_pH_DO <- RSD_FC_cohen_pub[,c("Identifier_R","InChI.Key.y","PubChem", "SMILES","pH:DO_p.value","logFC_LL")]

colnames(chemrich_input_DO) <- c("Compound Name", "InChiKeys",	"Pubchem ID",	"SMILES", "pvalue",	"foldchange")
colnames(chemrich_input_pH) <- c("Compound Name", "InChiKeys",	"Pubchem ID",	"SMILES", "pvalue",	"foldchange")
colnames(chemrich_input_pH_DO) <- c("Compound Name", "InChiKeys",	"Pubchem ID",	"SMILES", "pvalue",	"foldchange")


write.table(chemrich_input_DO, "~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/chemrich_input_HL.csv", sep = "\t",quote = FALSE, row.names = FALSE)
write.table(chemrich_input_pH, "~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/chemrich_input_LH.csv", sep = "\t",quote = FALSE, row.names = FALSE)
write.table(chemrich_input_pH_DO, "~/Documents/GitHub/pteropod_pHxDO_metabolomics/MetamappAnalysis/data/chemrich_input_LL.csv", sep = "\t",quote = FALSE, row.names = FALSE)



########################################################
####################### STEP 5 #########################
########################################################

## Run MetaMapp

# Imported .tsv file into excel
# removed rows that had no SMILES IDs (there were 3 rows removed)
# Copied and pasted table to MetaMapp web interface
# Downloaded .sif and .tsv files
# moved files into project directory

########################################################
####################### STEP 6 #########################
########################################################

# Add more stats data to node attribute file

#read in MetaMapp node attribute output file
node_attr <- read.table("data/node_attributes_chemsim_krp_07-1.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE)

#merge with RSD, fold change, and cohen D data
node_attr <- merge(node_attr[,1:4], RSD_FC_cohen_pub[,c(3:20,22,25:28,30)], by ="SMILES")

#create column denoting lm significance
node_attr$sig <- "1"
for (i in 1:nrow(node_attr)){
  if(node_attr$`Pr(>Chisq)_overall`[i] < 0.05){
    node_attr$sig[i] <- "5"
  }
}

write.table(node_attr, "data/node_attributes_chemsim_krp_07-1.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


#############################

#get pubchem IDs for AOV and PLSDA sig compounds

#read in data 
PLSDA_aov_cmpds <- read.csv("data/PLSDA_AOV_sig_cmpds.csv", stringsAsFactors = F)


#merge with pubchem IDs

pubchem_PLSDA_aov_cmpds <- merge(RSD_FC_cohen_pub,PLSDA_aov_cmpds[,c(1,5,24:30)], by = "Identifier")

#order by sig compounds
pubchem_PLSDA_aov_cmpds <- pubchem_PLSDA_aov_cmpds[order(pubchem_PLSDA_aov_cmpds$method),]

#wrtie out file so pubchem IDs can be uplodaded to metaboAnalyst enrichment analysis
write.csv(pubchem_PLSDA_aov_cmpds, "pubchem_PLSDA_aov_cmpds.csv", quote = F, row.names = F)

####################################################
############### OLD ANALYSIS #######################
####################################################
###METAMAPP
#copied and pasted arbitrary FC values and uploaded to metamapp

##need to resolve duplicated compound names (see email from Dinesh)

#print duplicated compound names
#dups <- df[which(duplicated(df$Compound_Name)),"Compound_Name"]
#dups <- df[which(duplicated(df$PubChem)),"PubChem"]

#View(df[which(df$PubChem %in% dups),])

#removed duplicates manually and tried to keep the ones that had KEGG IDs

#went through and manually removed duplicated pubchem IDs (most of these were from compounds detected in both GC/MS and LC/MS with different compound names but same pubchem IDs)
#saved this file as metamapp_input_nodups.csv

#copied and pasted header from Metamapp example data into metamapp_input_nodups.csv
#Uploaded to MetaMapp and downloaded sif and attributes file
