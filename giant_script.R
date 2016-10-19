# The giant script consolidating codes from start to finish without the extra codes that were written for testing.
# Redundant comments are also removed.

################################################### START ####################################################
# Processing all gene data in the ExAC dataset for use in downstream scripts

setwd("~/ExAC dataset project/csv_records")
# Load libraries
library(plyr)
library(dplyr)
library(biomaRt)
library(goseq)
marthuman = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
martchimp = useEnsembl(biomart="ensembl", dataset="ptroglodytes_gene_ensembl", GRCh=37)
martorang = useEnsembl(biomart="ensembl", dataset="pabelii_gene_ensembl", GRCh=37)

exacgenes <- read.table("exac_r03.txt", header = TRUE, sep="")
write.csv(exacgenes, "allexacgenes.csv", row.names=FALSE) # Keeping a .csv record

exacgenes <- read.csv("allexacgenes.csv", header = TRUE, sep =",")
# Turn the transcript column into characters so that the next line will work
exacgenes$transcript <- as.character(exacgenes$transcript)
# Removing the ".x" from the ends of the transcript IDs
exacgenes$transcript <- sub("\\..*", "\\", exacgenes$transcript)
# Only keep the transcript, gene and pLI columns
exacgenes <- subset(exacgenes, select=c("transcript", "gene", "pLI"))

# Get ensembl gene IDs
geneid <- getBM(attributes = c("ensembl_gene_id"), filters = "ensembl_transcript_id", values = exacgenes$transcript, mart = marthuman, uniqueRows = FALSE)
# Append to original data frame and rename columns
exacgenes <- cbind(geneid, exacgenes)
colnames(exacgenes) <- c("Ensembl_gene_ID", "Ensembl_transcript_ID", "Gene_name", "pLI")
# Order by ascending values of pLI
exacgenes <- exacgenes[order(exacgenes$pLI),]
head(exacgenes)
write.csv(exacgenes, "exacgenes.csv", row.names=FALSE)

############################################ Getting dN/dS values ############################################

allexacgenes <- read.csv("exacgenes.csv", header = TRUE, sep =",")

# Filter out all human ribosomal protein genes
geneFamilies <- getBM(attributes=c("ensembl_gene_id", "family", "family_description", "external_gene_name"), mart=marthuman)
badRibosomes <- geneFamilies[grepl("MRP|RPS|RPL|RP L|RP S", geneFamilies$family_description),1]
moreBadRibos <- geneFamilies[grepl("MRP|RPS|RPL|RP L|RP S", geneFamilies$external_gene_name),1]
allBadRibos <- c(badRibosomes, moreBadRibos)
# Record ribo IDs that will be removed from dataset
Genes_lost_humanribos <- allexacgenes[allexacgenes$Ensembl_gene_ID %in% allBadRibos,]
write.csv(Genes_lost_humanribos, "Lost_humanribos.csv", row.names=FALSE) # 206
# Remove ribosomal protein genes from list of genes
geneidx <- allexacgenes[!allexacgenes$Ensembl_gene_ID %in% allBadRibos,]

# Extracting dn and ds values
# Human vs chimp
dndshumchimp <- getBM(attributes = c("ensembl_gene_id", "ptroglodytes_homolog_dn", "ptroglodytes_homolog_ds"),
                      filters = "ensembl_gene_id", values = geneidx$Ensembl_gene_ID, mart = marthuman)
colnames(dndshumchimp) <- c("Ensembl_gene_ID", "dnhumchimp", "dshumchimp")
humchimpdndsr <- c(dndshumchimp$dnhumchimp / dndshumchimp$dshumchimp)
dndshumchimp <- cbind(dndshumchimp, humchimpdndsr)
write.csv(dndshumchimp, "dndsHC.csv", row.names=FALSE)

# Human vs orangutan
dndshumorang <- getBM(attributes = c("ensembl_gene_id", "pabelii_homolog_dn", "pabelii_homolog_ds"),
                      filters = "ensembl_gene_id", values = geneidx$Ensembl_gene_ID, mart = marthuman)
colnames(dndshumorang) <- c("Ensembl_gene_ID", "dnhumorang", "dshumorang")
humorangdndsr <- c(dndshumorang$dnhumorang / dndshumorang$dshumorang)
dndshumorang <- cbind(dndshumorang, humorangdndsr)
write.csv(dndshumorang, "dndsHO.csv", row.names=FALSE)

# Human vs macaca
dndshummac <- getBM(attributes = c("ensembl_gene_id", "mmulatta_homolog_dn", "mmulatta_homolog_ds"),
                    filters = "ensembl_gene_id", values = geneidx$Ensembl_gene_ID, mart = marthuman)
colnames(dndshummac) <- c("Ensembl_gene_ID", "dnhummac", "dshummac")
hummacdndsr <- c(dndshummac$dnhummac / dndshummac$dshummac)
dndshummac <- cbind(dndshummac, hummacdndsr)
write.csv(dndshummac, "dndsHM.csv", row.names=FALSE)

# To get dN and dS values for chimp comparisons with other species, 
# need to convert to chimp IDs to avoid problems with converting IDs between multiple species orthologues
# Filter out all chimp ribosomal proteins
geneFamilies <- getBM(attributes=c("ensembl_gene_id", "family", "family_description", "external_gene_name"), mart=martchimp)
badCRibosomes <- geneFamilies[grepl("MRP|RPS|RPL|RP L|RP S", geneFamilies$family_description),1]
moreBadCRibos <- geneFamilies[grepl("MRP|RPS|RPL|RP L|RP S", geneFamilies$external_gene_name),1]
allBadCRibos <- c(badCRibosomes, moreBadCRibos) # These are in chimp IDs

# Convert human ensembl gene IDs to chimpanzee IDS
chimpgeneid <- getBM(attributes=c("ptroglodytes_homolog_ensembl_gene", "ensembl_gene_id", "ptroglodytes_homolog_orthology_type"),
                     filters = "ensembl_gene_id", values = geneidx$Ensembl_gene_ID, mart = marthuman, uniqueRows = FALSE)
chimpgeneidx <- chimpgeneid[!chimpgeneid$ptroglodytes_homolog_ensembl_gene %in% allBadCRibos,]
Genes_lost_chimpribos <- chimpgeneid[chimpgeneid$ptroglodytes_homolog_ensembl_gene %in% allBadCRibos,]
write.csv(Genes_lost_chimpribos, "Lost_chimpribos.csv", row.names=FALSE) # 1
# Not all human IDs can be converted to chimp IDs. Change all blanks to NA for kicking out.
chimpgeneidx[chimpgeneidx == ""] <- NA
Genes_lost_chimpNA <- chimpgeneidx[is.na(chimpgeneidx$ptroglodytes_homolog_ensembl_gene), ]
write.csv(Genes_lost_chimpNA, "Lost_chimpIDNA.csv", row.names=FALSE) # 1194
chimpgeneidy <- chimpgeneidx[!is.na(chimpgeneidx$ptroglodytes_homolog_ensembl_gene), ]
# Some human IDs give multiple chimp IDs. Check for repeated chimp IDs.
chimpdup <- chimpgeneidy[(duplicated(chimpgeneidy$ptroglodytes_homolog_ensembl_gene) == TRUE),]
# Record the repeats which will be kicked out. These are many human : 1 chimp IDs.
chimpduprec <- chimpgeneidy[(chimpgeneidy$ptroglodytes_homolog_ensembl_gene %in% chimpdup$ptroglodytes_homolog_ensembl_gene), ]
write.csv(chimpduprec, "Lost_many2one_chimp.csv", row.names=FALSE) # 408
chimpgeneidy[(chimpgeneidy$ptroglodytes_homolog_ensembl_gene %in% chimpdup$ptroglodytes_homolog_ensembl_gene), ] <- NA
chimpgeneidy <- chimpgeneidy[!is.na(chimpgeneidy$ptroglodytes_homolog_ensembl_gene), ]

# Chimp vs orangutan
dndschimporang <- getBM(attributes = c("ensembl_gene_id", "pabelii_homolog_dn", "pabelii_homolog_ds"),
                        filters = "ensembl_gene_id", values = chimpgeneidy, mart = martchimp)
colnames(dndschimporang) <- c("Ensembl_gene_ID_chimp", "dnchimporang", "dschimporang")
chimporangdndsr <- c(dndschimporang$dnchimporang / dndschimporang$dschimporang)
dndschimporang <- cbind(dndschimporang, chimporangdndsr)
write.csv(dndschimporang, "dndsCO.csv", row.names=FALSE)

# Chimp vs macaca
dndschimpmac <- getBM(attributes = c("ensembl_gene_id", "mmulatta_homolog_dn", "mmulatta_homolog_ds"),
                      filters = "ensembl_gene_id", values = chimpgeneidy, mart = martchimp)
colnames(dndschimpmac) <- c("Ensembl_gene_ID_chimp", "dnchimpmac", "dschimpmac")
chimpmacdndsr <- c(dndschimpmac$dnchimpmac / dndschimpmac$dschimpmac)
dndschimpmac <- cbind(dndschimpmac, chimpmacdndsr)
summary(dndschimpmac)
write.csv(dndschimpmac, "dndsCM.csv", row.names=FALSE)

# Convert human ensembl gene IDs to orangutan IDs
oranggeneid <- getBM(attributes=c("pabelii_homolog_ensembl_gene", "ensembl_gene_id", "pabelii_homolog_orthology_type"),
                     filters = "ensembl_gene_id", values = geneidx$Ensembl_gene_ID, mart = marthuman, uniqueRows = FALSE)
geneFamilies <- getBM(attributes=c("ensembl_gene_id", "family", "family_description", "external_gene_name"), mart=martorang)
badORibosomes <- geneFamilies[grepl("MRP|RPS|RPL|RP L|RP S", geneFamilies$family_description),1]
moreBadORibos <- geneFamilies[grepl("MRP|RPS|RPL|RP L|RP S", geneFamilies$external_gene_name),1]
allBadORibos <- c(badORibosomes, moreBadORibos) # These are in orangutan IDs
Genes_lost_orangribos <- oranggeneid[oranggeneid$pabelii_homolog_ensembl_gene %in% allBadORibos,]
write.csv(Genes_lost_orangribos, "Lost_orangribos.csv", row.names=FALSE)
oranggeneidx <- oranggeneid[!oranggeneid$pabelii_homolog_ensembl_gene %in% allBadORibos,] # No genes lost here
oranggeneidx[oranggeneidx == ""] <- NA
Genes_lost_orangNA <- oranggeneidx[is.na(oranggeneidx$pabelii_homolog_ensembl_gene), ]
write.csv(Genes_lost_orangNA, "Lost_orangIDNA.csv", row.names=FALSE) # 1250
oranggeneidy <- oranggeneidx[!is.na(oranggeneidx$pabelii_homolog_ensembl_gene), ]
orangdup <- oranggeneidy[(duplicated(oranggeneidy$pabelii_homolog_ensembl_gene) == TRUE),]
orangduprec <- oranggeneidy[(oranggeneidy$pabelii_homolog_ensembl_gene %in% orangdup$pabelii_homolog_ensembl_gene), ] # 703
write.csv(orangduprec, "Lost_many2one_orang.csv", row.names=FALSE)
oranggeneidy[(oranggeneidy$pabelii_homolog_ensembl_gene %in% orangdup$pabelii_homolog_ensembl_gene), ] <- NA
oranggeneidy <- oranggeneidy[!is.na(oranggeneidy$pabelii_homolog_ensembl_gene), ]

#Orangutan vs macaca
dndsorangmac <- getBM(attributes = c("ensembl_gene_id", "mmulatta_homolog_dn", "mmulatta_homolog_ds"),
                      filters = "ensembl_gene_id", values = oranggeneidy, mart = martorang)
colnames(dndsorangmac) <- c("Ensembl_gene_ID_orang", "dnorangmac", "dsorangmac")
orangmacdndsr <- c(dndsorangmac$dnorangmac / dndsorangmac$dsorangmac)
dndsorangmac <- cbind(dndsorangmac, orangmacdndsr)
write.csv(dndsorangmac, "dndsOM.csv", row.names=FALSE)

# Records of human IDs with the corresponding chimp or orang IDs
write.csv(chimpgeneidy, "Gene_IDs_Human_Chimp.csv", row.names=FALSE)
write.csv(oranggeneidy, "Gene_IDs_Human_Orang.csv", row.names=FALSE)

########################################### Dataframe manipulations ##########################################

# Removing all rows with NA in dndsr column. Remove dN/dS values >10 (gets rid of Inf values).
# Rearranging data frames by ascending dn/ds ratios.
# All comparisons with gorilla dropped because all dN/dS values in human vs gorilla are NA.
# Script starts with .csv loading to allow bypassing of the whole dN/dS value retrival part above.

# Human vs chimp
clnhumchimp <- read.csv("dndsHC.csv", header=TRUE, sep=",") # 18225
# Record NA's
geneslost <- clnhumchimp[is.na(clnhumchimp$humchimpdndsr), ] # 1721
write.csv(geneslost, "Lost_dndsNA_HC.csv", row.names=FALSE)
clnhumchimp <- clnhumchimp[!is.na(clnhumchimp$humchimpdndsr), ] # 16504
lost <- subset(clnhumchimp, clnhumchimp$humchimpdndsr >10) # 557
write.csv(lost, "Lost_highcutHC.csv", row.names=FALSE)
clnhumchimp <- subset(clnhumchimp, !clnhumchimp$humchimpdndsr >10) 
clnhumchimp <- clnhumchimp[order(clnhumchimp$humchimpdndsr), ]
write.csv(clnhumchimp, "clnhumchimp.csv", row.names=FALSE) # 15947

# Human vs orangutan
clnhumorang <- read.csv("dndsHO.csv", header=TRUE, sep=",") # 18495
geneslost <- clnhumorang[is.na(clnhumorang$humorangdndsr), ] # 1330
write.csv(geneslost, "Lost_dndsNA_HO.csv, row.names=FALSE")
clnhumorang <- clnhumorang[!is.na(clnhumorang$humorangdndsr), ] # 17165
lost <- subset(clnhumorang, clnhumorang$humorangdndsr >10) # 66
write.csv(lost, "Lost_highcutHO.csv", row.names=FALSE)
clnhumorang <- subset(clnhumorang, !clnhumorang$humorangdndsr >10) 
clnhumorang <- clnhumorang[order(clnhumorang$humorangdndsr), ]
write.csv(clnhumorang, "clnhumorang.csv", row.names=FALSE) # 17009

# Human vs macaca
clnhummac <- read.csv("dndsHM.csv", header=TRUE, sep=",") # 19007
geneslost <- clnhummac[is.na(clnhummac$hummacdndsr), ] # 1192
write.csv(geneslost, "Lost_dndsNA_HM.csv", row.names=FALSE)
clnhummac <- clnhummac[!is.na(clnhummac$hummacdndsr), ] # 17815
lost <- subset(clnhummac, clnhummac$hummacdndsr >10) # 21
write.csv(lost, "Lost_highcutHM.csv", row.names=FALSE)
clnhummac <- subset(clnhummac, !clnhummac$hummacdndsr >10) 
clnhummac <- clnhummac[order(clnhummac$hummacdndsr), ] 
write.csv(clnhummac, "clnhummac.csv", row.names=FALSE) # 17794

# Chimp vs orangutan
clnchimporang <- read.csv("dndsCO.csv", header=TRUE, sep=",") # 16970
geneslost <- clnchimporang[is.na(clnchimporang$chimporangdndsr),] # 895
write.csv(geneslost, "Lost_dndsNA_CO.csv", row.names=FALSE)
clnchimporang <- clnchimporang[!is.na(clnchimporang$chimporangdndsr), ] # 16075
lost <- subset(clnchimporang, clnchimporang$chimporangdndsr >10) # 80
write.csv(lost, "Lost_highcutCO.csv", row.names=FALSE)
clnchimporang <- subset(clnchimporang, !clnchimporang$chimporangdndsr >10) 
clnchimporang <- clnchimporang[order(clnchimporang$chimporangdndsr), ]
write.csv(clnchimporang, "clnchimporang.csv", row.names=FALSE) # 15995

# Chimp vs macaca
clnchimpmac <- read.csv("dndsCM.csv", header=TRUE, sep=",") # 17475
geneslost <- clnchimpmac[is.na(clnchimpmac$chimpmacdndsr), ] # 821
write.csv(geneslost, "Lost_dndsNA_CM.csv", row.names=FALSE) 
clnchimpmac <- clnchimpmac[!is.na(clnchimpmac$chimpmacdndsr), ] # 16654
lost <-  subset(clnchimpmac, clnchimpmac$chimpmacdndsr >10) # 18
write.csv(lost, "Lost_highcutCM.csv", row.names=FALSE)
clnchimpmac <- subset(clnchimpmac, !clnchimpmac$chimpmacdndsr >10)
clnchimpmac <- clnchimpmac[order(clnchimpmac$chimpmacdndsr), ]
write.csv(clnchimpmac, "clnchimpmac.csv", row.names=FALSE) # 16636

# Orangutan vs macaca
clnorangmac <- read.csv("dndsOM.csv", header=TRUE, sep=",") # 15995
geneslost <- clnorangmac[is.na(clnorangmac$orangmacdndsr), ] # 708
write.csv(geneslost, "Lost_dndsNA_OM.csv", row.names=FALSE)
clnorangmac <- clnorangmac[!is.na(clnorangmac$orangmacdndsr), ] # 16683
lost <- subset(clnorangmac, clnorangmac$orangmacdndsr >10) # 31
write.csv(lost, "Lost_highcutOM.csv", row.names=FALSE)
clnorangmac <- subset(clnorangmac, !clnorangmac$orangmacdndsr >10) 
clnorangmac <- clnorangmac[order(clnorangmac$orangmacdndsr), ]
write.csv(clnorangmac, "clnorangmac.csv", row.names=FALSE) # 16652

############################################# Consolidate repeats ############################################

# Getting rid of the ID repeats that plague my life
# Remove repeated gene IDs from data frames by consolidating dN/dS values and taking means
# Reading in .csv files are for bypassing the above section.

# Read in Human and chimp ID list
IDhumanchimp <- read.csv("Gene_IDs_Human_Chimp.csv", header = TRUE, sep = ",")
IDhumanchimp$ptroglodytes_homolog_orthology_type <- NULL # Remove orthology type
colnames(IDhumanchimp) <- c("Ensembl_gene_ID_chimp", "Ensembl_gene_ID") # Rename columns for merging later

# Human vs chimp
hc <- read.csv("clnhumchimp.csv", header=TRUE, sep=",") # 15947
# Find human ID repeats
geneslost <- hc[(duplicated(hc$Ensembl_gene_ID) == TRUE),]
# Consolidate into just one list with gene names to make some sense of what is lost (not really helpful, actually)
geneslostx <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id", values = geneslost$Ensembl_gene_ID, mart = marthuman) #131
write.csv(geneslostx, "Lost_mean_HC.csv", row.names=FALSE)
# Find repeated human IDs and consolidate dN/dS values into just one mean value
sghc <- ddply(hc, "Ensembl_gene_ID", numcolwise(mean)) # 15755
# Check for repeated IDs was done again. No repeats.

# Human vs orangutan
ho <- read.csv("clnhumorang.csv", header=TRUE, sep=",") # 17099
geneslost <- ho[(duplicated(ho$Ensembl_gene_ID) == TRUE),]
geneslostx <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id", values = geneslost$Ensembl_gene_ID, mart = marthuman) # 300
write.csv(geneslostx, "Lost_mean_HO.csv", row.names=FALSE)
sgho <- ddply(ho, "Ensembl_gene_ID", numcolwise(mean)) # 16632
# Check for repeated IDs was done again. No repeats.

# Human vs macaca
hm <- read.csv("clnhummac.csv", header=TRUE, sep=",") # 17794
geneslost <- hm[(duplicated(hm$Ensembl_gene_ID) == TRUE),]
geneslostx <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id", values = geneslost$Ensembl_gene_ID, mart = marthuman) #981
write.csv(geneslostx, "Lost_mean_HM.csv", row.names=FALSE)
sghm <- ddply(hm, "Ensembl_gene_ID", numcolwise(mean)) # 16813
# Check for repeated IDs was done again. No repeats.

# Chimp vs orangutan
co <- read.csv("clnchimporang.csv", header=TRUE, sep=",") # 15995
# Find chimp ID dups
geneslost <- co[(duplicated(co$Ensembl_gene_ID_chimp) == TRUE),]
geneslostx <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id", values = geneslost$Ensembl_gene_ID_chimp, mart = martchimp) #240
write.csv(geneslostx, "Lost_mean_CO.csv", row.names=FALSE)
sgco <- ddply(co, "Ensembl_gene_ID_chimp", numcolwise(mean)) # 15677
# Now match back to human IDs
sgcomatch <- merge(IDhumanchimp, sgco, by="Ensembl_gene_ID_chimp") # 15677
# Check for repeated IDs
duphum <- sgcomatch[(duplicated(sgcomatch$Ensembl_gene_ID) == TRUE), ]
# There are repeats. These are "1 human to many chimp" orthologue cases. Take means again.
geneslosty <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id", values = duphum$Ensembl_gene_ID, mart = marthuman) # 80
write.csv(geneslosty, "Lost_one2many_CO.csv", row.names=FALSE)
# Check for repeated IDs was done again. No repeats.
# Drop chimp ID column
sgcomatch$Ensembl_gene_ID_chimp <- NULL
# Take means
sgco <- ddply(sgcomatch, "Ensembl_gene_ID", numcolwise(mean)) # 15588 left

# Chimp vs macaca
cm <- read.csv("clnchimpmac.csv", header=TRUE, sep=",") # 16636
geneslost <- cm[(duplicated(cm$Ensembl_gene_ID_chimp) == TRUE),]
geneslostx <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id", values = geneslost$Ensembl_gene_ID_chimp, mart = martchimp) # 558
write.csv(geneslostx, "Lost_mean_CM.csv", row.names=FALSE)
sgcm <- ddply(cm, "Ensembl_gene_ID_chimp", numcolwise(mean)) # 15812
sgcmmatch <- merge(IDhumanchimp, sgcm, by="Ensembl_gene_ID_chimp") # 15812
duphum <- sgcmmatch[(duplicated(sgcmmatch$Ensembl_gene_ID) == TRUE), ]
# There are repeats. These are "1 human to many chimp" orthologue cases. Take means again.
geneslosty <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id", values = duphum$Ensembl_gene_ID, mart = marthuman) # 77
write.csv(geneslosty, "Lost_one2many_CM.csv", row.names=FALSE)
# Check for repeated IDs was done again. No repeats.
sgcmmatch$Ensembl_gene_ID_chimp <- NULL
sgcm <- ddply(sgcmmatch, "Ensembl_gene_ID", numcolwise(mean)) # 15728

# Read in Human and Orangutan ID list
IDhumanorang <- read.csv("Gene_IDs_Human_Orang.csv", header = TRUE, sep = ",")
IDhumanorang$pabelii_homolog_orthology_type <- NULL
colnames(IDhumanorang) <- c("Ensembl_gene_ID_orang", "Ensembl_gene_ID")

# Orang vs macaca
om <- read.csv("clnorangmac.csv", header=TRUE, sep=",") # 16652
geneslost <- om[(duplicated(om$Ensembl_gene_ID_orang) == TRUE),]
geneslostx <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id", values = geneslost$Ensembl_gene_ID_orang, mart = martorang) #533
write.csv(geneslostx, "Lost_mean_OM.csv", row.names=FALSE)
sgom <- ddply(om, "Ensembl_gene_ID_orang", numcolwise(mean))# 15847
sgommatch <- merge(IDhumanorang, sgom, by="Ensembl_gene_ID_orang") # 15847
duphum <- sgommatch[(duplicated(sgommatch$Ensembl_gene_ID) == TRUE), ]
# There are repeats. These are "1 human to many orangutan" orthologue cases. Take means again.
geneslosty <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id", values = duphum$Ensembl_gene_ID, mart = marthuman) # 179
write.csv(geneslosty, "Lost_one2many_OM.csv", row.names=FALSE)
# Check for repeated IDs was done again. No repeats.
sgcmmatch$Ensembl_gene_ID_orang <- NULL
sgom <- ddply(sgommatch, "Ensembl_gene_ID", numcolwise(mean)) # 15648

############################################### goseq analysis ###############################################

# Run gene list through goseq() to determine enrichment or underrepresentation
# of various GO categories with dN/dS >1 as "DE"

# Subset genes with dN/dS >1
# Human vs chimp
dndshighHC <- subset(sghc, sghc$humchimpdndsr >1) #1064
# Human vs orangutan
dndshighHO <- subset(sgho, sgho$humorangdndsr >1) #579
# Human vs macaca
dndshighHM <- subset(sghm, sghm$hummacdndsr >1) #321
# Chimp vs orangutan
dndshighCO <- subset(sgco, sgco$chimporangdndsr >1) #472
# Chimp vs macaca
dndshighCM <- subset(sgcm, sgcm$chimpmacdndsr >1) #279
# Orang vs macaca
dndshighOM <- subset(sgom, sgom$orangmacdndsr >1) #264

# Human vs chimp
# Constructing an input table
# "Differentially expressed" genes are genes with dN/dS >1
x <- data.frame(dndshighHC$Ensembl_gene_ID, DE=1) # Assign "1" to the "DE" genes
colnames(x) <- c("Ensembl_gene_ID", "DE")
labeled <- merge(sghc, x, all.x = TRUE, by="Ensembl_gene_ID")
# Convert NA's to 0
labeled[is.na(labeled)] <- 0
# Keep only necessary columns
labeled <- labeled[c("Ensembl_gene_ID", "DE")]
# Convert to goseq inputtable format
genelist <- as.integer(labeled$DE)
names(genelist) <- labeled$Ensembl_gene_ID
# For getlength to work:
genelistx <- names(genelist)
genelength <- getlength(genelistx, "hg19", "ensGene")
# Use source("http://bioconductor.org/biocLite.R") then biocLite("org.Hs.eg.db")
# to manually download the database if getlength() doesn't work
# Create pwf
pwf <- nullp(genelist, "hg19", "ensGene", bias.data=genelength, plot=FALSE)
genetocat <- getgo(genelistx, "hg19", "ensGene")
rownames(pwf) <- names(genelist) # Somehow this makes it work
# Running goseq()
output <- goseq(pwf, "hg19", genelist, gene2cat = genetocat, 
                test.cats=c("GO:BP","GO:MF","GO:CC"),
                method = "Wallenius",
                use_genes_without_cat=FALSE) # Excluded 806 genes without categories
# Multiple testing correction; FDR correction
enriched <- output$category[p.adjust(output$over_represented_pvalue, method="BH") < 0.05] # 21 categories
# Keep only truly enriched ones
enrichedoutput <- output$category %in% enriched
outpute <- output[enrichedoutput,]
outpute$under_represented_pvalue <- NULL
colnames(outpute) <- c("GO category","p-value","No. of DE genes","No. of genes in cat", "GO term","Ontology")
write.csv(outpute, "enrichedHC.csv", row.names=FALSE)
# Reversemapping (getting genes out of enriched categories)
cattogenes <- goseq:::reversemapping(genetocat) # Convert categories to genes
# Get genes for every outpute$category
relevantcats <- cattogenes[outpute$`GO category`] # This cuts down the CATEGORIES to only those in outpute
# Cut down GENES to only those present in the DE group
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighHC$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
# enddf is now a long dataframe of each gene with their associated category
write.csv(enddf, "enriched_HC_genecat.csv", row.names=FALSE)

# Human vs orang
x <- data.frame(dndshighHO$Ensembl_gene_ID, DE=1)
colnames(x) <- c("Ensembl_gene_ID", "DE")
labeled <- merge(sgho, x, all.x = TRUE, by="Ensembl_gene_ID")
labeled[is.na(labeled)] <- 0
labeled <- labeled[c("Ensembl_gene_ID", "DE")]
genelist <- as.integer(labeled$DE)
names(genelist) <- labeled$Ensembl_gene_ID
genelistx <- names(genelist)
genelength <- getlength(genelistx, "hg19", "ensGene")
pwf <- nullp(genelist, "hg19", "ensGene", bias.data=genelength, plot=FALSE)
genetocat <- getgo(genelistx, "hg19", "ensGene")
rownames(pwf) <- names(genelist)
output <- goseq(pwf, "hg19", genelist, gene2cat = genetocat, 
                test.cats=c("GO:BP","GO:MF","GO:CC"),
                method = "Wallenius", use_genes_without_cat=FALSE) # Excluded 839 genes without categories
enriched <- output$category[p.adjust(output$over_represented_pvalue, method="BH") < 0.05] # 7 categories
enrichedoutput <- output$category %in% enriched
outpute <- output[enrichedoutput,]
outpute$under_represented_pvalue <- NULL
colnames(outpute) <- c("GO category","p-value","No. of DE genes","No. of genes in cat", "GO term","Ontology")
write.csv(outpute, "enrichedHO.csv", row.names=FALSE)
cattogenes <- goseq:::reversemapping(genetocat)
relevantcats <- cattogenes[outpute$`GO category`]
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighHO$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
write.csv(enddf, "enriched_HO_genecat.csv", row.names=FALSE)

# Human vs macaca
x <- data.frame(dndshighHM$Ensembl_gene_ID, DE=1)
colnames(x) <- c("Ensembl_gene_ID", "DE")
labeled <- merge(sghm, x, all.x = TRUE, by="Ensembl_gene_ID")
labeled[is.na(labeled)] <- 0
labeled <- labeled[c("Ensembl_gene_ID", "DE")]
genelist <- as.integer(labeled$DE)
names(genelist) <- labeled$Ensembl_gene_ID
genelistx <- names(genelist)
genelength <- getlength(genelistx, "hg19", "ensGene")
pwf <- nullp(genelist, "hg19", "ensGene", bias.data=genelength, plot=FALSE)
genetocat <- getgo(genelistx, "hg19", "ensGene")
rownames(pwf) <- names(genelist)
output <- goseq(pwf, "hg19", genelist, gene2cat = genetocat, 
                test.cats=c("GO:BP","GO:MF","GO:CC"),
                method = "Wallenius", use_genes_without_cat=FALSE) # Excluded 799 genes without categories
enriched <- output$category[p.adjust(output$over_represented_pvalue, method="BH") < 0.05] # 17 categories
enrichedoutput <- output$category %in% enriched
outpute <- output[enrichedoutput,]
outpute$under_represented_pvalue <- NULL
colnames(outpute) <- c("GO category","p-value","No. of DE genes","No. of genes in cat", "GO term","Ontology")
write.csv(outpute, "enrichedHM.csv", row.names=FALSE)
cattogenes <- goseq:::reversemapping(genetocat)
relevantcats <- cattogenes[outpute$`GO category`] 
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighHM$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
write.csv(enddf, "enriched_HM_genecat.csv", row.names=FALSE)

# Chimp vs orang
x <- data.frame(dndshighCO$Ensembl_gene_ID, DE=1)
colnames(x) <- c("Ensembl_gene_ID", "DE")
labeled <- merge(sgco, x, all.x = TRUE, by="Ensembl_gene_ID")
labeled[is.na(labeled)] <- 0
labeled <- labeled[c("Ensembl_gene_ID", "DE")]
genelist <- as.integer(labeled$DE)
names(genelist) <- labeled$Ensembl_gene_ID
genelistx <- names(genelist)
genelength <- getlength(genelistx, "hg19", "ensGene")
pwf <- nullp(genelist, "hg19", "ensGene", bias.data=genelength, plot=FALSE)
genetocat <- getgo(genelistx, "hg19", "ensGene")
rownames(pwf) <- names(genelist)
output <- goseq(pwf, "hg19", genelist, gene2cat = genetocat, 
                test.cats=c("GO:BP","GO:MF","GO:CC"),
                method = "Wallenius", use_genes_without_cat=FALSE) # Excluded 686 genes without categories
enriched <- output$category[p.adjust(output$over_represented_pvalue, method="BH") < 0.05] # 6 categories
enrichedoutput <- output$category %in% enriched
outpute <- output[enrichedoutput,]
outpute$under_represented_pvalue <- NULL
colnames(outpute) <- c("GO category","p-value","No. of DE genes","No. of genes in cat", "GO term","Ontology")
write.csv(outpute, "enrichedCO.csv", row.names=FALSE)
cattogenes <- goseq:::reversemapping(genetocat)
relevantcats <- cattogenes[outpute$`GO category`]
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighCO$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
write.csv(enddf, "enriched_CO_genecat.csv", row.names=FALSE)

# Chimp vs macaca
x <- data.frame(dndshighCM$Ensembl_gene_ID, DE=1)
colnames(x) <- c("Ensembl_gene_ID", "DE")
labeled <- merge(sgcm, x, all.x = TRUE, by="Ensembl_gene_ID")
labeled[is.na(labeled)] <- 0
labeled <- labeled[c("Ensembl_gene_ID", "DE")]
genelist <- as.integer(labeled$DE)
names(genelist) <- labeled$Ensembl_gene_ID
genelistx <- names(genelist)
genelength <- getlength(genelistx, "hg19", "ensGene")
pwf <- nullp(genelist, "hg19", "ensGene", bias.data=genelength, plot=FALSE)
genetocat <- getgo(genelistx, "hg19", "ensGene")
rownames(pwf) <- names(genelist)
output <- goseq(pwf, "hg19", genelist, gene2cat = genetocat, 
                test.cats=c("GO:BP","GO:MF","GO:CC"),
                method = "Wallenius", use_genes_without_cat=FALSE) # Excluded 657 genes without categories
enriched <- output$category[p.adjust(output$over_represented_pvalue, method="BH") < 0.05] # 20 categories
enrichedoutput <- output$category %in% enriched
outpute <- output[enrichedoutput,]
outpute$under_represented_pvalue <- NULL
colnames(outpute) <- c("GO category","p-value","No. of DE genes","No. of genes in cat", "GO term","Ontology")
write.csv(outpute, "enrichedCM.csv", row.names=FALSE)
cattogenes <- goseq:::reversemapping(genetocat)
relevantcats <- cattogenes[outpute$`GO category`]
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighCM$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
write.csv(enddf, "enriched_CM_genecat.csv", row.names=FALSE)

# Orang vs macaca
x <- data.frame(dndshighOM$Ensembl_gene_ID, DE=1)
colnames(x) <- c("Ensembl_gene_ID", "DE")
labeled <- merge(sgom, x, all.x = TRUE, by="Ensembl_gene_ID")
labeled[is.na(labeled)] <- 0
labeled <- labeled[c("Ensembl_gene_ID", "DE")]
genelist <- as.integer(labeled$DE)
names(genelist) <- labeled$Ensembl_gene_ID
genelistx <- names(genelist)
genelength <- getlength(genelistx, "hg19", "ensGene")
pwf <- nullp(genelist, "hg19", "ensGene", bias.data=genelength, plot=FALSE)
genetocat <- getgo(genelistx, "hg19", "ensGene")
rownames(pwf) <- names(genelist)
output <- goseq(pwf, "hg19", genelist, gene2cat = genetocat, 
                test.cats=c("GO:BP","GO:MF","GO:CC"),
                method = "Wallenius", use_genes_without_cat=FALSE) # Excluded 625 genes without categories
enriched <- output$category[p.adjust(output$over_represented_pvalue, method="BH") < 0.05] # 13 categories
enrichedoutput <- output$category %in% enriched
outpute <- output[enrichedoutput,]
outpute$under_represented_pvalue <- NULL
colnames(outpute) <- c("GO category","p-value","No. of DE genes","No. of genes in cat", "GO term","Ontology")
write.csv(outpute, "enrichedOM.csv", row.names=FALSE)
cattogenes <- goseq:::reversemapping(genetocat)
relevantcats <- cattogenes[outpute$`GO category`]
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighOM$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
write.csv(enddf, "enriched_OM_genecat.csv", row.names=FALSE)


# Create long GO category list
eHC <- read.csv("enrichedHC.csv",header=TRUE,sep=",",stringsAsFactors = FALSE)
eHO <- read.csv("enrichedHO.csv",header=TRUE,sep=",",stringsAsFactors = FALSE)
eHM <- read.csv("enrichedHM.csv",header=TRUE,sep=",",stringsAsFactors = FALSE)
eCO <- read.csv("enrichedCO.csv",header=TRUE,sep=",",stringsAsFactors = FALSE)
eCM <- read.csv("enrichedCM.csv",header=TRUE,sep=",",stringsAsFactors = FALSE)
eOM <- read.csv("enrichedOM.csv",header=TRUE,sep=",",stringsAsFactors = FALSE)
altogether <- rbind(eHC,eHO,eHM,eCO,eCM,eOM)

# Finding GO categories unique to each comparison
# Human vs chimp
# Starts with 21
eHCx <- eHC[!eHC$GO.category %in% eHO$GO.category,]
eHCx <- eHCx[!eHCx$GO.category %in% eHM$GO.category,]
eHCx <- eHCx[!eHCx$GO.category %in% eCO$GO.category,]
eHCx <- eHCx[!eHCx$GO.category %in% eCM$GO.category,]
eHCx <- eHCx[!eHCx$GO.category %in% eOM$GO.category,] # 15 categories
write.csv(eHCx, "uniqueGO_HC.csv",row.names=FALSE)

# Human vs orangutan
# Starts with 7
eHOx <- eHO[!eHO$GO.category %in% eHC$GO.category,]
eHOx <- eHOx[!eHOx$GO.category %in% eHM$GO.category,]
eHOx <- eHOx[!eHOx$GO.category %in% eCO$GO.category,]
eHOx <- eHOx[!eHOx$GO.category %in% eCM$GO.category,]
eHOx <- eHOx[!eHOx$GO.category %in% eOM$GO.category,] # 0 categories

# Human vs macaca
# Starts with 17
eHMx <- eHM[!eHM$GO.category %in% eHC$GO.category,]
eHMx <- eHMx[!eHMx$GO.category %in% eHO$GO.category,]
eHMx <- eHMx[!eHMx$GO.category %in% eCO$GO.category,]
eHMx <- eHMx[!eHMx$GO.category %in% eCM$GO.category,]
eHMx <- eHMx[!eHMx$GO.category %in% eOM$GO.category,] # 1 category
write.csv(eHMx, "uniqueGO_HM.csv",row.names=FALSE)

# Chimp vs orang
# Starts with 6
eCOx <- eCO[!eCO$GO.category %in% eHC$GO.category,]
eCOx <- eCOx[!eCOx$GO.category %in% eHO$GO.category,]
eCOx <- eCOx[!eCOx$GO.category %in% eHM$GO.category,]
eCOx <- eCOx[!eCOx$GO.category %in% eCM$GO.category,]
eCOx <- eCOx[!eCOx$GO.category %in% eOM$GO.category,] # 0 categories

# Chimp vs macaca
# Starts with 20
eCMx <- eCM[!eCM$GO.category %in% eHC$GO.category,]
eCMx <- eCMx[!eCMx$GO.category %in% eHO$GO.category,]
eCMx <- eCMx[!eCMx$GO.category %in% eHM$GO.category,]
eCMx <- eCMx[!eCMx$GO.category %in% eCO$GO.category,]
eCMx <- eCMx[!eCMx$GO.category %in% eOM$GO.category,] # 4 categories
write.csv(eCMx, "uniqueGO_CM.csv",row.names=FALSE)

# Orang vs macaca
# Starts with 13
eOMx <- eOM[!eOM$GO.category %in% eHC$GO.category,]
eOMx <- eOMx[!eOMx$GO.category %in% eHO$GO.category,]
eOMx <- eOMx[!eOMx$GO.category %in% eHM$GO.category,]
eOMx <- eOMx[!eOMx$GO.category %in% eCO$GO.category,]
eOMx <- eOMx[!eOMx$GO.category %in% eCM$GO.category,] # 0 categories

# Finding common GO categories across comparisons
allcatsandterms <- altogether[c("GO.category","GO.term")]
countrepeats <- ddply(allcatsandterms,c("GO.category","GO.term"),nrow)
commonGOcats <- subset(countrepeats, countrepeats$V1=="6") # 2 categories
write.csv(commonGOcats,"common_GOcategories.csv",row.names=FALSE)
allhuman <- rbind(eHC,eHM,eHO)
allcatsandtermsH <- allhuman[c("GO.category","GO.term")]
countrepeatsH <- ddply(allhuman,c("GO.category","GO.term"),nrow)
commonGOcatsH <- subset(countrepeatsH, countrepeatsH$V1=="3") # 2 categories
write.csv(commonGOcatsH,"common_GOcategoriesH.csv",row.names=FALSE)
allchimp <- rbind(eHC,eCM,eCO)
allcatsandtermsC <- allchimp[c("GO.category","GO.term")]
countrepeatsC <- ddply(allchimp,c("GO.category","GO.term"),nrow)
commonGOcatsC <- subset(countrepeatsC, countrepeatsC$V1=="3") # 2 categories
write.csv(commonGOcatsC,"common_GOcategoriesC.csv",row.names=FALSE)
allorang <- rbind(eHO,eCO,eOM)
allcatsandtermsO <- allorang[c("GO.category","GO.term")]
countrepeatsO <- ddply(allorang,c("GO.category","GO.term"),nrow)
commonGOcatsO <- subset(countrepeatsO, countrepeatsO$V1=="3") # 3 categories
write.csv(commonGOcatsO,"common_GOcategoriesO.csv",row.names=FALSE)
allmac <- rbind(eHM,eCM,eOM)
allcatsandtermsM <- allorang[c("GO.category","GO.term")]
countrepeatsM <- ddply(allmac,c("GO.category","GO.term"),nrow)
commonGOcatsM <- subset(countrepeatsM, countrepeatsM$V1=="3") # 13 categories
write.csv(commonGOcatsM,"common_GOcategoriesM.csv",row.names=FALSE)

# Finding unique GO categories across comparisons
uniqueGOcats <- subset(countrepeats, countrepeats$V1=="1") # 20 unique categories
write.csv(uniqueGOcats,"unique_GOcategories.csv",row.names=FALSE)
uniqueGOcatsH <- subset(countrepeatsH, countrepeatsH$V1=="1") # 29 unique categories
write.csv(uniqueGOcatsH,"unique_GOcategoriesH.csv",row.names=FALSE)
uniqueGOcatsC <- subset(countrepeatsC, countrepeatsC$V1=="1") # 33 categories
write.csv(uniqueGOcatsC,"unique_GOcategoriesC.csv",row.names=FALSE)
uniqueGOcatsO <- subset(countrepeatsO, countrepeatsO$V1=="1") # 11 categories
write.csv(uniqueGOcatsO,"unique_GOcategoriesO.csv",row.names=FALSE)
uniqueGOcatsM <- subset(countrepeatsM, countrepeatsM$V1=="1") # 5 categories
write.csv(uniqueGOcatsM,"unique_GOcategoriesM.csv",row.names=FALSE)

#Finding GO cats unique to a certain species comparison
uhuman <- uniqueGOcatsH[!uniqueGOcatsH$GO.category %in% uniqueGOcatsC$GO.category,]
uhuman <- uhuman[!uhuman$GO.category %in% uniqueGOcatsO$GO.category,]
uhuman <- uhuman[!uhuman$GO.category %in% uniqueGOcatsM$GO.category,] # 3 categories (all taste-related)

uchimp <- uniqueGOcatsC[!uniqueGOcatsC$GO.category %in% uniqueGOcatsH$GO.category,]
uchimp <- uchimp[!uchimp$GO.category %in% uniqueGOcatsO$GO.category,]
uchimp <- uchimp[!uchimp$GO.category %in% uniqueGOcatsM$GO.category,] # No categories

uorang <- uniqueGOcatsO[!uniqueGOcatsO$GO.category %in% uniqueGOcatsH$GO.category,]
uorang <- uorang[!uorang$GO.category %in% uniqueGOcatsC$GO.category,]
uorang <- uorang[!uorang$GO.category %in% uniqueGOcatsM$GO.category,] # No categories

umac <- uniqueGOcatsM[!uniqueGOcatsM$GO.category %in% uniqueGOcatsH$GO.category,]
umac <- umac[!umac$GO.category %in% uniqueGOcatsC$GO.category,]
umac <- umac[!umac$GO.category %in% uniqueGOcatsO$GO.category,] # No categories

######################################### Identify genes in category #########################################

# Condense gene lists to only non-redundant ones and drop GO categories
# Also append corresponding pLI scores

# Human vs chimp
genecatlist <- read.csv("enriched_HC_genecat.csv", header=TRUE, sep=",")
# Remove repeated genes and drop category names
genecatlist <- genecatlist[duplicated(genecatlist$gene)==FALSE,]
genecatlist <- genecatlist["gene"] # 99
# What genes are these?
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                     filters = "ensembl_gene_id", values = genecatlist$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
# Append pLI scores
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
# Append dN/dS values
all.dndsr <- dndshighHC[c("Ensembl_gene_ID","humchimpdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
write.csv(knowgene, "genelist_HC.csv", row.names=FALSE)

# Human vs orang
genecatlist <- read.csv("enriched_HO_genecat.csv", header=TRUE, sep=",")
genecatlist <- genecatlist[duplicated(genecatlist$gene)==FALSE,]
genecatlist <- genecatlist["gene"] # 64
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = genecatlist$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighHO[c("Ensembl_gene_ID","humorangdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
write.csv(knowgene, "genelist_HO.csv", row.names=FALSE)

# Human vs Macaca
genecatlist <- read.csv("enriched_HM_genecat.csv", header=TRUE, sep=",")
genecatlist <- genecatlist[duplicated(genecatlist$gene)==FALSE,]
genecatlist <- genecatlist["gene"] # 101
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = genecatlist$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighHM[c("Ensembl_gene_ID","hummacdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
write.csv(knowgene, "genelist_HM.csv", row.names=FALSE)

# Chimp vs orang
genecatlist <- read.csv("enriched_CO_genecat.csv", header=TRUE, sep=",")
genecatlist <- genecatlist[duplicated(genecatlist$gene)==FALSE,]
genecatlist <- genecatlist["gene"] # 12
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = genecatlist$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighCO[c("Ensembl_gene_ID","chimporangdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
write.csv(knowgene, "genelist_CO.csv", row.names=FALSE)

# Chimp vs macaca
genecatlist <- read.csv("enriched_CM_genecat.csv", header=TRUE, sep=",")
genecatlist <- genecatlist[duplicated(genecatlist$gene)==FALSE,]
genecatlist <- genecatlist["gene"] # 93
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = genecatlist$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighCM[c("Ensembl_gene_ID","chimpmacdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
write.csv(knowgene, "genelist_CM.csv", row.names=FALSE)

# Orang vs macaca
genecatlist <- read.csv("enriched_OM_genecat.csv", header=TRUE, sep=",")
genecatlist <- genecatlist[duplicated(genecatlist$gene)==FALSE,]
genecatlist <- genecatlist["gene"] # 69
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = genecatlist$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighOM[c("Ensembl_gene_ID","orangmacdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
write.csv(knowgene, "genelist_OM.csv", row.names=FALSE)

################################## Finding genes unique to each comparison ###################################

genelisthc <- read.csv("genelist_HC.csv",header=TRUE,sep=",")
genelistho <- read.csv("genelist_HO.csv",header=TRUE,sep=",")
genelisthm <- read.csv("genelist_HM.csv",header=TRUE,sep=",")
genelistco <- read.csv("genelist_CO.csv",header=TRUE,sep=",")
genelistcm <- read.csv("genelist_CM.csv",header=TRUE,sep=",")
genelistom <- read.csv("genelist_OM.csv",header=TRUE,sep=",")

# Human vs chimp
# Starts with 99
genelisthcx <- genelisthc[!genelisthc$Ensembl_gene_ID %in% genelisthm$Ensembl_gene_ID,]
genelisthcx <- genelisthcx[!genelisthcx$Ensembl_gene_ID %in% genelistho$Ensembl_gene_ID,]
genelisthcx <- genelisthcx[!genelisthcx$Ensembl_gene_ID %in% genelistco$Ensembl_gene_ID,]
genelisthcx <- genelisthcx[!genelisthcx$Ensembl_gene_ID %in% genelistcm$Ensembl_gene_ID,]
genelisthcx <- genelisthcx[!genelisthcx$Ensembl_gene_ID %in% genelistom$Ensembl_gene_ID,] # 79 genes
write.csv(genelisthcx,"uniquegenesHC.csv",row.names=FALSE)

# Human vs orang
# Starts with 17
genelisthox <- genelistho[!genelistho$Ensembl_gene_ID %in% genelisthc$Ensembl_gene_ID,]
genelisthox <- genelisthox[!genelisthox$Ensembl_gene_ID %in% genelisthm$Ensembl_gene_ID,]
genelisthox <- genelisthox[!genelisthox$Ensembl_gene_ID %in% genelistco$Ensembl_gene_ID,]
genelisthox <- genelisthox[!genelisthox$Ensembl_gene_ID %in% genelistcm$Ensembl_gene_ID,]
genelisthox <- genelisthox[!genelisthox$Ensembl_gene_ID %in% genelistom$Ensembl_gene_ID,] # 0 genes

# Human vs macaca
# Starts with 101
genelisthmx <- genelisthm[!genelisthm$Ensembl_gene_ID %in% genelisthc$Ensembl_gene_ID,]
genelisthmx <- genelisthmx[!genelisthmx$Ensembl_gene_ID %in% genelistho$Ensembl_gene_ID,]
genelisthmx <- genelisthmx[!genelisthmx$Ensembl_gene_ID %in% genelistco$Ensembl_gene_ID,]
genelisthmx <- genelisthmx[!genelisthmx$Ensembl_gene_ID %in% genelistcm$Ensembl_gene_ID,]
genelisthmx <- genelisthmx[!genelisthmx$Ensembl_gene_ID %in% genelistom$Ensembl_gene_ID,] # 24 genes
write.csv(genelisthmx,"uniquegenesHM.csv",row.names=FALSE)

# Chimp vs orang
# Starts with 12
genelistcox <- genelistco[!genelistco$Ensembl_gene_ID %in% genelisthc$Ensembl_gene_ID,]
genelistcox <- genelistcox[!genelistcox$Ensembl_gene_ID %in% genelistho$Ensembl_gene_ID,]
genelistcox <- genelistcox[!genelistcox$Ensembl_gene_ID %in% genelisthm$Ensembl_gene_ID,]
genelistcox <- genelistcox[!genelistcox$Ensembl_gene_ID %in% genelistcm$Ensembl_gene_ID,]
genelistcox <- genelistcox[!genelistcox$Ensembl_gene_ID %in% genelistom$Ensembl_gene_ID,] # 0 genes

# Chimp vs macaca
# Starts with 93
genelistcmx <- genelistcm[!genelistcm$Ensembl_gene_ID %in% genelisthc$Ensembl_gene_ID,]
genelistcmx <- genelistcmx[!genelistcmx$Ensembl_gene_ID %in% genelistho$Ensembl_gene_ID,]
genelistcmx <- genelistcmx[!genelistcmx$Ensembl_gene_ID %in% genelisthm$Ensembl_gene_ID,]
genelistcmx <- genelistcmx[!genelistcmx$Ensembl_gene_ID %in% genelistco$Ensembl_gene_ID,]
genelistcmx <- genelistcmx[!genelistcmx$Ensembl_gene_ID %in% genelistom$Ensembl_gene_ID,] # 16 genes
write.csv(genelistcmx,"uniquegenesCM.csv",row.names=FALSE)

# Orang vs macaca
# Starts with 69
genelistomx <- genelistom[!genelistom$Ensembl_gene_ID %in% genelisthc$Ensembl_gene_ID,]
genelistomx <- genelistomx[!genelistomx$Ensembl_gene_ID %in% genelistho$Ensembl_gene_ID,]
genelistomx <- genelistomx[!genelistomx$Ensembl_gene_ID %in% genelisthm$Ensembl_gene_ID,]
genelistomx <- genelistomx[!genelistomx$Ensembl_gene_ID %in% genelistco$Ensembl_gene_ID,]
genelistomx <- genelistomx[!genelistomx$Ensembl_gene_ID %in% genelistcm$Ensembl_gene_ID,] # 20 genes
write.csv(genelistomx,"uniquegenesOM.csv",row.names=FALSE)

# Finding genes unique to comparisons with a certain species did not yield any genes
genelistH <- rbind(genelisthc, genelisthm, genelistho)
genelistC <- rbind(genelisthc, genelistco, genelistcm)
genelistO <- rbind(genelistho, genelistco, genelistom)
genelistM <- rbind(genelisthm, genelistcm, genelistom)

# Human
genelistHx <- genelistH[!genelistH$Ensembl_gene_ID %in% genelistC$Ensembl_gene_ID,]
genelistHx <- genelistHx[!genelistHx$Ensembl_gene_ID %in% genelistO$Ensembl_gene_ID,]
genelistHx <- genelistHx[!genelistHx$Ensembl_gene_ID %in% genelistM$Ensembl_gene_ID,] # 0 genes

# Chimp
genelistCx <- genelistC[!genelistC$Ensembl_gene_ID %in% genelistH$Ensembl_gene_ID,]
genelistCx <- genelistCx[!genelistCx$Ensembl_gene_ID %in% genelistO$Ensembl_gene_ID,]
genelistCx <- genelistCx[!genelistCx$Ensembl_gene_ID %in% genelistM$Ensembl_gene_ID,] # 0 genes

# Orang
genelistOx <- genelistO[!genelistO$Ensembl_gene_ID %in% genelistH$Ensembl_gene_ID,]
genelistOx <- genelistOx[!genelistOx$Ensembl_gene_ID %in% genelistC$Ensembl_gene_ID,]
genelistOx <- genelistOx[!genelistOx$Ensembl_gene_ID %in% genelistM$Ensembl_gene_ID,] # 0 genes

# Macaca
genelistMx <- genelistM[!genelistM$Ensembl_gene_ID %in% genelistH$Ensembl_gene_ID,]
genelistMx <- genelistMx[!genelistMx$Ensembl_gene_ID %in% genelistC$Ensembl_gene_ID,]
genelistMx <- genelistMx[!genelistMx$Ensembl_gene_ID %in% genelistO$Ensembl_gene_ID,] # 0 genes

######################################### Getting sequences of genes #########################################

# Use gene list from enriched GO categories unique to pairwise comparisons
# Retrieve protein sequences for usein PAML
# Orangutan comparions yielded no genes except for orang vs macaca

HCgenes <- read.csv("uniquegenesHC.csv",header=TRUE,sep=",")
HMgenes <- read.csv("uniquegenesHM.csv",header=TRUE,sep=",")
CMgenes <- read.csv("uniquegenesCM.csv",header=TRUE,sep=",")
OMgenes <- read.csv("uniquegenesOM.csv",header=TRUE,sep=",")

seqHC <- getSequence(id=HCgenes$Gene_symbol, type="hgnc_symbol", seqType="cdna", mart=marthuman)

##############################################################################################################