# The giant script consolidating codes from start to finish without the extra codes that were written for testing.
# Redundant comments are also removed.

################################################### START ####################################################
# Processing all gene data in the ExAC dataset for use in downstream scripts

setwd("~/ExAC_dataset_project/csv_records")
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
exacgenes <- exacgenes[c("transcript", "gene", "pLI")]

# Get ensembl gene IDs
geneid <- getBM(attributes = c("ensembl_gene_id","ensembl_transcript_id"), filters = "ensembl_transcript_id", values = exacgenes$transcript, mart = marthuman, uniqueRows = FALSE)
# Append to original data frame and rename columns
colnames(geneid) <- c("Ensembl_gene_ID", "Ensembl_transcript_ID")
colnames(exacgenes) <- c("Ensembl_transcript_ID", "Gene_symbol", "pLI")
exacgenesx <- merge(geneid, exacgenes, by="Ensembl_transcript_ID")
write.csv(exacgenesx, "exacgenes.csv", row.names=FALSE)

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

# No orangutan ribosomal proteins present.

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
write.csv(orangduprec, "Lost_many2one_orang.csv", row.names=FALSE) # 703
oranggeneidy[(oranggeneidy$pabelii_homolog_ensembl_gene %in% orangdup$pabelii_homolog_ensembl_gene), ] <- NA
oranggeneidy <- oranggeneidy[!is.na(oranggeneidy$pabelii_homolog_ensembl_gene), ]

# Orangutan vs macaca
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
write.csv(geneslost, "Lost_dndsNA_HO.csv", row.names=FALSE)
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

write.csv(sghc,"sghc.csv",row.names=F)
write.csv(sgho,"sgho.csv",row.names=F)
write.csv(sghm,"sghm.csv",row.names=F)
write.csv(sgco,"sgco.csv",row.names=F)
write.csv(sgcm,"sgcm.csv",row.names=F)
write.csv(sgom,"sgom.csv",row.names=F)

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

write.csv(dndshighHC, "dndshighHC.csv", row.names=FALSE)
write.csv(dndshighHO, "dndshighHO.csv", row.names=FALSE)
write.csv(dndshighHM, "dndshighHM.csv", row.names=FALSE)
write.csv(dndshighCO, "dndshighCO.csv", row.names=FALSE)
write.csv(dndshighCM, "dndshighCM.csv", row.names=FALSE)
write.csv(dndshighOM, "dndshighOM.csv", row.names=FALSE)

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
genetocatHC <- genetocat # For future use
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
colnames(enddf) <- c("GO category", "Ensembl_gene_ID") # enddf is now a long dataframe of each gene with their associated category
# Append other useful information
enddfHC <- merge(enddf, outpute, by="GO category")
# Find out what genes are these
whatgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = enddfHC$Ensembl_gene_ID, mart = marthuman) # 176 genes
colnames(whatgene) <- c("Ensembl_gene_ID","Gene_symbol","Description")
enddfHC <- merge(whatgene, enddfHC, by="Ensembl_gene_ID")
enddfHC  <- enddfHC[,c("GO category","GO term","Ensembl_gene_ID","Gene_symbol","Description","p-value","Ontology")]
write.csv(enddfHC, "enriched_HC_genecat.csv", row.names=FALSE)
whatisunique <- unique(enddfHC$Ensembl_gene_ID)
length(whatisunique) # 176 unique genes

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
genetocatHO <- genetocat # For future use
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
colnames(enddf) <- c("GO category", "Ensembl_gene_ID")
enddfHO <- merge(enddf, outpute, by="GO category")
whatgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = enddfHO$Ensembl_gene_ID, mart = marthuman) # 41 genes
colnames(whatgene) <- c("Ensembl_gene_ID","Gene_symbol","Description")
enddfHO <- merge(whatgene, enddfHO, by="Ensembl_gene_ID")
enddfHO  <- enddfHO[,c("GO category","GO term","Ensembl_gene_ID","Gene_symbol","Description","p-value","Ontology")]
write.csv(enddfHO, "enriched_HO_genecat.csv", row.names=FALSE)
whatisunique <- unique(enddfHO$Ensembl_gene_ID)
length(whatisunique) # 41 unique genes

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
genetocatHM <- genetocat # For future use
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
colnames(enddf) <- c("GO category", "Ensembl_gene_ID")
enddfHM <- merge(enddf, outpute, by="GO category")
whatgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = enddfHM$Ensembl_gene_ID, mart = marthuman) # 101 genes
colnames(whatgene) <- c("Ensembl_gene_ID","Gene_symbol","Description")
enddfHM <- merge(whatgene, enddfHM, by="Ensembl_gene_ID")
enddfHM  <- enddfHM[,c("GO category","GO term","Ensembl_gene_ID","Gene_symbol","Description","p-value","Ontology")]
write.csv(enddfHM, "enriched_HM_genecat.csv", row.names=FALSE)
whatisunique <- unique(enddfHM$Ensembl_gene_ID)
length(whatisunique) # 101 unique genes

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
genetocatCO <- genetocat # For future use
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
colnames(enddf) <- c("GO category", "Ensembl_gene_ID")
enddfCO <- merge(enddf, outpute, by="GO category")
whatgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = enddfCO$Ensembl_gene_ID, mart = marthuman) # 12 genes
colnames(whatgene) <- c("Ensembl_gene_ID","Gene_symbol","Description")
enddfCO <- merge(whatgene, enddfCO, by="Ensembl_gene_ID")
enddfCO  <- enddfCO[,c("GO category","GO term","Ensembl_gene_ID","Gene_symbol","Description","p-value","Ontology")]
write.csv(enddfCO, "enriched_CO_genecat.csv", row.names=FALSE)
whatisunique <- unique(enddfCO$Ensembl_gene_ID)
length(whatisunique) # 12 unique genes

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
genetocatCM <- genetocat # For future use
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
colnames(enddf) <- c("GO category", "Ensembl_gene_ID")
enddfCM <- merge(enddf, outpute, by="GO category")
whatgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = enddfCM$Ensembl_gene_ID, mart = marthuman) # 93 genes
colnames(whatgene) <- c("Ensembl_gene_ID","Gene_symbol","Description")
enddfCM <- merge(whatgene, enddfCM, by="Ensembl_gene_ID")
enddfCM  <- enddfCM[,c("GO category","GO term","Ensembl_gene_ID","Gene_symbol","Description","p-value","Ontology")]
write.csv(enddfCM, "enriched_CM_genecat.csv", row.names=FALSE)
whatisunique <- unique(enddfCM$Ensembl_gene_ID)
length(whatisunique) # 93 unique genes

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
genetocatOM <- genetocat # For future use
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
colnames(enddf) <- c("GO category", "Ensembl_gene_ID")
enddfOM <- merge(enddf, outpute, by="GO category")
whatgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = enddfOM$Ensembl_gene_ID, mart = marthuman) # 69 genes
colnames(whatgene) <- c("Ensembl_gene_ID","Gene_symbol","Description")
enddfOM <- merge(whatgene, enddfOM, by="Ensembl_gene_ID")
enddfOM  <- enddfOM[,c("GO category","GO term","Ensembl_gene_ID","Gene_symbol","Description","p-value","Ontology")]
write.csv(enddfOM, "enriched_OM_genecat.csv", row.names=FALSE)
whatisunique <- unique(enddfOM$Ensembl_gene_ID)
length(whatisunique) # 69 unique genes

########################################## Comparing GO categories ###########################################

# Create long GO category list
eHC <- read.csv("enrichedHC.csv",header=TRUE,sep=",",stringsAsFactors = FALSE)
eHO <- read.csv("enrichedHO.csv",header=TRUE,sep=",",stringsAsFactors = FALSE)
eHM <- read.csv("enrichedHM.csv",header=TRUE,sep=",",stringsAsFactors = FALSE)
eCO <- read.csv("enrichedCO.csv",header=TRUE,sep=",",stringsAsFactors = FALSE)
eCM <- read.csv("enrichedCM.csv",header=TRUE,sep=",",stringsAsFactors = FALSE)
eOM <- read.csv("enrichedOM.csv",header=TRUE,sep=",",stringsAsFactors = FALSE)
altogether <- rbind(eHC,eHO,eHM,eCO,eCM,eOM)

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

######################################### Identify genes in category #########################################

# Find genes within common enriched GO categories and manually remove redundant GO categories
# Genes are all assigned to only one GO category (no gene repeats)

# Human
# Human vs chimp
cattogenes <- goseq:::reversemapping(genetocatHC) # Convert categories to genes
Hcat <- read.csv("common_GOcategoriesH.csv", header=TRUE, sep=",",stringsAsFactors=FALSE)
# From HC:
# Get genes for every common category
relevantcats <- cattogenes[Hcat$GO.category] # This cuts down the CATEGORIES to only those in outpute
# Cut down GENES to only those initially present in the DE group
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighHC$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
geneHC <- enddf
write.csv(enddf, "commonGO_H_HC_genecat.csv", row.names=FALSE)
# What genes are these?
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = geneHC$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
# Append pLI scores
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
# Append dN/dS values
all.dndsr <- dndshighHC[c("Ensembl_gene_ID","humchimpdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
# With GO cats
colnames(geneHC) <- c("Category","Ensembl_gene_ID")
withterms <- eHC[c("GO.category","GO.term")]
colnames(withterms) <- c("Category","GO_term")
withterms <- merge(geneHC, withterms,by="Category")
knowgene <- merge(withterms, knowgene, by="Ensembl_gene_ID")
commongenesHC <- knowgene
# All genes in taste receptor activity category are the same as those in bitter taste receptor category
# Remove "taste receptor activity" category
commongenesHCx <- commongenesHC[commongenesHC$GO_term != "taste receptor activity",]

# Human vs orang
cattogenes <- goseq:::reversemapping(genetocatHO)
relevantcats <- cattogenes[Hcat$GO.category]
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighHO$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
geneHO <- enddf
write.csv(enddf, "commonGO_H_HO_genecat.csv", row.names=FALSE)
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = geneHO$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighHO[c("Ensembl_gene_ID","humorangdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
colnames(geneHO) <- c("Category","Ensembl_gene_ID")
withterms <- eHO[c("GO.category","GO.term")]
colnames(withterms) <- c("Category","GO_term")
withterms <- merge(geneHO, withterms,by="Category")
knowgene <- merge(withterms, knowgene, by="Ensembl_gene_ID")
commongenesHO <- knowgene
# All genes in taste receptor activity category are the same as those in bitter taste receptor category
# Except for the PIP gene, the sensory perception of bitter taste category is also redundant
# Remove redundant ones
commongenesHOx <- commongenesHO[commongenesHO$GO_term != "taste receptor activity",]
commongenesHOx <- subset(commongenesHOx, commongenesHOx$Gene_symbol=="PIP"|commongenesHOx$GO_term != "sensory perception of bitter taste")

# Human vs macaca
cattogenes <- goseq:::reversemapping(genetocatHM)
relevantcats <- cattogenes[Hcat$GO.category]
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighHM$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
geneHM <- enddf
write.csv(enddf, "commonGO_H_HM_genecat.csv", row.names=FALSE)
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = geneHM$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighHM[c("Ensembl_gene_ID","hummacdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
colnames(geneHM) <- c("Category","Ensembl_gene_ID")
withterms <- eHM[c("GO.category","GO.term")]
colnames(withterms) <- c("Category","GO_term")
withterms <- merge(geneHM, withterms,by="Category")
knowgene <- merge(withterms, knowgene, by="Ensembl_gene_ID")
commongenesHM <- knowgene
# All genes in taste receptor activity category are the same as those in bitter taste receptor category
# Remove "taste receptor activity" category
commongenesHMx <- commongenesHM[commongenesHM$GO_term != "taste receptor activity",]

commonH <- rbind(commongenesHCx,commongenesHOx,commongenesHMx)
# Need to take means of dN/dS values
meandnds <- ddply(commonH, "Ensembl_gene_ID", numcolwise(mean))
commonH$dndsr <- NULL
commonH$pLI <- NULL
commonH <- merge(commonH, meandnds, by="Ensembl_gene_ID")
commonH <- subset(commonH, duplicated(commonH$Ensembl_gene_ID)==FALSE)
write.csv(commonH, "common_genes_H.csv",row.names=FALSE)

# Chimpanzee
# Human vs chimp
Ccat <- read.csv("common_GOcategoriesC.csv", header=TRUE, sep=",",stringsAsFactors=FALSE)
cattogenes <- goseq:::reversemapping(genetocatHC) # Convert categories to genes
# From HC:
# Get genes for every common category
relevantcats <- cattogenes[Ccat$GO.category] # This cuts down the CATEGORIES to only those in outpute
# Cut down GENES to only those initially present in the DE group
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighHC$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
geneHC <- enddf
# What genes are these?
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = geneHC$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
# Append pLI scores
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
# Append dN/dS values
all.dndsr <- dndshighHC[c("Ensembl_gene_ID","humchimpdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
# With GO cats
colnames(geneHC) <- c("Category","Ensembl_gene_ID")
withterms <- eHC[c("GO.category","GO.term")]
colnames(withterms) <- c("Category","GO_term")
withterms <- merge(geneHC, withterms,by="Category")
knowgene <- merge(withterms, knowgene, by="Ensembl_gene_ID")
commongenesHC <- knowgene
# All genes in taste receptor activity category are the same as those in bitter taste receptor category
# Remove "taste receptor activity" category
commongenesHCy <- commongenesHC[commongenesHC$GO_term != "taste receptor activity",]

# Chimp vs orang
cattogenes <- goseq:::reversemapping(genetocatCO)
relevantcats <- cattogenes[Ccat$GO.category]
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighCO$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
geneCO <- enddf
write.csv(enddf, "commonGO_C_CO_genecat.csv", row.names=FALSE)
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = geneCO$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighCO[c("Ensembl_gene_ID","chimporangdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
colnames(geneCO) <- c("Category","Ensembl_gene_ID")
withterms <- eCO[c("GO.category","GO.term")]
colnames(withterms) <- c("Category","GO_term")
withterms <- merge(geneCO, withterms,by="Category")
knowgene <- merge(withterms, knowgene, by="Ensembl_gene_ID")
commongenesCO <- knowgene
commongenesCOx <- commongenesCO[commongenesCO$GO_term != "taste receptor activity",]

# Chimp vs macaca
cattogenes <- goseq:::reversemapping(genetocatCM)
relevantcats <- cattogenes[Ccat$GO.category]
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighCM$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
geneCM <- enddf
write.csv(enddf, "commonGO_C_CM_genecat.csv", row.names=FALSE)
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = geneCM$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighCM[c("Ensembl_gene_ID","chimpmacdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
colnames(geneCM) <- c("Category","Ensembl_gene_ID")
withterms <- eCM[c("GO.category","GO.term")]
colnames(withterms) <- c("Category","GO_term")
withterms <- merge(geneCM, withterms,by="Category")
knowgene <- merge(withterms, knowgene, by="Ensembl_gene_ID")
commongenesCM <- knowgene
# All genes in taste receptor activity category are the same as those in bitter taste receptor category
# Remove "taste receptor activity" category
commongenesCMx <- commongenesCM[commongenesCM$GO_term != "taste receptor activity",]

commonC <- rbind(commongenesHCy,commongenesCOx,commongenesCMx)
meandnds <- ddply(commonC, "Ensembl_gene_ID", numcolwise(mean))
commonC$dndsr <- NULL
commonC$pLI <- NULL
commonC <- merge(commonC, meandnds, by="Ensembl_gene_ID")
commonC <- subset(commonC, duplicated(commonC$Ensembl_gene_ID)==FALSE)
write.csv(commonC, "common_genes_C.csv",row.names=FALSE)

# Orangutan
# Human vs orang
Ocat <- read.csv("common_GOcategoriesO.csv", header=TRUE, sep=",",stringsAsFactors=FALSE)
cattogenes <- goseq:::reversemapping(genetocatHO)
relevantcats <- cattogenes[Ocat$GO.category]
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighHO$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
geneHO <- enddf
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = geneHO$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighHO[c("Ensembl_gene_ID","humorangdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
colnames(geneHO) <- c("Category","Ensembl_gene_ID")
withterms <- eHO[c("GO.category","GO.term")]
colnames(withterms) <- c("Category","GO_term")
withterms <- merge(geneHO, withterms,by="Category")
knowgene <- merge(withterms, knowgene, by="Ensembl_gene_ID")
commongenesHO <- knowgene
# All genes in taste receptor activity category are the same as those in bitter taste receptor category
# Remove "taste receptor activity" category
commongenesHOy <- commongenesHO[commongenesHO$GO_term != "taste receptor activity",]
commongenesHOy <- subset(commongenesHOy, commongenesHOy$Gene_symbol=="PIP"|commongenesHOy$GO_term != "sensory perception of bitter taste")

# Chimp vs orang
cattogenes <- goseq:::reversemapping(genetocatCO)
relevantcats <- cattogenes[Ocat$GO.category]
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighCO$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
geneCO <- enddf
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = geneCO$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighCO[c("Ensembl_gene_ID","chimporangdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
colnames(geneCO) <- c("Category","Ensembl_gene_ID")
withterms <- eCO[c("GO.category","GO.term")]
colnames(withterms) <- c("Category","GO_term")
withterms <- merge(geneCO, withterms,by="Category")
knowgene <- merge(withterms, knowgene, by="Ensembl_gene_ID")
commongenesCO <- knowgene
# All genes in taste receptor activity category are the same as those in bitter taste receptor category
# Except for the PIP and TAS2R60 genes, the sensory perception of bitter taste category is also redundant
# Remove redundant ones
commongenesCOy <- commongenesCO[commongenesCO$GO_term != "taste receptor activity",]
commongenesCOy <- subset(commongenesCOy, 
                         commongenesCOy$Gene_symbol=="PIP"|commongenesCOy$GO_term != "sensory perception of bitter taste"|commongenesCOy$Gene_symbol=="TAS2R60")

# Orang vs macaca
cattogenes <- goseq:::reversemapping(genetocatOM)
relevantcats <- cattogenes[Ocat$GO.category]
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighOM$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
geneOM <- enddf
write.csv(enddf, "commonGO_O_OM_genecat.csv", row.names=FALSE)
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = geneOM$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighOM[c("Ensembl_gene_ID","orangmacdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
colnames(geneOM) <- c("Category","Ensembl_gene_ID")
withterms <- eOM[c("GO.category","GO.term")]
colnames(withterms) <- c("Category","GO_term")
withterms <- merge(geneOM, withterms,by="Category")
knowgene <- merge(withterms, knowgene, by="Ensembl_gene_ID")
commongenesOM <- knowgene
# All genes in taste receptor activity category are the same as those in bitter taste receptor category
# Except for the PIP and TAS2R60 genes, the sensory perception of bitter taste category is also redundant
# Remove redundant ones
commongenesOMx <- commongenesOM[commongenesOM$GO_term != "taste receptor activity",]
commongenesOMx <- subset(commongenesOMx, 
                         commongenesOMx$Gene_symbol=="PIP"|commongenesOMx$GO_term != "sensory perception of bitter taste")

commonO <- rbind(commongenesHOy,commongenesCOy,commongenesOMx)
meandnds <- ddply(commonO, "Ensembl_gene_ID", numcolwise(mean))
commonO$dndsr <- NULL
commonO$pLI <- NULL
commonO <- merge(commonO, meandnds, by="Ensembl_gene_ID")
commonO <- subset(commonO, duplicated(commonO$Ensembl_gene_ID)==FALSE)
write.csv(commonO, "common_genes_O.csv",row.names=FALSE)

# Macaca
# Human vs macaca
Mcat <- read.csv("common_GOcategoriesM.csv", header=TRUE, sep=",",stringsAsFactors=FALSE)
cattogenes <- goseq:::reversemapping(genetocatHM)
relevantcats <- cattogenes[Mcat$GO.category]
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighHM$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
geneHM <- enddf
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = geneHM$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighHM[c("Ensembl_gene_ID","hummacdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
colnames(geneHM) <- c("Category","Ensembl_gene_ID")
withterms <- eHM[c("GO.category","GO.term")]
colnames(withterms) <- c("Category","GO_term")
withterms <- merge(geneHM, withterms,by="Category")
knowgene <- merge(withterms, knowgene, by="Ensembl_gene_ID")
commongenesHM <- knowgene
# All genes in taste receptor activity category are the same as those in bitter taste receptor category
# Remove "taste receptor activity" category
commongenesHMz <- commongenesHM[commongenesHM$GO_term != "taste receptor activity",]
commongenesHMz <- subset(commongenesHMz,
                         commongenesHMz$Gene_symbol=="PIP"|commongenesHMz$Gene_symbol=="IFI44"|commongenesHMz$GO_term=="immune response"|commongenesHMz$GO_term=="defense response"|commongenesHMz$GO_term=="bitter taste receptor activity")
dups <- commongenesHMz[duplicated(commongenesHMz$Gene_symbol)==TRUE,]
hold <- commongenesHMz[!((commongenesHMz$Ensembl_gene_ID %in% dups$Ensembl_gene_ID) & (commongenesHMz$Category=="GO:0006955")),]
hold <- hold[hold$Category!="GO:0009607",]
hold <- hold[hold$Category!="GO:0051707",]
hold <- hold[hold$Category!="GO:0050909",]
hold <- hold[hold$Category!="GO:0001580",]
hold <- hold[hold$Category!="GO:0050912",]
commongenesHMz <- hold

# Chimp vs macaca
cattogenes <- goseq:::reversemapping(genetocatCM)
relevantcats <- cattogenes[Mcat$GO.category]
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighCM$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
geneCM <- enddf
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = geneCM$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighCM[c("Ensembl_gene_ID","chimpmacdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
colnames(geneCM) <- c("Category","Ensembl_gene_ID")
withterms <- eCM[c("GO.category","GO.term")]
colnames(withterms) <- c("Category","GO_term")
withterms <- merge(geneCM, withterms,by="Category")
knowgene <- merge(withterms, knowgene, by="Ensembl_gene_ID")
commongenesCM <- knowgene
# All genes in taste receptor activity category are the same as those in bitter taste receptor category
# Remove "taste receptor activity" category
commongenesCMz <- commongenesCM[commongenesCM$GO_term != "taste receptor activity",]
dups <- commongenesCMz[duplicated(commongenesCMz$Gene_symbol)==TRUE,]
commongenesCMz <- commongenesCMz[!((commongenesCMz$Ensembl_gene_ID %in% dups$Ensembl_gene_ID) & (commongenesCMz$Category=="GO:0006955")),]
commongenesCMz <- commongenesCMz[!((commongenesCMz$Ensembl_gene_ID %in% dups$Ensembl_gene_ID) & (commongenesCMz$Category=="GO:0009607")),]
commongenesCMz <- commongenesCMz[!((commongenesCMz$Ensembl_gene_ID %in% dups$Ensembl_gene_ID) & (commongenesCMz$Category=="GO:0051707")),]
commongenesCMz <- commongenesCMz[!((commongenesCMz$Ensembl_gene_ID %in% dups$Ensembl_gene_ID) & (commongenesCMz$Category=="GO:0050912")),]
dups <- commongenesCMz[duplicated(commongenesCMz$Gene_symbol)==TRUE,]
commongenesCMz <- commongenesCMz[!((commongenesCMz$Ensembl_gene_ID %in% dups$Ensembl_gene_ID) & (commongenesCMz$Category=="GO:0098542")),]
commongenesCMz <- commongenesCMz[!((commongenesCMz$Ensembl_gene_ID %in% dups$Ensembl_gene_ID) & (commongenesCMz$Category=="GO:0050909")),]
commongenesCMz <- commongenesCMz[!((commongenesCMz$Ensembl_gene_ID %in% dups$Ensembl_gene_ID) & (commongenesCMz$Category=="GO:0001580")),]
commongenesCMz <- commongenesCMz[commongenesCMz$Category!="GO:0050913",]
dups <- commongenesCMz[duplicated(commongenesCMz$Gene_symbol)==TRUE,]
commongenesCMz <- commongenesCMz[!((commongenesCMz$Ensembl_gene_ID %in% dups$Ensembl_gene_ID) & (commongenesCMz$Category!="GO:0006952")),]

# Orang vs macaca
cattogenes <- goseq:::reversemapping(genetocatOM)
relevantcats <- cattogenes[Mcat$GO.category]
enddf <- data.frame()
catsnames <- names(relevantcats)
for(i in catsnames){
  y <- relevantcats[[i]]
  z <- y[y %in% dndshighOM$Ensembl_gene_ID]
  someholder <- data.frame(i, z)
  enddf <- rbind(enddf, someholder)
}
colnames(enddf) <- c("category", "gene")
geneOM <- enddf
knowgene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name","description"),
                  filters = "ensembl_gene_id", values = geneOM$gene, mart = marthuman)
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_Symbol", "Description")
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
knowgene.pli <- merge(knowgene, pli, by="Ensembl_gene_ID")
all.dndsr <- dndshighOM[c("Ensembl_gene_ID","orangmacdndsr")]
knowgene <- merge(knowgene.pli, all.dndsr, by="Ensembl_gene_ID")
colnames(knowgene) <- c("Ensembl_gene_ID", "Gene_symbol", "Description", "pLI", "dndsr")
colnames(geneOM) <- c("Category","Ensembl_gene_ID")
withterms <- eOM[c("GO.category","GO.term")]
colnames(withterms) <- c("Category","GO_term")
withterms <- merge(geneOM, withterms,by="Category")
knowgene <- merge(withterms, knowgene, by="Ensembl_gene_ID")
commongenesOM <- knowgene
# Remove redundant ones
commongenesOMy <- commongenesOM[commongenesOM$GO_term != "taste receptor activity",]
dups <- commongenesOMy[duplicated(commongenesOMy$Gene_symbol)==TRUE,]
wantedcats <- c("GO:0050913","GO:0006955","GO:0006952","GO:0033038")
commongenesOMy <- commongenesOMy[!((commongenesOMy$Ensembl_gene_ID %in% dups$Ensembl_gene_ID) & !(commongenesOMy$Category %in% wantedcats)),]
dups <- commongenesOMy[duplicated(commongenesOMy$Gene_symbol)==TRUE,]
commongenesOMy <- commongenesOMy[!((commongenesOMy$Ensembl_gene_ID %in% dups$Ensembl_gene_ID) & (commongenesOMy$Category=="GO:0006955")),]
dups <- commongenesOMy[duplicated(commongenesOMy$Gene_symbol)==TRUE,]
commongenesOMy <- commongenesOMy[!((commongenesOMy$Ensembl_gene_ID %in% dups$Ensembl_gene_ID) & (commongenesOMy$Category=="GO:0050913")),]
# IFI44 was left out. Put it back in.
commongenesOMy <- rbind(commongenesOMy,commongenesOM[commongenesOM$Category=="GO:0043207"&commongenesOM$Gene_symbol=="IFI44",])

commonM <- rbind(commongenesHMz,commongenesCMz,commongenesOMy)
meandnds <- ddply(commonM, "Ensembl_gene_ID", numcolwise(mean))
commonM$dndsr <- NULL
commonM$pLI <- NULL
commonM <- merge(commonM, meandnds, by="Ensembl_gene_ID")
commonM <- subset(commonM, duplicated(commonM$Ensembl_gene_ID)==FALSE)
write.csv(commonM, "common_genes_M.csv",row.names=FALSE)

# Finding common genes from lists
Hgenes <- read.csv("common_genes_H.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
Cgenes <- read.csv("common_genes_C.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
Ogenes <- read.csv("common_genes_O.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
Mgenes <- read.csv("common_genes_M.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

# Find matches
Hgenesx <- Hgenes
Cgenesx <- Cgenes
Ogenesx <- Ogenes
Mgenesx <- Mgenes
Hgenes$dndsr <- NULL
Cgenes$dndsr <- NULL
Ogenes$dndsr <- NULL
Mgenes$dndsr <- NULL
commongenes <- semi_join(Ogenes, Hgenes)
commongenes <- semi_join(commongenes, Cgenes)
commongenes <- semi_join(commongenes, Mgenes)

# Append means of dN/dS values
Hgenesy <- subset(Hgenesx, Hgenes$Ensembl_gene_ID %in% commongenes$Ensembl_gene_ID)
Cgenesy <- subset(Cgenesx, Cgenes$Ensembl_gene_ID %in% commongenes$Ensembl_gene_ID)
Ogenesy <- subset(Ogenesx, Ogenes$Ensembl_gene_ID %in% commongenes$Ensembl_gene_ID)
Mgenesy <- subset(Mgenesx, Mgenes$Ensembl_gene_ID %in% commongenes$Ensembl_gene_ID)
commongenesx <- rbind(Hgenesy, Cgenesy, Ogenesy, Mgenesy)
commongenesy <- ddply(commongenesx, "Ensembl_gene_ID", numcolwise(mean))
commongenesx$pLI <- NULL
commongenesx$dndsr <- NULL
commongenes <- merge(commongenesx, commongenesy, by="Ensembl_gene_ID")
commongenes <- subset(commongenes, duplicated(commongenes$Ensembl_gene_ID)==FALSE)

# 10 genes are common to all comparisons
write.csv(commongenes,"commonGO_commongenes.csv",row.names=FALSE)

############################################# Plot dN/dS vs pLI ##############################################

# Plot dN/dS vs pLI
# Append pLI and dN/dS values

# Original dataset
allexacgenes <- read.csv("exacgenes.csv", header = TRUE, sep =",")
hist(allexacgenes$pLI, main="All Genes", xlab="pLI", breaks=seq(0,1,0.1),ylim=(c(0,12000)))

# Enriched genes
allHC <- read.csv("enriched_HC_genecat.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
allHC <- allHC[,c("Ensembl_gene_ID","Gene_symbol")]
allHC <- subset(allHC, duplicated(allHC$Ensembl_gene_ID)==FALSE)
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
all.pli <- merge(allHC, pli, by="Ensembl_gene_ID")
dndsr <- dndshighHC[c("Ensembl_gene_ID","humchimpdndsr")]
toplot <- merge(dndsr, all.pli, by="Ensembl_gene_ID") # 176 genes
# Take out the most constrained genes (pLI >= 0.9)
constrainedHC <- subset(toplot, toplot$pLI >= 0.9) # 2 genes
rownames(constrainedHC) <- 1:nrow(constrainedHC)
colnames(constrainedHC) <- c("Ensembl_gene_ID","dN/dS","Gene_symbol","pLI")
constrainedHC <- constrainedHC[c("Ensembl_gene_ID","Gene_symbol","dN/dS","pLI")]
smoothScatter(toplot$pLI, toplot$humchimpdndsr, main="Human vs Chimp",xlab="pLI",ylab="dN/dS",ylim=(c(1,10)))
hist(toplot$pLI, main="Human vs Chimp", xlab="pLI", breaks=seq(0,1,0.1))

allHO <- read.csv("enriched_HO_genecat.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
allHO <- allHO[,c("Ensembl_gene_ID","Gene_symbol")]
allHO <- subset(allHO, duplicated(allHO$Ensembl_gene_ID)==FALSE)
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
all.pli <- merge(allHO, pli, by="Ensembl_gene_ID")
dndsr <- dndshighHO[c("Ensembl_gene_ID","humorangdndsr")]
toplot <- merge(dndsr, all.pli, by="Ensembl_gene_ID") # 41 genes
# Take out the most constrained genes (pLI >= 0.9)
constrainedHO <- subset(toplot, toplot$pLI >= 0.9) # 0 genes
rownames(constrainedHO) <- 1:nrow(constrainedHO)
colnames(constrainedHO) <- c("Ensembl_gene_ID","dN/dS","Gene_symbol","pLI")
constrainedHO <- constrainedHO[c("Ensembl_gene_ID","Gene_symbol","dN/dS","pLI")]
smoothScatter(toplot$pLI, toplot$humorangdndsr, main="Human vs Orangutan",xlab="pLI",ylab="dN/dS",ylim=(c(1,10)))
hist(toplot$pLI, main="Human vs Orangutan", xlab="pLI", breaks=seq(0,1,0.1))

allHM <- read.csv("enriched_HM_genecat.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
allHM <- allHM[,c("Ensembl_gene_ID","Gene_symbol")]
allHM <- subset(allHM, duplicated(allHM$Ensembl_gene_ID)==FALSE)
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
all.pli <- merge(allHM, pli, by="Ensembl_gene_ID")
dndsr <- dndshighHM[c("Ensembl_gene_ID","hummacdndsr")]
toplot <- merge(dndsr, all.pli, by="Ensembl_gene_ID") # 101 genes
# Take out the most constrained genes (pLI >= 0.9)
constrainedHM <- subset(toplot, toplot$pLI >= 0.9) #0
rownames(constrainedHM) <- 1:nrow(constrainedHM)
colnames(constrainedHM) <- c("Ensembl_gene_ID","dN/dS","Gene_symbol","pLI")
constrainedHM <- constrainedHM[c("Ensembl_gene_ID","Gene_symbol","dN/dS","pLI")]
smoothScatter(toplot$pLI, toplot$hummacdndsr, main="Human vs Macaque",xlab="pLI",ylab="dN/dS",ylim=(c(1,10)))
hist(toplot$pLI, main="Human vs Macaque", xlab="pLI", breaks=seq(0,1,0.1))

allCO <- read.csv("enriched_CO_genecat.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
allCO <- allCO[,c("Ensembl_gene_ID","Gene_symbol")]
allCO <- subset(allCO, duplicated(allCO$Ensembl_gene_ID)==FALSE)
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
all.pli <- merge(allCO, pli, by="Ensembl_gene_ID")
dndsr <- dndshighCO[c("Ensembl_gene_ID","chimporangdndsr")]
toplot <- merge(dndsr, all.pli, by="Ensembl_gene_ID") # 12 genes
# Take out the most constrained genes (pLI >= 0.9)
constrainedCO <- subset(toplot, toplot$pLI >= 0.9) # 0
rownames(constrainedCO) <- 1:nrow(constrainedCO)
colnames(constrainedCO) <- c("Ensembl_gene_ID","dN/dS","Gene_symbol","pLI")
constrainedCO <- constrainedCO[c("Ensembl_gene_ID","Gene_symbol","dN/dS","pLI")]
smoothScatter(toplot$pLI, toplot$chimporangdndsr, main="Chimpanzee vs Orangutan",xlab="pLI",ylab="dN/dS",ylim=(c(1,10)))
hist(toplot$pLI, main="Chimpanzee vs Orangutan", xlab="pLI", breaks=seq(0,1,0.1))

allCM <- read.csv("enriched_CM_genecat.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
allCM <- allCM[,c("Ensembl_gene_ID","Gene_symbol")]
allCM <- subset(allCM, duplicated(allCM$Ensembl_gene_ID)==FALSE)
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
all.pli <- merge(allCM, pli, by="Ensembl_gene_ID")
dndsr <- dndshighCM[c("Ensembl_gene_ID","chimpmacdndsr")]
toplot <- merge(dndsr, all.pli, by="Ensembl_gene_ID") # 93 genes
# Take out the most constrained genes (pLI >= 0.9)
constrainedCM <- subset(toplot, toplot$pLI >= 0.9) # 1
rownames(constrainedCM) <- 1:nrow(constrainedCM)
colnames(constrainedCM) <- c("Ensembl_gene_ID","dN/dS","Gene_symbol","pLI")
constrainedCM <- constrainedCM[c("Ensembl_gene_ID","Gene_symbol","dN/dS","pLI")]
smoothScatter(toplot$pLI, toplot$chimpmacdndsr, main="Chimpanzee vs Macaque",xlab="pLI",ylab="dN/dS",ylim=(c(1,10)))
hist(toplot$pLI, main="Chimpanzee vs Macaque", xlab="pLI", breaks=seq(0,1,0.1))

allOM <- read.csv("enriched_OM_genecat.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
allOM <- allOM[,c("Ensembl_gene_ID","Gene_symbol")]
allOM <- subset(allOM, duplicated(allOM$Ensembl_gene_ID)==FALSE)
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
all.pli <- merge(allOM, pli, by="Ensembl_gene_ID")
dndsr <- dndshighOM[c("Ensembl_gene_ID","orangmacdndsr")]
toplot <- merge(dndsr, all.pli, by="Ensembl_gene_ID") # 69 genes
# Take out the most constrained genes (pLI >= 0.9)
constrainedOM <- subset(toplot, toplot$pLI >= 0.9) # 0
rownames(constrainedOM) <- 1:nrow(constrainedOM)
colnames(constrainedOM) <- c("Ensembl_gene_ID","dN/dS","Gene_symbol","pLI")
constrainedOM <- constrainedOM[c("Ensembl_gene_ID","Gene_symbol","dN/dS","pLI")]
smoothScatter(toplot$pLI, toplot$orangmacdndsr, main="Orangutan vs Macaque",xlab="pLI",ylab="dN/dS",ylim=(c(1,10)))
hist(toplot$pLI, main="Orangutan vs Macaque", xlab="pLI", breaks=seq(0,1,0.1))

# Random sampling test for genes with high dN/dS and constraint scores
allexacgenes <- read.csv("exacgenes.csv", header = TRUE, sep =",")
rsample <- function(x){
  rs <- sample(allexacgenes$pLI, x, replace=F)
  highdnds <- subset(rs, rs >=0.9)
  return(length(highdnds))
}
pv <- function(a,b){
  test <- replicate(1000000, rsample(a))
  counts <- sum(test > b)
  p <- counts/1000000 
  return(p)
}
pvless <- function(a,b){
  test <- replicate(1000000, rsample(a))
  counts <- sum(test <= b)
  p <- counts/1000000 
  return(p)
}

# For goseq data
# Human vs Chimp
pvless(176,2) # < 1e-06

# Human vs Orang
pvless(41,0) # 0.000327

# Human vs Macaca
pvless(101,0) # < 1e-06

# Chimp vs Orang
pvless(12,0) # 0.096282

# Chimp vs Macaca
pvless(93,1) # 0

# Orang vs Macaca
pvless(69,0) # 1e-06

# Genes common among all comparisons
allconstrained <- rbind(constrainedHC, constrainedHO, constrainedHM, constrainedCO, constrainedCM, constrainedOM)
repeats <- ddply(allconstrained,c("Ensembl_gene_ID","Gene_symbol"),nrow)
allcommon <- subset(repeats, repeats$V1=="6")
rownames(allcommon) <- 1:nrow(allcommon)

# Common among all human
allhuman <- rbind(constrainedHC, constrainedHM, constrainedHO)
repeats <- ddply(allhuman,c("Ensembl_gene_ID","Gene_symbol"),nrow)
allhuman <- subset(repeats, repeats$V1=="3")
rownames(allhuman) <- 1:nrow(allhuman)

# Common among all chimp
allchimp <- rbind(constrainedHC, constrainedCM, constrainedCO)
repeats <- ddply(allchimp,c("Ensembl_gene_ID","Gene_symbol"),nrow)
allchimp <- subset(repeats, repeats$V1=="3")
rownames(allchimp) <- 1:nrow(allchimp)

# Common among all orangutan
allorang <- rbind(constrainedHO, constrainedOM, constrainedCO)
repeats <- ddply(allorang,c("Ensembl_gene_ID","Gene_symbol"),nrow)
allorang <- subset(repeats, repeats$V1=="3")
rownames(allorang) <- 1:nrow(allorang)

# Common among all macaque
allmac <- rbind(constrainedHM, constrainedOM, constrainedCM)
repeats <- ddply(allmac,c("Ensembl_gene_ID","Gene_symbol"),nrow)
allmac <- subset(repeats, repeats$V1=="3")
rownames(allmac) <- 1:nrow(allmac)

############################################# Attempting PANTHER #############################################
setwd("~/ExAC_dataset_project/csv_records")

dndshighHC <- read.csv("dndshighHC.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
dndshighHO <- read.csv("dndshighHO.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
dndshighHM <- read.csv("dndshighHM.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
dndshighCO <- read.csv("dndshighCO.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
dndshighCM <- read.csv("dndshighCM.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
dndshighOM <- read.csv("dndshighOM.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

# Input files (the "DE"s)
dndshighHC <- dndshighHC["Ensembl_gene_ID"]
colnames(dndshighHC) <- NULL
dndshighHO <- dndshighHO["Ensembl_gene_ID"]
colnames(dndshighHO) <- NULL
dndshighHM <- dndshighHM["Ensembl_gene_ID"]
colnames(dndshighHM) <- NULL
dndshighCO <- dndshighCO["Ensembl_gene_ID"]
colnames(dndshighCO) <- NULL
dndshighCM <- dndshighCM["Ensembl_gene_ID"]
colnames(dndshighCM) <- NULL
dndshighOM <- dndshighOM["Ensembl_gene_ID"]
colnames(dndshighOM) <- NULL

setwd("~/ExAC_dataset_project/PANTHER")
write.table(dndshighHC, "dndshighHC_ID.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(dndshighHO, "dndshighHO_ID.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(dndshighHM, "dndshighHM_ID.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(dndshighCO, "dndshighCO_ID.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(dndshighCM, "dndshighCM_ID.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(dndshighOM, "dndshighOM_ID.txt", quote=FALSE, row.names=FALSE, sep="\t")

# Background file
sghcID <- sghc["Ensembl_gene_ID"]
colnames(sghcID) <- NULL
sghoID <- sgho["Ensembl_gene_ID"]
colnames(sghoID) <- NULL
sghmID <- sghm["Ensembl_gene_ID"]
colnames(sghmID) <- NULL
sgcoID <- sgco["Ensembl_gene_ID"]
colnames(sgcoID) <- NULL
sgcmID <- sgcm["Ensembl_gene_ID"]
colnames(sgcmID) <- NULL
sgomID <- sgom["Ensembl_gene_ID"]
colnames(sgomID) <- NULL

write.table(sghcID, "IDsghc.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(sghoID, "IDsgho.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(sghmID, "IDsghm.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(sgcoID, "IDsgco.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(sgcmID, "IDsgcm.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(sgomID, "IDsgom.txt", quote=FALSE, row.names=FALSE, sep="\t")

# Human vs Chimpanzee
# Taste
pantherHClist1 <- read.table("PANTHER_HC_list1.txt", header=FALSE, sep="\t")
pantherHClist1 <- data.frame(pantherHClist1)
colnames(pantherHClist1) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherHClist1 <- pantherHClist1[c("Ensembl_gene_ID", "Gene_symbol")]
pantherHClist1$Gene_symbol <- as.character(pantherHClist1$Gene_symbol)
# Removing the ";ortholog" from the ends
pantherHClist1$Gene_symbol <- sub("\\;ortholog", "\\", pantherHClist1$Gene_symbol)
# Keep only the gene symbol
pantherHClist1$Gene_symbol <- sub(".*;", "", pantherHClist1$Gene_symbol)

# Olfactory
pantherHClist2 <- read.table("PANTHER_HC_list2.txt", header=FALSE, sep="\t")
pantherHClist2 <- data.frame(pantherHClist2)
colnames(pantherHClist2) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherHClist2 <- pantherHClist2[c("Ensembl_gene_ID", "Gene_symbol")]
pantherHClist2$Gene_symbol <- as.character(pantherHClist2$Gene_symbol)
pantherHClist2$Gene_symbol <- sub("\\;ortholog", "\\", pantherHClist2$Gene_symbol)
pantherHClist2$Gene_symbol <- sub(".*;", "", pantherHClist2$Gene_symbol)

pantherHClist3 <- read.table("PANTHER_HC_list3.txt", header=FALSE, sep="\t")
pantherHClist3 <- data.frame(pantherHClist3)
colnames(pantherHClist3) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherHClist3 <- pantherHClist3[c("Ensembl_gene_ID", "Gene_symbol")]
pantherHClist3$Gene_symbol <- as.character(pantherHClist3$Gene_symbol)
pantherHClist3$Gene_symbol <- sub("\\;ortholog", "\\", pantherHClist3$Gene_symbol)
pantherHClist3$Gene_symbol <- sub(".*;", "", pantherHClist3$Gene_symbol)
pantherHClist3 <- anti_join(pantherHClist3, pantherHClist1, by="Ensembl_gene_ID")
pantherHClist3 <- anti_join(pantherHClist3, pantherHClist2, by="Ensembl_gene_ID")

totalHC <- rbind(pantherHClist1,pantherHClist2,pantherHClist3)

# Human vs Orangutan
# Taste
pantherHOlist1 <- read.table("PANTHER_HO_list1.txt", header=FALSE, sep="\t")
pantherHOlist1 <- data.frame(pantherHOlist1)
colnames(pantherHOlist1) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherHOlist1 <- pantherHOlist1[c("Ensembl_gene_ID", "Gene_symbol")]
pantherHOlist1$Gene_symbol <- as.character(pantherHOlist1$Gene_symbol)
pantherHOlist1$Gene_symbol <- sub("\\;ortholog", "\\", pantherHOlist1$Gene_symbol)
pantherHOlist1$Gene_symbol <- sub(".*;", "", pantherHOlist1$Gene_symbol)

# Immunity
pantherHOlist2 <- read.table("PANTHER_HO_list2.txt", header=FALSE, sep="\t")
pantherHOlist2 <- data.frame(pantherHOlist2)
colnames(pantherHOlist2) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherHOlist2 <- pantherHOlist2[c("Ensembl_gene_ID", "Gene_symbol")]
pantherHOlist2$Gene_symbol <- as.character(pantherHOlist2$Gene_symbol)
pantherHOlist2$Gene_symbol <- sub("\\;ortholog", "\\", pantherHOlist2$Gene_symbol)
pantherHOlist2$Gene_symbol <- sub(".*;", "", pantherHOlist2$Gene_symbol)

pantherHOlist4 <- read.table("PANTHER_HO_list4.txt", header=FALSE, sep="\t")
pantherHOlist4 <- data.frame(pantherHOlist4)
colnames(pantherHOlist4) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherHOlist4 <- pantherHOlist4[c("Ensembl_gene_ID", "Gene_symbol")]
pantherHOlist4$Gene_symbol <- as.character(pantherHOlist4$Gene_symbol)
pantherHOlist4$Gene_symbol <- sub("\\;ortholog", "\\", pantherHOlist4$Gene_symbol)
pantherHOlist4$Gene_symbol <- sub(".*;", "", pantherHOlist4$Gene_symbol)

pantherHOlist5 <- read.table("PANTHER_HO_list5.txt", header=FALSE, sep="\t")
pantherHOlist5 <- data.frame(pantherHOlist5)
colnames(pantherHOlist5) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherHOlist5 <- pantherHOlist5[c("Ensembl_gene_ID", "Gene_symbol")]
pantherHOlist5$Gene_symbol <- as.character(pantherHOlist5$Gene_symbol)
pantherHOlist5$Gene_symbol <- sub("\\;ortholog", "\\", pantherHOlist5$Gene_symbol)
pantherHOlist5$Gene_symbol <- sub(".*;", "", pantherHOlist5$Gene_symbol)

pantherHOlistI <- rbind(pantherHOlist2, pantherHOlist4, pantherHOlist5)
pantherHOlistI <- subset(pantherHOlistI, duplicated(pantherHOlistI$Ensembl_gene_ID)==FALSE)
rownames(pantherHOlistI) <- 1:nrow(pantherHOlistI)

# Olfactory receptors
pantherHOlist3 <- read.table("PANTHER_HO_list3.txt", header=FALSE, sep="\t")
pantherHOlist3 <- data.frame(pantherHOlist3)
colnames(pantherHOlist3) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherHOlist3 <- pantherHOlist3[c("Ensembl_gene_ID", "Gene_symbol")]
pantherHOlist3$Gene_symbol <- as.character(pantherHOlist3$Gene_symbol)
pantherHOlist3$Gene_symbol <- sub("\\;ortholog", "\\", pantherHOlist3$Gene_symbol)
pantherHOlist3$Gene_symbol <- sub(".*;", "", pantherHOlist3$Gene_symbol)

totalHO <- rbind(pantherHOlist1,pantherHOlistI,pantherHOlist3)

# Human vs Macaque
# Taste
pantherHMlist1 <- read.table("PANTHER_HM_list1.txt", header=FALSE, sep="\t")
pantherHMlist1 <- data.frame(pantherHMlist1)
colnames(pantherHMlist1) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherHMlist1 <- pantherHMlist1[c("Ensembl_gene_ID", "Gene_symbol")]
pantherHMlist1$Gene_symbol <- as.character(pantherHMlist1$Gene_symbol)
pantherHMlist1$Gene_symbol <- sub("\\;ortholog", "\\", pantherHMlist1$Gene_symbol)
pantherHMlist1$Gene_symbol <- sub(".*;", "", pantherHMlist1$Gene_symbol)

# Energy
pantherHMlist3 <- read.table("PANTHER_HM_list3.txt", header=FALSE, sep="\t")
pantherHMlist3 <- data.frame(pantherHMlist3)
colnames(pantherHMlist3) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherHMlist3 <- pantherHMlist3[c("Ensembl_gene_ID", "Gene_symbol")]
pantherHMlist3$Gene_symbol <- as.character(pantherHMlist3$Gene_symbol)
pantherHMlist3$Gene_symbol <- sub("\\;ortholog", "\\", pantherHMlist3$Gene_symbol)
pantherHMlist3$Gene_symbol <- sub(".*;", "", pantherHMlist3$Gene_symbol)

# Immunity genes
pantherHMlist2 <- read.table("PANTHER_HM_list2.txt", header=FALSE, sep="\t")
pantherHMlist2 <- data.frame(pantherHMlist2)
colnames(pantherHMlist2) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherHMlist2 <- pantherHMlist2[c("Ensembl_gene_ID", "Gene_symbol")]
pantherHMlist2$Gene_symbol <- as.character(pantherHMlist2$Gene_symbol)
pantherHMlist2$Gene_symbol <- sub("\\;ortholog", "\\", pantherHMlist2$Gene_symbol)
pantherHMlist2$Gene_symbol <- sub(".*;", "", pantherHMlist2$Gene_symbol)

pantherHMlist4 <- read.table("PANTHER_HM_list4.txt", header=FALSE, sep="\t")
pantherHMlist4 <- data.frame(pantherHMlist4)
colnames(pantherHMlist4) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherHMlist4 <- pantherHMlist4[c("Ensembl_gene_ID", "Gene_symbol")]
pantherHMlist4$Gene_symbol <- as.character(pantherHMlist4$Gene_symbol)
pantherHMlist4$Gene_symbol <- sub("\\;ortholog", "\\", pantherHMlist4$Gene_symbol)
pantherHMlist4$Gene_symbol <- sub(".*;", "", pantherHMlist4$Gene_symbol)

pantherHMlist5 <- read.table("PANTHER_HM_list5.txt", header=FALSE, sep="\t")
pantherHMlist5 <- data.frame(pantherHMlist5)
colnames(pantherHMlist5) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherHMlist5 <- pantherHMlist5[c("Ensembl_gene_ID", "Gene_symbol")]
pantherHMlist5$Gene_symbol <- as.character(pantherHMlist5$Gene_symbol)
pantherHMlist5$Gene_symbol <- sub("\\;ortholog", "\\", pantherHMlist5$Gene_symbol)
pantherHMlist5$Gene_symbol <- sub(".*;", "", pantherHMlist5$Gene_symbol)

pantherHMlistI <- rbind(pantherHMlist2, pantherHMlist4, pantherHMlist5)
pantherHMlistI <- subset(pantherHMlistI, duplicated(pantherHMlistI$Ensembl_gene_ID)==FALSE)
rownames(pantherHMlistI) <- 1:nrow(pantherHMlistI)

totalHM <- rbind(pantherHMlist1,pantherHMlist3,pantherHMlistI)

# Chimp vs Orangutan
# Taste
pantherCOlist1 <- read.table("~/ExAC_dataset_project/PANTHER/PANTHER_CO_list1.txt", header=FALSE, sep="\t")
pantherCOlist1 <- data.frame(pantherCOlist1)
colnames(pantherCOlist1) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherCOlist1 <- pantherCOlist1[c("Ensembl_gene_ID", "Gene_symbol")]
pantherCOlist1$Gene_symbol <- as.character(pantherCOlist1$Gene_symbol)
pantherCOlist1$Gene_symbol <- sub("\\;ortholog", "\\", pantherCOlist1$Gene_symbol)
pantherCOlist1$Gene_symbol <- sub(".*;", "", pantherCOlist1$Gene_symbol)

# Immunity
pantherCOlist2 <- read.table("~/ExAC_dataset_project/PANTHER/PANTHER_CO_list2.txt", header=FALSE, sep="\t")
pantherCOlist2 <- data.frame(pantherCOlist2)
colnames(pantherCOlist2) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherCOlist2 <- pantherCOlist2[c("Ensembl_gene_ID", "Gene_symbol")]
pantherCOlist2$Gene_symbol <- as.character(pantherCOlist2$Gene_symbol)
pantherCOlist2$Gene_symbol <- sub("\\;ortholog", "\\", pantherCOlist2$Gene_symbol)
pantherCOlist2$Gene_symbol <- sub(".*;", "", pantherCOlist2$Gene_symbol)

pantherCOlist3 <- read.table("~/ExAC_dataset_project/PANTHER/PANTHER_CO_list3.txt", header=FALSE, sep="\t")
pantherCOlist3 <- data.frame(pantherCOlist3)
colnames(pantherCOlist3) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherCOlist3 <- pantherCOlist3[c("Ensembl_gene_ID", "Gene_symbol")]
pantherCOlist3$Gene_symbol <- as.character(pantherCOlist3$Gene_symbol)
pantherCOlist3$Gene_symbol <- sub("\\;ortholog", "\\", pantherCOlist3$Gene_symbol)
pantherCOlist3$Gene_symbol <- sub(".*;", "", pantherCOlist3$Gene_symbol)

pantherCOlistI <- rbind(pantherCOlist2, pantherCOlist3)
pantherCOlistI <- subset(pantherCOlistI, duplicated(pantherCOlistI$Ensembl_gene_ID)==FALSE)
rownames(pantherCOlistI) <- 1:nrow(pantherCOlistI)

totalCO <- rbind(pantherCOlist1,pantherCOlistI)

# Chimp vs Macaque
# Taste
pantherCMlist1 <- read.table("~/ExAC_dataset_project/PANTHER/PANTHER_CM_list1.txt", header=FALSE, sep="\t")
pantherCMlist1 <- data.frame(pantherCMlist1)
colnames(pantherCMlist1) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherCMlist1 <- pantherCMlist1[c("Ensembl_gene_ID", "Gene_symbol")]
pantherCMlist1$Gene_symbol <- as.character(pantherCMlist1$Gene_symbol)
pantherCMlist1$Gene_symbol <- sub("\\;ortholog", "\\", pantherCMlist1$Gene_symbol)
pantherCMlist1$Gene_symbol <- sub(".*;", "", pantherCMlist1$Gene_symbol)

# Immunity
pantherCMlist2 <- read.table("~/ExAC_dataset_project/PANTHER/PANTHER_CM_list2.txt", header=FALSE, sep="\t")
pantherCMlist2 <- data.frame(pantherCMlist2)
colnames(pantherCMlist2) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherCMlist2 <- pantherCMlist2[c("Ensembl_gene_ID", "Gene_symbol")]
pantherCMlist2$Gene_symbol <- as.character(pantherCMlist2$Gene_symbol)
pantherCMlist2$Gene_symbol <- sub("\\;ortholog", "\\", pantherCMlist2$Gene_symbol)
pantherCMlist2$Gene_symbol <- sub(".*;", "", pantherCMlist2$Gene_symbol)

pantherCMlist3 <- read.table("~/ExAC_dataset_project/PANTHER/PANTHER_CM_list3.txt", header=FALSE, sep="\t")
pantherCMlist3 <- data.frame(pantherCMlist3)
colnames(pantherCMlist3) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherCMlist3 <- pantherCMlist3[c("Ensembl_gene_ID", "Gene_symbol")]
pantherCMlist3$Gene_symbol <- as.character(pantherCMlist3$Gene_symbol)
pantherCMlist3$Gene_symbol <- sub("\\;ortholog", "\\", pantherCMlist3$Gene_symbol)
pantherCMlist3$Gene_symbol <- sub(".*;", "", pantherCMlist3$Gene_symbol)

pantherCMlistI <- rbind(pantherCMlist2, pantherCMlist3)
pantherCMlistI <- subset(pantherCMlistI, duplicated(pantherCMlistI$Ensembl_gene_ID)==FALSE)
rownames(pantherCMlistI) <- 1:nrow(pantherCMlistI)

totalCM <- rbind(pantherCMlist1,pantherCMlistI)

# Orangutan vs Macaque
# Taste
pantherOMlist2 <- read.table("~/ExAC_dataset_project/PANTHER/PANTHER_OM_list2.txt", header=FALSE, sep="\t")
pantherOMlist2 <- data.frame(pantherOMlist2)
colnames(pantherOMlist2) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherOMlist2 <- pantherOMlist2[c("Ensembl_gene_ID", "Gene_symbol")]
pantherOMlist2$Gene_symbol <- as.character(pantherOMlist2$Gene_symbol)
pantherOMlist2$Gene_symbol <- sub("\\;ortholog", "\\", pantherOMlist2$Gene_symbol)
pantherOMlist2$Gene_symbol <- sub(".*;", "", pantherOMlist2$Gene_symbol)

# Immunity
pantherOMlist1 <- read.table("~/ExAC_dataset_project/PANTHER/PANTHER_OM_list1.txt", header=FALSE, sep="\t")
pantherOMlist1 <- data.frame(pantherOMlist1)
colnames(pantherOMlist1) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherOMlist1 <- pantherOMlist1[c("Ensembl_gene_ID", "Gene_symbol")]
pantherOMlist1$Gene_symbol <- as.character(pantherOMlist1$Gene_symbol)
pantherOMlist1$Gene_symbol <- sub("\\;ortholog", "\\", pantherOMlist1$Gene_symbol)
pantherOMlist1$Gene_symbol <- sub(".*;", "", pantherOMlist1$Gene_symbol)

pantherOMlist3 <- read.table("~/ExAC_dataset_project/PANTHER/PANTHER_OM_list3.txt", header=FALSE, sep="\t")
pantherOMlist3 <- data.frame(pantherOMlist3)
colnames(pantherOMlist3) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherOMlist3 <- pantherOMlist3[c("Ensembl_gene_ID", "Gene_symbol")]
pantherOMlist3$Gene_symbol <- as.character(pantherOMlist3$Gene_symbol)
pantherOMlist3$Gene_symbol <- sub("\\;ortholog", "\\", pantherOMlist3$Gene_symbol)
pantherOMlist3$Gene_symbol <- sub(".*;", "", pantherOMlist3$Gene_symbol)

pantherOMlist5 <- read.table("~/ExAC_dataset_project/PANTHER/PANTHER_OM_list5.txt", header=FALSE, sep="\t")
pantherOMlist5 <- data.frame(pantherOMlist5)
colnames(pantherOMlist5) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherOMlist5 <- pantherOMlist5[c("Ensembl_gene_ID", "Gene_symbol")]
pantherOMlist5$Gene_symbol <- as.character(pantherOMlist5$Gene_symbol)
pantherOMlist5$Gene_symbol <- sub("\\;ortholog", "\\", pantherOMlist5$Gene_symbol)
pantherOMlist5$Gene_symbol <- sub(".*;", "", pantherOMlist5$Gene_symbol)

pantherOMlist6 <- read.table("~/ExAC_dataset_project/PANTHER/PANTHER_OM_list6.txt", header=FALSE, sep="\t")
pantherOMlist6 <- data.frame(pantherOMlist6)
colnames(pantherOMlist6) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherOMlist6 <- pantherOMlist6[c("Ensembl_gene_ID", "Gene_symbol")]
pantherOMlist6$Gene_symbol <- as.character(pantherOMlist6$Gene_symbol)
pantherOMlist6$Gene_symbol <- sub("\\;ortholog", "\\", pantherOMlist6$Gene_symbol)
pantherOMlist6$Gene_symbol <- sub(".*;", "", pantherOMlist6$Gene_symbol)

pantherOMlistI <- rbind(pantherOMlist1, pantherOMlist3, pantherOMlist5, pantherOMlist6)
pantherOMlistI <- subset(pantherOMlistI, duplicated(pantherOMlistI$Ensembl_gene_ID)==FALSE)
rownames(pantherOMlistI) <- 1:nrow(pantherOMlistI)

# Energy
pantherOMlist4 <- read.table("~/ExAC_dataset_project/PANTHER/PANTHER_OM_list4.txt", header=FALSE, sep="\t")
pantherOMlist4 <- data.frame(pantherOMlist4)
colnames(pantherOMlist4) <- c("UniProt", "Ensembl_gene_ID", "Gene_symbol")
pantherOMlist4 <- pantherOMlist4[c("Ensembl_gene_ID", "Gene_symbol")]
pantherOMlist4$Gene_symbol <- as.character(pantherOMlist4$Gene_symbol)
pantherOMlist4$Gene_symbol <- sub("\\;ortholog", "\\", pantherOMlist4$Gene_symbol)
pantherOMlist4$Gene_symbol <- sub(".*;", "", pantherOMlist4$Gene_symbol)

totalOM <- rbind(pantherOMlist2,pantherOMlistI,pantherOMlist4)

# Common genes
# Genes that consistently show up as enriched in all comparisons
Giantlist <- rbind(totalHC,totalHO,totalHM,totalCO,totalCM,totalOM)
countrepeats <- ddply(Giantlist,c("Ensembl_gene_ID","Gene_symbol"),nrow)
allcomgenes <- subset(countrepeats, countrepeats$V1=="6")
rownames(allcomgenes) <- 1:nrow(allcomgenes)

Humanlist <- rbind(totalHC,totalHO,totalHM)
Chimplist <- rbind(totalHC,totalCO,totalCM)
Oranglist <- rbind(totalHO,totalCO,totalOM)
Maclist <- rbind(totalHM,totalCM,totalOM)

# Genes only in human comparisons
Humancomgenes <- anti_join(Humanlist,Chimplist)
Humancomgenes <- anti_join(Humancomgenes,Oranglist)
Humancomgenes <- anti_join(Humancomgenes,Maclist) # None

# Genes only in chimp comparisons
Chimpcomgenes <- anti_join(Chimplist,Humanlist)
Chimpcomgenes <- anti_join(Chimpcomgenes,Oranglist)
Chimpcomgenes <- anti_join(Chimpcomgenes,Maclist) # None

# Genes only in orangutan comparisons
Orangcomgenes <- anti_join(Oranglist,Humanlist)
Orangcomgenes <- anti_join(Orangcomgenes,Chimplist)
Orangcomgenes <- anti_join(Orangcomgenes,Maclist) # None

# Genes only in macaca comparisons
Maccomgenes <- anti_join(Maclist,Humanlist)
Maccomgenes <- anti_join(Maccomgenes,Chimplist)
Maccomgenes <- anti_join(Maccomgenes,Oranglist) # None

# Genes common to human comparisons
countrepeats <- ddply(Humanlist,c("Ensembl_gene_ID","Gene_symbol"),nrow)
allhumgenes <- subset(countrepeats, countrepeats$V1=="3")
rownames(allhumgenes) <- 1:nrow(allhumgenes)

# Genes common to chimp comparisons
countrepeats <- ddply(Chimplist,c("Ensembl_gene_ID","Gene_symbol"),nrow)
allchimpgenes <- subset(countrepeats, countrepeats$V1=="3")
rownames(allchimpgenes) <- 1:nrow(allchimpgenes)

# Genes common to orang comparisons
countrepeats <- ddply(Oranglist,c("Ensembl_gene_ID","Gene_symbol"),nrow)
alloranggenes <- subset(countrepeats, countrepeats$V1=="3")
rownames(alloranggenes) <- 1:nrow(alloranggenes)

# Genes common to macaca comparisons
countrepeats <- ddply(Maclist,c("Ensembl_gene_ID","Gene_symbol"),nrow)
allmacgenes <- subset(countrepeats, countrepeats$V1=="3")
rownames(allmacgenes) <- 1:nrow(allmacgenes)

# Plot dN/dS against pLI scores for these genes
setwd("~/ExAC_dataset_project/csv_records")
# Human vs chimp
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
totalHC.pli <- merge(totalHC, pli, by="Ensembl_gene_ID")
sghc <- read.csv("sghc.csv",header=TRUE,stringsAsFactors = FALSE,sep=",")
dnds <- sghc[c("Ensembl_gene_ID","humchimpdndsr")]
toplot <- merge(totalHC.pli,dnds,by="Ensembl_gene_ID") # 125
allHC <- toplot
# Take out the most constrained genes (pLI >= 0.9)
constrainedHC <- subset(toplot, toplot$pLI >= 0.9) # 21
rownames(constrainedHC) <- 1:nrow(constrainedHC)
colnames(constrainedHC) <- c("Ensembl_gene_ID","Gene_symbol","pLI","dN/dS")
constrainedHC <- constrainedHC[c("Ensembl_gene_ID","Gene_symbol","dN/dS","pLI")]
smoothScatter(toplot$pLI, toplot$humchimpdndsr, main="Human vs Chimpanzee",xlab="pLI",ylab="dN/dS",ylim=(c(1,10)))
hist(toplot$pLI, main="Human vs Chimpanzee", xlab="pLI", breaks=seq(0,1,0.1))

# Human vs Orangutan
totalHO.pli <- merge(totalHO, pli, by="Ensembl_gene_ID")
sgho <- read.csv("sgho.csv",header=TRUE,stringsAsFactors = FALSE,sep=",")
dnds <- sgho[c("Ensembl_gene_ID","humorangdndsr")]
toplot <- merge(totalHO.pli,dnds,by="Ensembl_gene_ID") # 93
# Take out the most constrained genes (pLI >= 0.9)
constrainedHO <- subset(toplot, toplot$pLI >= 0.9) # 14
rownames(constrainedHO) <- 1:nrow(constrainedHO)
colnames(constrainedHO) <- c("Ensembl_gene_ID","Gene_symbol","pLI","dN/dS")
constrainedHO <- constrainedHO[c("Ensembl_gene_ID","Gene_symbol","dN/dS","pLI")]
smoothScatter(toplot$pLI, toplot$humorangdndsr, main="Human vs Orangutan",xlab="pLI",ylab="dN/dS",ylim=(c(1,10)))
hist(toplot$pLI, main="Human vs Orangutan", xlab="pLI", breaks=seq(0,1,0.1))

# Human vs Macaque
totalHM.pli <- merge(totalHM, pli, by="Ensembl_gene_ID")
sghm <- read.csv("sghm.csv",header=TRUE,stringsAsFactors = FALSE,sep=",")
dnds <- sghm[c("Ensembl_gene_ID","hummacdndsr")]
toplot <- merge(totalHM.pli,dnds,by="Ensembl_gene_ID") # 59
# Take out the most constrained genes (pLI >= 0.9)
constrainedHM <- subset(toplot, toplot$pLI >= 0.9) # 6
rownames(constrainedHM) <- 1:nrow(constrainedHM)
colnames(constrainedHM) <- c("Ensembl_gene_ID","Gene_symbol","pLI","dN/dS")
constrainedHM <- constrainedHM[c("Ensembl_gene_ID","Gene_symbol","dN/dS","pLI")]
smoothScatter(toplot$pLI, toplot$hummacdndsr, main="Human vs Macaque",xlab="pLI",ylab="dN/dS",ylim=(c(1,10)))
hist(toplot$pLI, main="Human vs Macaque", xlab="pLI", breaks=seq(0,1,0.1))

# Chimp vs Orangutan
totalCO.pli <- merge(totalCO, pli, by="Ensembl_gene_ID")
sgco <- read.csv("sgco.csv",header=TRUE,stringsAsFactors = FALSE,sep=",")
dnds <- sgco[c("Ensembl_gene_ID","chimporangdndsr")]
toplot <- merge(totalCO.pli,dnds,by="Ensembl_gene_ID") # 57 genes
# Take out the most constrained genes (pLI >= 0.9)
constrainedCO <- subset(toplot, toplot$pLI >= 0.9) # 9
rownames(constrainedCO) <- 1:nrow(constrainedCO)
colnames(constrainedCO) <- c("Ensembl_gene_ID","Gene_symbol","pLI","dN/dS")
constrainedCO <- constrainedCO[c("Ensembl_gene_ID","Gene_symbol","dN/dS","pLI")]
smoothScatter(toplot$pLI, toplot$chimporangdndsr, main="Chimpanzee vs Orangutan",xlab="pLI",ylab="dN/dS",ylim=(c(1,10)))
hist(toplot$pLI, main="Chimpanzee vs Orangutan", xlab="pLI", breaks=seq(0,1,0.1))

# Chimp vs Macaque
totalCM.pli <- merge(totalCM, pli, by="Ensembl_gene_ID")
sgcm <- read.csv("sgcm.csv",header=TRUE,stringsAsFactors = FALSE,sep=",")
dnds <- sgcm[c("Ensembl_gene_ID","chimpmacdndsr")]
toplot <- merge(totalCM.pli,dnds,by="Ensembl_gene_ID") # 44 genes
# Take out the most constrained genes (pLI >= 0.9)
constrainedCM <- subset(toplot, toplot$pLI >= 0.9) # 4
rownames(constrainedCM) <- 1:nrow(constrainedCM)
colnames(constrainedCM) <- c("Ensembl_gene_ID","Gene_symbol","pLI","dN/dS")
constrainedCM <- constrainedCM[c("Ensembl_gene_ID","Gene_symbol","dN/dS","pLI")]
smoothScatter(toplot$pLI, toplot$chimpmacdndsr, main="Chimpanzee vs Macaque",xlab="pLI",ylab="dN/dS",ylim=(c(1,10)))
hist(toplot$pLI, main="Chimpanzee vs Macaque", xlab="pLI", breaks=seq(0,1,0.1))

# Orangutan vs Macaque
totalOM.pli <- merge(totalOM, pli, by="Ensembl_gene_ID")
sgom <- read.csv("sgom.csv",header=TRUE,stringsAsFactors = FALSE,sep=",")
dnds <- sgom[c("Ensembl_gene_ID","orangmacdndsr")]
toplot <- merge(totalOM.pli,dnds,by="Ensembl_gene_ID") # 48 genes
# Take out the most constrained genes (pLI >= 0.9)
constrainedOM <- subset(toplot, toplot$pLI >= 0.9) # 6
rownames(constrainedOM) <- 1:nrow(constrainedOM)
colnames(constrainedOM) <- c("Ensembl_gene_ID","Gene_symbol","pLI","dN/dS")
constrainedOM <- constrainedOM[c("Ensembl_gene_ID","Gene_symbol","dN/dS","pLI")]
smoothScatter(toplot$pLI, toplot$orangmacdndsr, main="Orangutan vs Macaque",xlab="pLI",ylab="dN/dS",ylim=(c(1,10)))
hist(toplot$pLI, main="Orangutan vs Macaque", xlab="pLI", breaks=seq(0,1,0.1))

# Permutation test for high dN/dS and high constraint score genes
allexacgenes <- read.csv("exacgenes.csv", header = TRUE, sep =",")

rsample <- function(x){
  rs <- sample(allexacgenes$pLI, x, replace=F)
  highdnds <- subset(rs, rs >=0.9)
  return(length(highdnds))
}
pv <- function(a,b){
  test <- replicate(1000000, rsample(a))
  counts <- sum(test > b)
  p <- counts/1000000 
  return(p)
}
pvless <- function(a,b){
  test <- replicate(1000000, rsample(a))
  counts <- sum(test <= b)
  p <- counts/1000000 
  return(p)
}
# For PANTHER data
# Human vs Chimp
pvless(125,1) # < 1e0-06

# Human vs Orang
pvless(93,1) # < 1e0-06

# Human vs Macaca
pvless(59,0) # 1.2e-05

# Chimp vs Orang
pvless(57,1) # 0.000178

# Chimp vs Macaca
pvless(44,0) # 0.000195

# Orang vs Macaca
pvless(48,0) # 7.7e-05

# Genes common among all comparisons
allconstrained <- rbind(constrainedHC, constrainedHO, constrainedHM, constrainedCO, constrainedCM, constrainedOM)
repeats <- ddply(allconstrained,c("Ensembl_gene_ID","Gene_symbol"),nrow)
allcommon <- subset(repeats, repeats$V1=="6")

# Common among all human
allhuman <- rbind(constrainedHC, constrainedHM, constrainedHO)
repeats <- ddply(allhuman,c("Ensembl_gene_ID","Gene_symbol"),nrow)
allhuman <- subset(repeats, repeats$V1=="3")

# Common among all chimp
allchimp <- rbind(constrainedHC, constrainedCM, constrainedCO)
repeats <- ddply(allchimp,c("Ensembl_gene_ID","Gene_symbol"),nrow)
allchimp <- subset(repeats, repeats$V1=="3")

# Common among all orangutan
allorang <- rbind(constrainedHO, constrainedOM, constrainedCO)
repeats <- ddply(allorang,c("Ensembl_gene_ID","Gene_symbol"),nrow)
allorang <- subset(repeats, repeats$V1=="3")
rownames(allorang) <- 1:nrow(allorang)

# Common among all macaque
allmac <- rbind(constrainedHM, constrainedOM, constrainedCM)
repeats <- ddply(allmac,c("Ensembl_gene_ID","Gene_symbol"),nrow)
allmac <- subset(repeats, repeats$V1=="3")
rownames(allmac) <- 1:nrow(allmac)

########################################## dN/dS distribution plots ##########################################

allexacgenes <- read.csv("exacgenes.csv", header = TRUE, sep =",")
sghc <- read.csv("sghc.csv",header=T,stringsAsFactors = F,sep=",")
sgho <- read.csv("sgho.csv",header=T,stringsAsFactors = F,sep=",")
sghm <- read.csv("sghm.csv",header=T,stringsAsFactors = F,sep=",")
sgco <- read.csv("sgco.csv",header=T,stringsAsFactors = F,sep=",")
sgcm <- read.csv("sgcm.csv",header=T,stringsAsFactors = F,sep=",")
sgom <- read.csv("sgom.csv",header=T,stringsAsFactors = F,sep=",")

# Plots of dN/dS distributions
# HC
hcd <- density(sghc$humchimpdndsr)
plot(hcd,main="Distribution of dN/dS values",xlab="dN/dS",xlim=c(0,1),ylim=c(0,3), col="red")
msghc <- median(sghc$humchimpdndsr) # 0.2317666
# HO
hod <- density(sgho$humorangdndsr)
lines(hod, col="blue")
msgho <- median(sgho$humorangdndsr) # 0.2053318
# HM
hmd <- density(sghm$hummacdndsr)
lines(hmd, col="dark green")
msghm <- median(sghm$hummacdndsr) # 0.197201
# CO
cod <- density(sgco$chimporangdndsr)
lines(cod, col="orange")
msgco <- median(sgco$chimporangdndsr) # 0.2102625
# CM
cmd <- density(sgcm$chimpmacdndsr)
lines(cmd, col="black")
msgcm <- median(sgcm$chimpmacdndsr) # 0.202085
# OM
omd <- density(sgom$orangmacdndsr)
lines(omd, col="purple")
msgom <- median(sgom$orangmacdndsr) # 0.2066116

# log10(dN/dS)
# HC
sghcx <- sghc
sghcx$humchimpdndsr <- log10(sghcx$humchimpdndsr)
hcd <- density(sghcx$humchimpdndsr)
plot(hcd,main="Distribution of dN/dS values",xlab="log10(dN/dS)",ylim=c(0,1), col="red")
# HO
sghox <- sgho
sghox$humorangdndsr <- log10(sghox$humorangdndsr)
hod <- density(sghox$humorangdndsr)
lines(hod, col="blue")
# HM
sghmx <- sghm
sghmx$hummacdndsr <- log10(sghmx$hummacdndsr)
hmd <- density(sghmx$hummacdndsr)
lines(hmd, col="dark green")
# CO
sgcox <- sgco
sgcox$chimporangdndsr <- log10(sgcox$chimporangdndsr)
cod <- density(sgcox$chimporangdndsr)
lines(cod, col="orange")
# CM
sgcmx <- sgcm
sgcmx$chimpmacdndsr <- log10(sgcmx$chimpmacdndsr)
cmd <- density(sgcmx$chimpmacdndsr)
lines(cmd, col="black")
# OM
sgomx <- sgom
sgomx$orangmacdndsr <- log10(sgomx$orangmacdndsr)
omd <- density(sgomx$orangmacdndsr)
lines(omd, col="purple")

# Plot dN/dS against pLI
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]
# HC
sghc.pli <- merge(sghc, pli, by="Ensembl_gene_ID")
sghc.pli <- sghc.pli[c("Ensembl_gene_ID","humchimpdndsr","pLI")]
colnames(sghc.pli) <- c("Ensembl_gene_ID","dndsr","pLI")
smoothScatter(sghc.pli$pLI, sghc.pli$dndsr, main="Human vs Chimp", xlab="pLI",ylab="dN/dS")

sgho.pli <- merge(sgho, pli, by="Ensembl_gene_ID")
sgho.pli <- sgho.pli[c("Ensembl_gene_ID","humorangdndsr","pLI")]
colnames(sgho.pli) <- c("Ensembl_gene_ID","dndsr","pLI")
smoothScatter(sgho.pli$pLI, sgho.pli$dndsr, main="Human vs Orangutan", xlab="pLI",ylab="dN/dS")

sghm.pli <- merge(sghm, pli, by="Ensembl_gene_ID")
sghm.pli <- sghm.pli[c("Ensembl_gene_ID","hummacdndsr","pLI")]
colnames(sghm.pli) <- c("Ensembl_gene_ID","dndsr","pLI")
smoothScatter(sghm.pli$pLI, sghm.pli$dndsr, main="Human vs Macaque", xlab="pLI",ylab="dN/dS")

sgco.pli <- merge(sgco, pli, by="Ensembl_gene_ID")
sgco.pli <- sgco.pli[c("Ensembl_gene_ID","chimporangdndsr","pLI")]
colnames(sgco.pli) <- c("Ensembl_gene_ID","dndsr","pLI")
smoothScatter(sgco.pli$pLI, sgco.pli$dndsr, main="Chimp vs Orangutan", xlab="pLI",ylab="dN/dS")

sgcm.pli <- merge(sgcm, pli, by="Ensembl_gene_ID")
sgcm.pli <- sgcm.pli[c("Ensembl_gene_ID","chimpmacdndsr","pLI")]
colnames(sgcm.pli) <- c("Ensembl_gene_ID","dndsr","pLI")
smoothScatter(sgcm.pli$pLI, sgcm.pli$dndsr, main="Chimp vs Macaque", xlab="pLI",ylab="dN/dS")

sgom.pli <- merge(sgom, pli, by="Ensembl_gene_ID")
sgom.pli <- sgom.pli[c("Ensembl_gene_ID","orangmacdndsr","pLI")]
colnames(sgom.pli) <- c("Ensembl_gene_ID","dndsr","pLI")
smoothScatter(sgom.pli$pLI, sgom.pli$dndsr, main="Orangutan vs Macaque", xlab="pLI",ylab="dN/dS")

# Plot dN/dS distribution in deciles of pLI scores
pli <- allexacgenes[c("Ensembl_gene_ID","pLI")]

getsubsets <- function(x, y){
  d0.1 <- subset(x, x$pLI <= 0.1) # 9024
  d0.2 <- subset(x, x$pLI <= 0.2 & x$pLI > 0.1) # 676
  d0.3 <- subset(x, x$pLI <= 0.3 & x$pLI > 0.2) # 488
  d0.4 <- subset(x, x$pLI <= 0.4 & x$pLI > 0.3) # 453
  d0.5 <- subset(x, x$pLI <= 0.5 & x$pLI > 0.4) # 441
  d0.6 <- subset(x, x$pLI <= 0.6 & x$pLI > 0.5) # 390
  d0.7 <- subset(x, x$pLI <= 0.7 & x$pLI > 0.6) # 444
  d0.8 <- subset(x, x$pLI <= 0.8 & x$pLI > 0.7) # 459
  d0.9 <- subset(x, x$pLI <= 0.9 & x$pLI > 0.8) # 612
  d1.0 <- subset(x, x$pLI <= 1.0 & x$pLI > 0.9) # 2768
  d0.1plot <- density(d0.1$dndsr)
  d0.2plot <- density(d0.2$dndsr)
  d0.3plot <- density(d0.3$dndsr)
  d0.4plot <- density(d0.4$dndsr)
  d0.5plot <- density(d0.5$dndsr)
  d0.6plot <- density(d0.6$dndsr)
  d0.7plot <- density(d0.7$dndsr)
  d0.8plot <- density(d0.8$dndsr)
  d0.9plot <- density(d0.9$dndsr)
  d1.0plot <- density(d1.0$dndsr)
  plot(d0.1plot, main=y, xlab="log10(dN/dS)",ylim=c(0,1.2), col="red")
  lines(d0.2plot, col="orange")
  lines(d0.3plot, col="yellow")
  lines(d0.4plot, col="green")
  lines(d0.5plot, col="dark green")
  lines(d0.6plot, col="cyan")
  lines(d0.7plot, col="blue")
  lines(d0.8plot, col="dark blue")
  lines(d0.9plot, col="purple")
  lines(d1.0plot, col="black")
}

# Plot each decile per graph
getsubsetsdecile <- function(a, b, c, d, e, f, firstnum, secondnum, k, q, w, xname, xlim1, xlim2){
  da <- subset(a, a$pLI <= firstnum & a$pLI > secondnum)
  db <- subset(b, b$pLI <= firstnum & b$pLI > secondnum)
  dc <- subset(c, c$pLI <= firstnum & c$pLI > secondnum)
  dd <- subset(d, d$pLI <= firstnum & d$pLI > secondnum)
  de <- subset(e, e$pLI <= firstnum & e$pLI > secondnum)
  df <- subset(f, f$pLI <= firstnum & f$pLI > secondnum)
  daplot <- density(da$dndsr)
  dbplot <- density(db$dndsr)
  dcplot <- density(dc$dndsr)
  ddplot <- density(dd$dndsr)
  deplot <- density(de$dndsr)
  dfplot <- density(df$dndsr)
  plot(daplot, main=k, xlab=xname,ylim=c(q,w), xlim=c(xlim1,xlim2), col="red")
  lines(dbplot, col="blue")
  lines(dcplot, col="green")
  lines(ddplot, col="orange")
  lines(deplot, col="black")
  lines(dfplot, col="purple")
}

sghc.pli <- merge(sghc, pli, by="Ensembl_gene_ID")
sghc.pli <- sghc.pli[c("Ensembl_gene_ID","humchimpdndsr","pLI")]
colnames(sghc.pli) <- c("Ensembl_gene_ID","dndsr","pLI")
sghc.pli$dndsr <- log10(sghc.pli$dndsr)
getsubsets(sghc.pli, "Human vs Chimp")

sgho.pli <- merge(sgho, pli, by="Ensembl_gene_ID")
sgho.pli <- sgho.pli[c("Ensembl_gene_ID","humorangdndsr","pLI")]
colnames(sgho.pli) <- c("Ensembl_gene_ID","dndsr","pLI")
sgho.pli$dndsr <- log10(sgho.pli$dndsr)
getsubsets(sgho.pli, "Human vs Orangutan")

sghm.pli <- merge(sghm, pli, by="Ensembl_gene_ID")
sghm.pli <- sghm.pli[c("Ensembl_gene_ID","hummacdndsr","pLI")]
colnames(sghm.pli) <- c("Ensembl_gene_ID","dndsr","pLI")
sghm.pli$dndsr <- log10(sghm.pli$dndsr)
getsubsets(sghm.pli, "Human vs Macaque")

sgco.pli <- merge(sgco, pli, by="Ensembl_gene_ID")
sgco.pli <- sgco.pli[c("Ensembl_gene_ID","chimporangdndsr","pLI")]
colnames(sgco.pli) <- c("Ensembl_gene_ID","dndsr","pLI")
sgco.pli$dndsr <- log10(sgco.pli$dndsr)
getsubsets(sgco.pli, "Chimpanzee vs Orangutan")

sgcm.pli <- merge(sgcm, pli, by="Ensembl_gene_ID")
sgcm.pli <- sgcm.pli[c("Ensembl_gene_ID","chimpmacdndsr","pLI")]
colnames(sgcm.pli) <- c("Ensembl_gene_ID","dndsr","pLI")
sgcm.pli$dndsr <- log10(sgcm.pli$dndsr)
getsubsets(sgcm.pli, "Chimp vs Macaque")

sgom.pli <- merge(sgom, pli, by="Ensembl_gene_ID")
sgom.pli <- sgom.pli[c("Ensembl_gene_ID","orangmacdndsr","pLI")]
colnames(sgom.pli) <- c("Ensembl_gene_ID","dndsr","pLI")
sgom.pli$dndsr <- log10(sgom.pli$dndsr)
getsubsets(sgom.pli, "Orangutan vs Macaque")

getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.1,0,"0-0.1 pLI score", 0, 1.2, "log10(dN/dS)", -4, 1)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.2,0.1,"0.1-0.2 pLI score", 0, 1.2, "log10(dN/dS)", -4, 1)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.3,0.2,"0.2-0.3 pLI score", 0, 1.2, "log10(dN/dS)", -4, 1)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.4,0.3,"0.3-0.4 pLI score", 0, 1.2, "log10(dN/dS)", -4, 1)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.5,0.4,"0.4-0.5 pLI score", 0, 1.2, "log10(dN/dS)", -4, 1)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.6,0.5,"0.5-0.6 pLI score", 0, 1.2, "log10(dN/dS)", -4, 1)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.7,0.6,"0.6-0.7 pLI score", 0, 1.2, "log10(dN/dS)", -4, 1)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.8,0.7,"0.7-0.8 pLI score", 0, 1.2, "log10(dN/dS)", -4, 1)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.9,0.8,"0.8-0.9 pLI score", 0, 1.2, "log10(dN/dS)", -4, 1)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,1.0,0.9,"0.9-1.0 pLI score", 0, 1.2, "log10(dN/dS)", -4, 1)

# No log10
sghc.pli <- merge(sghc, pli, by="Ensembl_gene_ID")
sghc.pli <- sghc.pli[c("Ensembl_gene_ID","humchimpdndsr","pLI")]
colnames(sghc.pli) <- c("Ensembl_gene_ID","dndsr","pLI") # SD = 0.4730302

sgho.pli <- merge(sgho, pli, by="Ensembl_gene_ID")
sgho.pli <- sgho.pli[c("Ensembl_gene_ID","humorangdndsr","pLI")]
colnames(sgho.pli) <- c("Ensembl_gene_ID","dndsr","pLI") # SD = 0.3478538

sghm.pli <- merge(sghm, pli, by="Ensembl_gene_ID")
sghm.pli <- sghm.pli[c("Ensembl_gene_ID","hummacdndsr","pLI")]
colnames(sghm.pli) <- c("Ensembl_gene_ID","dndsr","pLI") # SD = 0.2844124

sgco.pli <- merge(sgco, pli, by="Ensembl_gene_ID")
sgco.pli <- sgco.pli[c("Ensembl_gene_ID","chimporangdndsr","pLI")]
colnames(sgco.pli) <- c("Ensembl_gene_ID","dndsr","pLI") # SD = 0.3244207

sgcm.pli <- merge(sgcm, pli, by="Ensembl_gene_ID")
sgcm.pli <- sgcm.pli[c("Ensembl_gene_ID","chimpmacdndsr","pLI")]
colnames(sgcm.pli) <- c("Ensembl_gene_ID","dndsr","pLI") # SD = 0.2869026

sgom.pli <- merge(sgom, pli, by="Ensembl_gene_ID")
sgom.pli <- sgom.pli[c("Ensembl_gene_ID","orangmacdndsr","pLI")]
colnames(sgom.pli) <- c("Ensembl_gene_ID","dndsr","pLI") # SD = 0.271154

getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.1,0,"0-0.1 pLI score", 0, 3, "dN/dS", 0, 1.2)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.2,0.1,"0.1-0.2 pLI score", 0, 3, "dN/dS", 0, 1.2)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.3,0.2,"0.2-0.3 pLI score", 0, 3, "dN/dS", 0, 1.2)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.4,0.3,"0.3-0.4 pLI score", 0, 3, "dN/dS", 0, 1.2)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.5,0.4,"0.4-0.5 pLI score", 0, 3, "dN/dS", 0, 1.2)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.6,0.5,"0.5-0.6 pLI score", 0, 3, "dN/dS", 0, 1.2)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.7,0.6,"0.6-0.7 pLI score", 0, 3, "dN/dS", 0, 1.2)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.8,0.7,"0.7-0.8 pLI score", 0, 3, "dN/dS", 0, 1.2)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,0.9,0.8,"0.8-0.9 pLI score", 0, 3, "dN/dS", 0, 1.2)
getsubsetsdecile(sghc.pli,sgho.pli,sghm.pli,sgco.pli,sgcm.pli,sgom.pli,1.0,0.9,"0.9-1.0 pLI score", 0, 3, "dN/dS", 0, 1.2)

# Group comparisons by species
human <- data.frame(rbind(sghc.pli,sgho.pli,sghm.pli))
chimp <- data.frame(rbind(sghc.pli,sgco.pli,sgcm.pli))
orang <- data.frame(rbind(sgho.pli,sgco.pli,sgom.pli))
mac <- data.frame(rbind(sghm.pli,sgcm.pli,sgom.pli))
# Plot
dhuman <- density(human$dndsr)
dchimp <- density(chimp$dndsr)
dorang <- density(orang$dndsr)
dmac <- density(mac$dndsr)

plot(dhuman, col="red", main="Distribution of dN/dS by Species Groups", xlab="dN/dS", xlim=c(0,1.2),ylim=c(0,3))
lines(dchimp, col="blue")
lines(dorang, col="orange")
lines(dmac,col="dark green")

# log10 dN/dS
humanx <- human
chimpx <- chimp
orangx <- orang
macx <- mac

humanx$dndsr <- log10(humanx$dndsr)
chimpx$dndsr <- log10(chimpx$dndsr)
orangx$dndsr <- log10(orangx$dndsr)
macx$dndsr <- log10(macx$dndsr)

dhumanx <- density(humanx$dndsr)
dchimpx <- density(chimpx$dndsr)
dorangx <- density(orangx$dndsr)
dmacx <- density(macx$dndsr)

plot(dhumanx, col="red", main="Distribution of dN/dS by Species Groups", xlab="log10(dN/dS)", ylim=c(0,1))
lines(dchimpx, col="blue")
lines(dorangx, col="orange")
lines(dmacx,col="dark green")

############################################## PAML and ggtree ###############################################
# For PAML and ggtree, testing with TAS2R14

# Load libraries
marthuman = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
martchimp = useEnsembl(biomart="ensembl", dataset="ptroglodytes_gene_ensembl", GRCh=37)
martorang = useEnsembl(biomart="ensembl", dataset="pabelii_gene_ensembl", GRCh=37)
martmac = useEnsembl(biomart="ensembl", dataset="mmulatta_gene_ensembl", GRCh=37)

# Human
Hgene <- getSequence(id="ENSG00000212127", type="ensembl_gene_id", seqType="coding", mart=marthuman)
# Keep the longest sequence if repeated
Hgene <- cbind(nchar(Hgene$coding), Hgene)
colnames(Hgene) <- c("Seq_length", "Coding_seq", "Ensembl_gene_ID")
Hgene <- Hgene[order(Hgene$Seq_length, decreasing=TRUE),]
Hgene <- subset(Hgene, duplicated(Hgene$Ensembl_gene_ID)==FALSE)
rownames(Hgene) <- 1:nrow(Hgene)

# Chimp
Cgene <- getSequence(id="ENSPTRG00000032538", type="ensembl_gene_id", seqType="coding", mart=martchimp)
# Keep the longest sequence if repeated
Cgene <- cbind(nchar(Cgene$coding), Cgene)
colnames(Cgene) <- c("Seq_length", "Coding_seq", "Ensembl_gene_ID")
Cgene <- Cgene[order(Cgene$Seq_length, decreasing=TRUE),]
Cgene <- subset(Cgene, duplicated(Cgene$Ensembl_gene_ID)==FALSE)
rownames(Cgene) <- 1:nrow(Cgene)

# Orangutan
Ogene <- getSequence(id="ENSPPYG00000004289", type="ensembl_gene_id", seqType="coding", mart=martorang)
# Keep the longest sequence if repeated
Ogene <- cbind(nchar(Ogene$coding), Ogene)
colnames(Ogene) <- c("Seq_length", "Coding_seq", "Ensembl_gene_ID")
Ogene <- Ogene[order(Ogene$Seq_length, decreasing=TRUE),]
Ogene <- subset(Ogene, duplicated(Ogene$Ensembl_gene_ID)==FALSE)
rownames(Ogene) <- 1:nrow(Ogene)

# Macaca
Mgene <- getSequence(id="ENSMMUG00000019281", type="ensembl_gene_id", seqType="coding", mart=martmac)
# Keep the longest sequence if repeated
Mgene <- cbind(nchar(Mgene$coding), Mgene)
colnames(Mgene) <- c("Seq_length", "Coding_seq", "Ensembl_gene_ID")
Mgene <- Mgene[order(Mgene$Seq_length, decreasing=TRUE),]
Mgene <- subset(Mgene, duplicated(Mgene$Ensembl_gene_ID)==FALSE)
rownames(Mgene) <- 1:nrow(Mgene)

# Assemble sequences into FASTA format for alignment
fasfile <-  character(length = 8) # 1 gene * 4 species + 4 gene names
fasfile[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)] = sprintf("> %sH", "TAS2R14")
fasfile[c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)] = Hgene$Coding_seq
fasfile[c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)] = sprintf("> %sC", "TAS2R14")
fasfile[c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)] = Cgene$Coding_seq
fasfile[c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)] = sprintf("> %sO", "TAS2R14")
fasfile[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)] = Ogene$Coding_seq
fasfile[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)] = sprintf("> %sM", "TAS2R14")
fasfile[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)] = Mgene$Coding_seq
writeLines(fasfile, "tas2r14.fasta")

# ggtree
setwd("~/ExAC_dataset_project/csv_records/PAML")
library(ggtree)
# load mlc file (Output from CodeML. Copy mlc file into R working directory.)
mlc <- read.codeml_mlc("mlc")
p1 <- plot(mlc, annotation="dN_vs_dS", annotation.color="blue", ndigits=3) + xlim(0, 0.4)
p1
p2 <- plot(mlc, branch.length="dN_vs_dS", annotation="dN_vs_dS", annotation.color="blue", ndigits=3) + xlim(0, 12)
p2

# For the mlc file that used a tree with macaque as outgroup
mlc <- read.codeml_mlc("mlcout")
p1 <- plot(mlc, annotation="dN_vs_dS", annotation.color="blue", ndigits=3) + xlim(0, 0.25)
p1
p2 <- plot(mlc, branch.length="dN_vs_dS", annotation="dN_vs_dS", annotation.color="blue", ndigits=3) + xlim(0, 1000)
p2

# Using non-primate as an outgroup
# Add on kangeroo rat (Dipodomys ordii) to the alignment file
martkang = useEnsembl(biomart="ensembl", dataset="dordii_gene_ensembl", GRCh=37)
# Kangeroo rat
Kgene <- getSequence(id="ENSDORG00000004548", type="ensembl_gene_id", seqType="coding", mart=martkang)
# Keep the longest sequence if repeated
Kgene <- cbind(nchar(Kgene$coding), Kgene)
colnames(Kgene) <- c("Seq_length", "Coding_seq", "Ensembl_gene_ID")
Kgene <- Kgene[order(Kgene$Seq_length, decreasing=TRUE),]
Kgene <- subset(Kgene, duplicated(Kgene$Ensembl_gene_ID)==FALSE)
rownames(Kgene) <- 1:nrow(Kgene)

# Assemble sequences into FASTA format for alignment
fasfile <-  character(length = 10) # 1 gene * 5 species + 5 gene names
fasfile[c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)] = sprintf("> %sH", "TAS2R14")
fasfile[c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)] = Hgene$Coding_seq
fasfile[c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)] = sprintf("> %sC", "TAS2R14")
fasfile[c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)] = Cgene$Coding_seq
fasfile[c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)] = sprintf("> %sO", "TAS2R14")
fasfile[c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)] = Ogene$Coding_seq
fasfile[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)] = sprintf("> %sM", "TAS2R14")
fasfile[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE)] = Mgene$Coding_seq
fasfile[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)] = sprintf("> %sK", "TAS2R14")
fasfile[c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)] = Mgene$Coding_seq
writeLines(fasfile, "tas2r14k.fasta")

# ggplot for using kangeroo rat as outgroup
setwd("~/ExAC_dataset_project/csv_records/PAML")
# load mlc file
mlc <- read.codeml_mlc("mlckang")
p2 <- plot(mlc, branch.length="dN_vs_dS", annotation="dN_vs_dS", annotation.color="blue", ndigits=3) + xlim(0, 12)
p2