# QIIME2 to Phyloseq

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

setwd("~/Documents/3_Holobiont_Stability/Stats/QIIME2")
### load data ####
# load OTU table (family level, OTUs in columns, samples in columns)
OTU_table<-read.csv('otu_table5976.csv', header = TRUE, sep = ",", dec = ".", strip.white = TRUE)
head(OTU_table)
str(OTU_table)

OTU_table_T<-read.csv('otu_table5976_T.csv', header = TRUE, sep = ",", dec = ".", strip.white = TRUE)
head(OTU_table_T) # OTUs in rows
str(OTU_table_T)

METADATA<-read.csv("Metadata_2.csv", header = TRUE, sep = ",")
head(METADATA)
levels(METADATA$Host)

TAXADATA<-read.csv("taxonomy.csv", header = TRUE, sep = ",")
head(TAXADATA)

# join datasets
library(dplyr)
library(tidyr)
OTU_table_new<-right_join(METADATA, OTU_table, by="Sample_ID")
levels(OTU_table_new$Host)

### prep Phyloseq Table ####
# Taxa table
TAXA_TABLE<-as.data.frame(TAXADATA[,1])
rownames(TAXADATA)<-TAXADATA[,1]
TAXA_TABLE_split<-TAXADATA %>% separate(Taxon, c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep=";", remove=FALSE) 
TAXA_TABLE_matrix<-as.matrix(TAXA_TABLE_split)
TAXA_TABLE<-tax_table(TAXA_TABLE_matrix) # tax_table=matrix

# OTU table
rownames(OTU_table_T)<-OTU_table_T[,1]
OTU_table_T<-data.frame(OTU_table_T[,-1])

OTU_TABLE<-otu_table(OTU_table_T, taxa_are_rows = TRUE)
sample_names(OTU_TABLE)

# Sample data
SAMPLEDATA<-OTU_table_new[,1:14]
rownames(SAMPLEDATA)<-OTU_table_new$Sample_ID
SAMPLEDATA<-data.frame(SAMPLEDATA[,-1])
levels(SAMPLEDATA$Host)

SAMPLEDATA<-sample_data(SAMPLEDATA) #Sample_data = data.frame
sample_names(SAMPLEDATA)
SAMPLEDATA$SamplingTimepoint<-factor(SAMPLEDATA$SamplingTimepoint, levels=c("before disturbance", "24h after disturbance", "168h after disturbance"))
SAMPLEDATA$MicrobiomeDiversity<-factor(SAMPLEDATA$MicrobiomeDiversity, levels=c("VLMD","LMD","HMD","VHMD"))

#combine Sampledata, Taxadata, OTU table to phyloseq table
PHYLOSEQ_TABLE<-phyloseq(OTU_TABLE, SAMPLEDATA, TAXA_TABLE)
PHYLOSEQ_TABLE

save(PHYLOSEQ_TABLE, file="PHYLOSEQ_TABLE.RData")
