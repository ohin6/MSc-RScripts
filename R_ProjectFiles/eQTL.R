###############
# Description #
###############
#* This script fine maps APOE locus by co locating GWAS with eQTL data. This uses
#* script allows the user to select expression analysis data from a catalog of
#* eQTL taken from different tissue types which may be more relevant to the
#* biological pathway of the disease.
#* 
#* This uses Tabix function as an efficient method for extracting relevant 
#* information
#* 
#* 
#* Steps
#* 1. Set up link with online eQTL catalog
#* 2. Define QTL region
#* 3. Extract eQTL for different tissue types and plot

####################
# Install packages #
####################

require(tidyverse)
library("readr")
library("coloc")
library("seqminer")

###########################################
# 1. Set up link with online eQTL catalog #
###########################################

# Set link to online eQTL catalog
tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

# Get list of eQTL data from GTex study 
imported_tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths_imported.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

# Create function for extracting information from study using Tabix
import_eQTLCatalogue <- function(ftp_path, region, selected_gene_id, column_names, verbose = TRUE){
  
  if(verbose){
    print(ftp_path)
  }
  
  #Fetch summary statistics with seqminer
  fetch_table = seqminer::tabix.read.table(tabixFile = ftp_path, tabixRange = region, stringsAsFactors = FALSE) %>%
    dplyr::as_tibble()
  colnames(fetch_table) = column_names
  
  #Remove rsid duplicates and multi-allelic variant
  summary_stats = dplyr::filter(fetch_table, gene_id == selected_gene_id) %>%
    dplyr::select(-rsid) %>% 
    dplyr::distinct() %>% #rsid duplicates
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1) #Multialllics
}

########################
# 2. Define QTL region #
########################
#* Region of interest is the APOE QTL on chromosome 19 located at 4490864bp variant 
#* rs429538. The GTEx uses the GRCh38 build which again co localisation will be 
#* chrm19 4490864bp. The QTL region of interest will be 500kb either side of peak
#* SNP.
#* 

#define region of interest (APOE)
#=/- 500kb of rs429538 
region = "19:44408684-45408684"

#identify studies in GTEx
APOE_Frontal_cortex_df = dplyr::filter(tabix_paths, study == "GTEx")
#view tissue types 
unique(APOE_Frontal_cortex_df$tissue_label)


#######################################################
# 3. Extract eQTL for different tissue types and plot #
#######################################################

#1) Frontal cortex
APOE_Frontal_cortex_df = dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "brain (DLPFC)")
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE_Frontal_cortex_df$ftp_path, n_max = 1))
#Import summary statistics
summary_stats1 = import_eQTLCatalogue(APOE_Frontal_cortex_df$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats1, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "brain frontal cortex") +
  geom_point()+
  ylim(0,5.5)

#2) brain (anterior cingulate cortex)
APOE_brain_cingulate = dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "brain (anterior cingulate cortex)")
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE_brain_cingulate$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE_brain_cingulate$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "brain (anterior cingulate cortex)") +
  geom_point()+
  ylim(0,5.5)

#3) brain cerebellum
APOE_cerebellum = dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "brain (cerebellum)")
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE_cerebellum$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE_cerebellum$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "brain cerebellum") +
  geom_point()+
  ylim(0,5.5)

#4) brain hypothalamus
APOE_hypothalamus= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "brain (hypothalamus)")
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE_hypothalamus$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE_hypothalamus$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "brain (hypothalamus)") +
  geom_point()+
  ylim(0,5.5)

#5) "brain (putamen)" 
APOE= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "brain (putamen)")
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "brain (putamen)") +
  geom_point()+
  ylim(0,5.5)


#6) "brain (substantia nigra)"  
APOE= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "brain (substantia nigra)" )
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "brain (substantia nigra)" ) +
  geom_point()+
  ylim(0,5.5)


#7) brain (amygdala)"  
APOE= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "brain (amygdala)" )
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "brain (amygdala)") +
  geom_point()+
  ylim(0,5.5)

#8) "brain (caudate)"  
APOE= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "brain (caudate)" )
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "brain (caudate)") +
  geom_point()+
  ylim(0,5.5)


#9) "brain (cortex)"  
APOE= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "brain (cortex)")
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "brain (cortex)") +
  geom_point()+
  ylim(0,5.5)

#10) "brain (hippocampus)"  
APOE= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "brain (hippocampus)")
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "brain (hippocampus)") +
  geom_point()+
  ylim(0,5.5)

#11) "brain (nucleus accumbens)"  
APOE= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "brain (nucleus accumbens)")
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "brain (nucleus accumbens)") +
  geom_point()+
  ylim(0,5.5)

#12) ""brain (spinal cord)""  
APOE= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "brain (spinal cord)")
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "brain (spinal cord)") +
  geom_point()+
  ylim(0,5.5)


#13) LCL  
APOE= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "LCL")
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "LCL") +
  geom_point() +
  ylim(0,5.5)


#14) Blood  
APOE= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "blood")
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "blood") +
  geom_point() +
  ylim(0,5.5)

#15) "adipose (visceral)"  
APOE= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "adipose (visceral)")
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "adipose (visceral)") +
  geom_point() +
  ylim(0,5.5)

#16) "adipose"  
APOE= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "adipose")
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "adipose") +
  geom_point() +
  ylim(0,5.5)

#16) "kidney (cortex)"  
APOE= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "kidney (cortex)" )
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "kidney (cortex)" ) +
  geom_point() +
  ylim(0,5.5)

#18) "liver"  
APOE= dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "liver" )
#Extract column names from first file
column_names = colnames(readr::read_tsv(APOE$ftp_path, n_max = 1))
#Import summary statistics
summary_stats2 = import_eQTLCatalogue(APOE$ftp_path, region, selected_gene_id = "ENSG00000130203", column_names)
ggplot(summary_stats2, aes(x = position, y = -log(pvalue, 10))) +
  labs(title = "liver" ) +
  geom_point() +
  ylim(0,5.5)
