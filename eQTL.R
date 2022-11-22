###############
# Description #
###############
#* The purpose of this script is to fine map a QTL by integrating eQTL with GWAS
#* data. Candidate SNPs are identified by co localisation of GWAS and eQTL data.
#* 
#* This script allows you extract eQTL information from different cell types and
#* integrate it with GWAS at a specific region.
#* 
#* This Script involves
#* 
#* 1. Set online directory with eQTL catalog
#* 2. Define QTL region
#* 3. Extract eQTL

####################
# Install packages #
####################

require(tidyverse)
require(readr)
require(coloc)
require(GenomicRanges)
require(seqminer)


#############################################
# 1. Set online directory with eQTL catalog #
#############################################

# Create link to online eQTL catalogs
tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()
imported_tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths_imported.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% dplyr::as_tibble()

# Create function to extract eQTL information
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

# Define region of interest which is the APOE QTL on Chrm19
#* GTEX is in GRCh38 build --> focus on rs429538 which is APOE E4 variant which is at 44908684bp
#* =/- 500kb of rs429538

region = "19:44408684-45408684"

###################
# 3. Extract eQTL #
###################

# Identify studies in GTEx (gene expression)
APOE_Frontal_cortex_df = dplyr::filter(tabix_paths, study == "GTEx")

# view tissue types in GTEx studies
unique(APOE_Frontal_cortex_df$tissue_label)


# Get eQTL for differnt tissue types

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



#1) Frontal cortex
APOE_Frontal_cortex_df = dplyr::filter(tabix_paths, study == "GTEx", tissue_label == "brain (DLPFC)", quant_method == 'exon')


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
