###############
# Description #
###############
#* The purpose of this script is to fine map the APOE locus on Chromosome 19 and
#* highlight the distribution of the most likely candidate SNPs across the locus.
#* 
#* Data source: GWAS results taken from Marioni et al., 2018
#* The most likely candidates SNPs are based from echolactor software 
#* 
#* 
#* Steps
#* 1. Import GWAS dataset
#* 2. Isolate DPOE locus
#* 3. Create Manhatten plots
#* 4. Create QQplots for each chromosome
#* 4. Integrate GWAS with RegulomedB variant information
#* 5. Create Plots stratified by SNP pathogenic rank
#* 6. Highlight likely candidate variants

####################
# Install packages #
####################

#Install packages
library(tidyverse)

############################
# 1. Download GWAS dataset #
############################

# GWAS data from Marioni et al., 2018
df = read.table("../Raw_data/3_UKB_AD_parental_meta_summary_output_June2019.txt", header = T)

# Modify P-values that are 0 to prevent -log(P) = inf
df$P[abs(df$P)<1e-999]=1e-999

#########################
# 2. Isolate DPOE locus #
#########################

#get DPOE loci
APOE_GWAS = df %>% filter(CHR == 19)
#order by BP position
arrange(APOE_GWAS, BP)
#filter DPOE region --> loci determined by top SNPs from Echolocator
APOE_GWAS2 <- APOE_GWAS %>% filter(between(BP, 45411941-500000, 45411941+500000))

###############
# 3. Plot QTL #
###############

# plot QTL region
ggplot(APOE_GWAS2, aes(x = BP, y = -log(P))) + 
  geom_point()

## additional higlight top SNPs on plot
#* Top SNPs identified through EchoLocator  
target = c('rs11556505', "rs12721046", "rs12721051", "rs12972156", "rs142042446", "rs150966173", "rs157592", "rs41289512")
highlight_df = filter(APOE_GWAS2, SNP %in% target)

# plot QTL region with highlighted points
ggplot(APOE_GWAS2, aes(x = BP, y = -log(P))) + 
  geom_point() +
  geom_point(data=highlight_df, aes(x = BP, y = -log(P), color=SNP)) +
  geom_vline(xintercept = 45351516, linetype="dotdash", color = "darkgrey", size=1) +
  geom_vline(xintercept = 45424514, linetype="dotdash", color = "darkgrey", size=1)

