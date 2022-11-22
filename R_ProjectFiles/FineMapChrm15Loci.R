###############
# Description #
###############
#* The purpose of this is to fine map a novel LOAD locus identified on Chrm15 
#* which was identified, from Marioni et al., 2019 study in order to identify
#* potential candidate genes.
#* 
#* This involves
#* 
#* 1. Import GWAS dataset
#* 2. Filter loci on chrm15
#* 3. Plot loci in large
#* 4. Integrate GWAS with RegulomedB variant information
#* 5. Create Plots stratified by SNP pathogenic rank
#* 6. Highlight likely candidate variants

####################
# Install packages #
####################
install.packages("qqman")
require(qqman)
require(tidyverse)
require(ggpubr)
require(Gviz)
library(GenomicRanges)


############################
# 1. Downlaod GWAS dataset #
############################

# GWAS data from Marioni et al., 2018
df = read.table("../Raw_data/3_UKB_AD_parental_meta_summary_output_June2019.txt", header = T)

############################
# 2. Filter loci on chrm15 #
############################

# filter loci on chromosome 15
df15 = df %>% filter(CHR == 15 & between(BP,58000000,60000000))

################
# 3. Plot loci #
################

# plot region
ggplot(df15, aes(x = BP, y = -log(P))) + 
  geom_point() +
  geom_hline(yintercept=8, color = "blue") +
  geom_hline(yintercept=5, color = "red") +
  # highlight loci --> from table
  geom_vline(xintercept = 58873555, linetype="longdash", color = "dark grey", size=2) +
  geom_vline(xintercept = 59120077, linetype="longdash", color = "dark grey", size=2) +
  theme_bw()

#########################################################
# 4. Integrate GWAS with RegulomedB variant information #
#########################################################

# Get df of SNPs within reported loci
df15b = df %>% filter(CHR == 15 & between(BP,58873555,59120077))

# Import Regulomedb data of isolated region 
regAll = as_tibble(read_csv("../Raw_data/regulomedb_results_chr15.csv"))

# Create new column containing a simplified regulomedB pathogenic ranking groups
# 1 to 7
regAll$rank = 'na'
regAll = regAll %>%
  mutate(rank = if_else(ranking %in% str_subset(unique(regAll$ranking), '^1'), '1', rank)) %>%
  mutate(rank = if_else(ranking %in% str_subset(unique(regAll$ranking), '^2'), '2', rank)) %>%
  mutate(rank = if_else(ranking %in% str_subset(unique(regAll$ranking), '^3'), '3', rank)) %>%
  mutate(rank = if_else(ranking == '4', '4', rank)) %>%
  mutate(rank = if_else(ranking == '5', '5', rank)) %>%
  mutate(rank = if_else(ranking == '6', '6', rank)) %>%
  mutate(rank = if_else(ranking == '7', '7', rank))

# Change row variant id name to match both datasets (allow left join)
regAll = regAll %>%
  rename(SNP = rsids)

#left join regulomeDb info with GWAS data
regAll = left_join(regAll, df15b, by = "SNP")

# Subset data based on their regulomedB pathogenic ranking (1-7)
highlight_df1 = filter(regAll, rank %in% 1)
highlight_df2 = filter(regAll, rank %in% 2)
highlight_df3 = filter(regAll, rank %in% 3)
highlight_df4 = filter(regAll, rank %in% 4)
highlight_df5 = filter(regAll, rank %in% 5)
highlight_df6 = filter(regAll, rank %in% 6)
highlight_df7 = filter(regAll, rank %in% 7)

#####################################################
# 5. Create Plots stratified by SNP pathogenic rank #
#####################################################

#plot all SNPs
plt_all = ggplot(df15, aes(x = BP, y = -log(P))) + 
  geom_point() +
  geom_point(data = regAll, aes(x=BP,y=-log(P), colour=rank)) +
  scale_colour_manual(values = c("red", "orange", "yellow"," green", "blue", "darkblue", "purple")) +
  geom_hline(yintercept=8, color = "blue") +
  geom_hline(yintercept=5, color = "red") +
  #highlight loci -->from table
  geom_vline(xintercept = 58873555, linetype="longdash", color = "dark grey", size=1) +
  geom_vline(xintercept = 59120077, linetype="longdash", color = "dark grey", size=1) +
  theme(axis.text.x = element_text(size=10, angle=45))


# Stratify SNP by ranks 1-2
plt_1 = ggplot(df15, aes(x = BP, y = -log(P))) + 
  geom_point() +
  geom_point(data=highlight_df1, aes(x=BP,y=-log(P)), color='red') +
  geom_point(data=highlight_df2, aes(x=BP,y=-log(P)), color='orange') +
  geom_hline(yintercept=8, color = "blue") +
  geom_hline(yintercept=5, color = "red") +
  #highlight loci -->from table
  geom_vline(xintercept = 58873555, linetype="longdash", color = "dark grey", size=1) +
  geom_vline(xintercept = 59120077, linetype="longdash", color = "dark grey", size=1) +
  theme(axis.text.x = element_text(size=10, angle=45))


# Stratify SNP by ranks 3-5
plt_2= ggplot(df15, aes(x = BP, y = -log(P))) + 
  geom_point() +
  geom_point(data=highlight_df5, aes(x=BP,y=-log(P)), color='blue') +
  geom_point(data=highlight_df4, aes(x=BP,y=-log(P)), color='green') +
  geom_point(data=highlight_df3, aes(x=BP,y=-log(P)), color='yellow') +
  geom_hline(yintercept=8, color = "blue") +
  geom_hline(yintercept=5, color = "red") +
  #highlight loci -->from table
  geom_vline(xintercept = 58873555, linetype="longdash", color = "dark grey", size=1) +
  geom_vline(xintercept = 59120077, linetype="longdash", color = "dark grey", size=1) +
  theme(axis.text.x = element_text(size=10, angle=45))

# Stratify SNP by ranks 6 and 7
plt_3=ggplot(df15, aes(x = BP, y = -log(P))) + 
  geom_point() +
  geom_point(data=highlight_df7, aes(x=BP,y=-log(P)), color='purple') +
  geom_point(data=highlight_df6, aes(x=BP,y=-log(P)), color='darkblue') +
  geom_hline(yintercept=8, color = "blue") +
  geom_hline(yintercept=5, color = "red") +
  #highlight loci -->from table
  geom_vline(xintercept = 58873555, linetype="longdash", color = "dark grey", size=1) +
  geom_vline(xintercept = 59120077, linetype="longdash", color = "dark grey", size=1) +
  theme(axis.text.x = element_text(size=10, angle=45))

#Join plots
ggarrange(plt_all, plt_1, plt_2, plt_3, labels = c('a)', 'b)', 'c)', 'd)'), common.legend = TRUE, legend = 'bottom')

##########################################
# 6. Highlight likely candidate variants #
##########################################

# Get putative SNPs
#* Variant must be highly likely ot be pathogenic based on RegulomedB and have be
#* highly significant.
SNP_list = regAll %>% filter(rank <= 2 & -log(P) >= 5) %>%
SNP_list = SNP_list %>% arrange(ranking)
SNP_list = select(SNP_list, rsids, chrom, BP, ranking, DIR, P)
SNP_list



