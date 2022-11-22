###############
# Description #
###############
#* The purpose of this script is identify loci strongly associated with AD from
#* GWAS data, as well as performing quality control checks on data.
#* 
#* GWAS data is from Marioni et al., 2018
#* 
#* This involves
#* 
#* 1. Import GWAS dataset
#* 2. Data wrangling
#* 3. Create Manhatten plots
#* 4. Create QQplots for each chromosome
#* 4. Integrate GWAS with RegulomedB variant information
#* 5. Create Plots stratified by SNP pathogenic rank
#* 6. Highlight likely candidate variants

####################
# Install packages #
####################

#Install and open manhatten plot
install.packages("qqman")
library(qqman)
library(tidyverse)

############################
# 1. Download GWAS dataset #
############################

# GWAS data from Marioni et al., 2018
df = read.table("../Raw_data/3_UKB_AD_parental_meta_summary_output_June2019.txt", header = T)

#####################
# 2. data Wrangling #
#####################

# Modify P-values that are 0 to prevent -log10(P) = inf
df$P[abs(df$P)<1e-999]=1e-999

# check not getting inf values
df$nan = -log10(df$P)

#############################
# 3. Create Manhattan Plots #
#############################
# Warning this can take a while due to file size

# Make the Manhattan plot on the gwas Results dataset
manhattan(df, main = "Manhattan Plot", ylim = c(0, 20), cex = 0.6, cex.axis = 0.5, 
          col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F)


######################
# 4. Create QQ Plots #
######################

chrm1 = subset(df, CHR==1)
png("chrm1.png")
qq(subset(df, CHR == 1)$P)
png("chrm2.png")
qq(subset(df, CHR == 2)$P)
png("chrm3.png")
qq(subset(df, CHR == 3)$P)
png("chrm4.png")
qq(subset(df, CHR == 4)$P)
png("chrm5.png")
qq(subset(df, CHR == 5)$P)
png("chrm6.png")
qq(subset(df, CHR == 6)$P)
png("chrm7.png")
qq(subset(df, CHR == 7)$P)
png("chrm8.png")
qq(subset(df, CHR == 8)$P)
png("chrm9.png")
qq(subset(df, CHR == 9)$P)
png("chrm10.png")
qq(subset(df, CHR == 10)$P)
png("chrm11.png")
qq(subset(df, CHR == 11)$P)
png("chrm12.png")
qq(subset(df, CHR == 12)$P)
png("chrm13.png")
qq(subset(df, CHR == 13)$P)
png("chrm14.png")
qq(subset(df, CHR == 14)$P)
png("chrm15.png")
qq(subset(df, CHR == 15)$P)
png("chrm16.png")
qq(subset(df, CHR == 16)$P)
png("chrm17.png")
qq(subset(df, CHR == 17)$P)
png("chrm18.png")
qq(subset(df, CHR == 18)$P)
png("chrm19.png")
qq(subset(df, CHR == 19)$P)
png("chrm20.png")
qq(subset(df, CHR == 20)$P)
png("chrm21.png")
qq(subset(df, CHR == 21)$P)
png("chrm22.png")
qq(subset(df, CHR == 22)$P)

