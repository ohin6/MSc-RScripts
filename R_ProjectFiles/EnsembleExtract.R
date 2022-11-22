###############
# Description #
###############

#* Script extracts selected Ensemble information from website, at specifed
#* genomic region.
#* 
#* This involves 
#* 
#* 1. Specify genome build
#* 2. Connecting to Ensemble database
#* 3. Identify attributes stored in Ensemble db
#* 4. Identify potential filters when extracting data
#* 5. Build query for data extraction
#* 6. Tidydata
#* 
#* This script extracts:
    #* Gene Id
    #* transcript id
    #* Gene symbol
    #* chrm number
    #* transcript start site
    #* transcript end
    #* transcript length
    #* 
    #* 
#* Author: Owen Williams

############
# Packages #
############
BiocManager::install("biomaRt")
require(biomaRt)
require(tidyverse)

###########################
# 1. Specify genome Build #
###########################

listEnsembl()
listEnsembl(GRCh=37)

###################################
# 2. Connect to Ensemble database #
###################################

# Connect to Ensemble dataset using GRch37 build
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)

# Specify host region
listMarts(host="uswest.ensembl.org")
ensembl_us_west = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://uswest.ensembl.org")

###################################
# 3. Identify Ensemble attributes #
###################################

# List attributes that can be extracted
listAttributes(ensembl)

################################
# 4. Identify Ensemble Filters #
################################

# Identify from list what to filter for
listFilters(ensembl)

##################
# 5. Build query #
##################
#* Extract Ensemble gene and Transcript names around APOE locus
#* filter based gene name


APOE_Locus = getBM(attributes = c('ensembl_gene_id','ensembl_transcript_id', 'ensembl_exon_id',
                                  'hgnc_symbol','chromosome_name','transcript_start',
                                  'transcript_end', 'transcript_length'),
                   filters = c('start','end', 'chromosome_name'),
                   values = list(start = 45350000, end = 45430000, chromosome_name = 19),
                   mart = ensembl)

###############
# 6. Tidydata #
###############

APOE_Locus = APOE_Locus %>%
  select(chromosome_name, transcript_start, transcript_end, ensembl_gene_id,
         ensembl_transcript_id)
