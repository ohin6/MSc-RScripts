############
# Packages #
############

require(Gviz)
require(GenomicRanges)
require(tidyverse)
require(biomaRt)
require(mapsnp)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)

###################################
# 1. Connect to Ensemble database #
###################################

ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)

####################
# 2. Create tracks #
####################

# Gene Track
#* This extracts data from online ensembl db (need internet connection and set 
#* up link - see above)
biomTrack = BiomartGeneRegionTrack(genome = "hg19", chromosome = 'chr19', 
                                    start = 45400000, end = 45420000,
                                    name = "ENSEMBL", biomart = ensembl,
                                    transcriptAnnotation = "symbol")

# Genome axis Track
axisTrack = GenomeAxisTrack(showId = TRUE)

# Genome axis Track 2
ideoTrack = IdeogramTrack(genome = "hg19", chromosome = "chr19", showBandId = TRUE)

# Annotation track
#* plotting variants (rs variant ids) linked to AD on chr19
RsTrack = AnnotationTrack(start=c(45351516, 45387459, 45396144, 45421204, 45421254,
                            45422160,45424514),
                    end=c(45351516, 45387459, 45396144, 45421204, 45421254,
                          45422160,45424514),
                    chromosome = 'chr19',
                    id=c('rs41289512', 'rs12972156', 'rs11556505', 'rs150966173',
                         'rs12721046', 'rs12721051', 'rs157592'),
                    genome = 'hg19', stacking="squish", name="Variant")

##################
# 3. Plot tracks #
##################

plotTracks(list(ideoTrack, axisTrack, biomTrack, x),
           background.panel = "#FFFEDB", background.title = "darkblue")

SNP = dplyr::tibble(Chr = '19',
             id=c('rs41289512', 'rs12972156', 'rs11556505', 'rs150966173',
                       'rs12721046', 'rs12721051', 'rs157592'),
             bp=c(45351516, 45387459, 45396144, 45421204, 45421254,
                     45422160,45424514))

x = msb(M = SNP, start = 45340000, end = 45435000, showLab.gene = FALSE,
    showBandId = TRUE, rotation.item = 90, snpWd = 0.1, IDWd = 0.05,
    fontsize.chr = 10, ) 

displayPars(msb(M = SNP, start = 45340000, end = 45435000, showLab.gene = FALSE,
                showBandId = TRUE, rotation.item = 90, snpWd = 0.1, IDWd = 0.05,
                fontsize.chr = 10, cex = 0.5, stackHt.snp = 0.3))



