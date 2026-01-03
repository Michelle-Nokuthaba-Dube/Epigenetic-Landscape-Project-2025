--------------------------------------------------------------------------------
  ###EPIGENETIC LANDSCAPE OF STRESS-RELATED GENES IN MICE BRAIN TISSUE###
--------------------------------------------------------------------------------
  
# STEP 3: Read Quantification in Peak regions
# Completed by Michelle Dube
  
--------------------------------------------------------------------------------
# Checking working directory
getwd()
[1] "C:/Users/User/Desktop/ICBB Course/Group Project"

# Set working directory
setwd("C:/Users/User/Desktop/ICBB Course/Group Project/3. Read Quantification in Peak regions")

# Organising Project
File menu > New Project > New Directory > New Project > Epigenetic_Landscape_Project

# Install required packages
install.packages("BiocManager")
BiocManager::install("ChIPseeker")
BiocManager::install(c("GenomicFeatures", "rtracklayer"))

# Install annotation package
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("org.Mm.eg.db")

# Load packages
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

# Set file paths
broad_file <- "C:/Users/User/Desktop/ICBB Course/Group Project/2. Peak calling and Annotation/macs2_output/MACS2_output_broad_peaks.bed"
filtered_file <- "C:/Users/User/Desktop/ICBB Course/Group Project/2. Peak calling and Annotation/macs2_output/sample_peaks_filtered.bed"

file.exists(broad_file)
file.exists(filtered_file)

# Read the BED files
broad_peaks <- readPeakFile(broad_file)
filtered_peaks <- readPeakFile(filtered_file)

# Check the loaded data
broad_peaks
filtered_peaks

# Load TxDb and annotation database
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Annotate broad peaks
anno_broad <- annotatePeak(
  broad_peaks,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),  # ±3 kb around TSS
  annoDb = "org.Mm.eg.db"
)

# Annotate filtered peaks
anno_filtered <- annotatePeak(
  filtered_peaks,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Mm.eg.db"
)

# View annotation summary
# Pie chart of peak distribution
plotAnnoPie(anno_broad)
plotAnnoPie(anno_filtered)

# Barplot comparison
plotAnnoBar(list(Broad = anno_broad, Filtered = anno_filtered))

# Convert to data frame & export
df_broad <- as.data.frame(anno_broad)
df_filtered <- as.data.frame(anno_filtered)

write.csv(df_broad, "C:/Users/User/Desktop/ICBB Course/Group Project/3. Read Quantification in Peak regions/macs2_output/annotated_broad_peaks.csv", row.names = FALSE)
write.csv(df_filtered, "C:/Users/User/Desktop/ICBB Course/Group Project/3. Read Quantification in Peak regions/macs2_output/annotated_filtered_peaks.csv", row.names = FALSE)