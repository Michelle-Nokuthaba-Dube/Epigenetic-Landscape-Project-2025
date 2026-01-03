--------------------------------------------------------------------------------
###EPIGENETIC LANDSCAPE OF STRESS-RELATED GENES IN MICE BRAIN TISSUE###
--------------------------------------------------------------------------------
  
# STEP 4: Functional Enrichment Analysis
# Completed by Michelle Dube
  
--------------------------------------------------------------------------------
# Checking working directory
getwd()
[1] "C:/Users/User/Desktop/ICBB Course/Group Project"

# Gene Ontology (GO) Biological Processes for Histone modification

# Load clusterProfiler
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
library(clusterProfiler)

# Load packages
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

# Read MACS2 peak file
broad_file <- "C:/Users/User/Desktop/ICBB Course/Group Project/macs2_output/2. Peak calling and Annotation/MACS2_output_broad_peaks.bed"

# Broad peaks
broad_file <- "C:/Users/User/Desktop/ICBB Course/Group Project/2. Peak calling and Annotation/macs2_output/MACS2_output_broad_peaks.bed"

# Filtered peaks
filtered_file <- "C:/Users/User/Desktop/ICBB Course/Group Project/2. Peak calling and Annotation/macs2_output/sample_peaks_filtered.bed"


broad_peaks <- readPeakFile(broad_file)

# Check if loaded
broad_peaks

# Annotate peaks
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

anno_broad <- annotatePeak(
  broad_peaks,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Mm.eg.db"
)

# Check the annotation
anno_broad

# Convert annotated peaks to data frame
df_broad <- as.data.frame(anno_broad)
head(df_broad)

# Extract gene IDs
genes_broad <- df_broad$geneId

# Extract gene IDs from annotated peaks

# For broad peaks
genes_broad <- as.data.frame(anno_broad)$geneId

# Remove duplicates
genes_broad <- unique(genes_broad)

# Check the result
head(genes_broad)
length(genes_broad)

# Load org.Mm.eg.db

# Check
ls("package:org.Mm.eg.db")  # should list functions/data

# Run GO enrichment (Biological Process)
library(clusterProfiler)

ego_broad <- enrichGO(
  gene          = genes_broad,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",          # Biological Process
  pAdjustMethod = "BH",          # Benjamini-Hochberg correction
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE           # Converts Entrez IDs to gene symbols
)

# Identify top stress-related pathways
library(ggplot2)
library(dplyr)

ls()

# Check exact file names
list.files("C:/Users/User/Desktop/ICBB Course/Group Project/2. Peak calling and Annotation/macs2_output")

# Load the filtered peaks
filtered_file <- "C:/Users/User/Desktop/ICBB Course/Group Project/2. Peak calling and Annotation/macs2_output/sample_peaks_filtered.bed"
filtered_peaks <- readPeakFile(filtered_file)
filtered_peaks

# Annotate the filtered peaks
anno_filtered <- annotatePeak(filtered_peaks, 
                              TxDb = txdb, 
                              tssRegion = c(-1000, 1000),
                              annoDb = "org.Mm.eg.db")

# Convert to data frame
df_filtered <- as.data.frame(anno_filtered)

# Check
head(df_filtered)

# Extract genes IDs
genes_filtered <- unique(df_filtered$geneId)
length(genes_filtered)

# Run GO enrichment for filtered peaks
library(clusterProfiler)
library(org.Mm.eg.db)

ego_filtered <- enrichGO(
  gene          = genes_filtered,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)

# Load filtered peak
library(ChIPseeker)

filtered_file <- "C:/Users/User/Desktop/ICBB Course/Group Project/2. Peak calling and Annotation/macs2_output/sample_peaks_filtered.bed"

filtered_peaks <- readPeakFile(filtered_file)

# Annotate filtered peaks
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

anno_filtered <- annotatePeak(
  filtered_peaks,
  tssRegion = c(-1000, 1000),
  TxDb = txdb
)

# Extract genes from filtered peaks
df_filtered <- as.data.frame(anno_filtered)

genes_filtered <- unique(df_filtered$geneId)

length(genes_filtered)   # just to check

# Run GO for filtered peaks
library(clusterProfiler)
library(org.Mm.eg.db)

ego_filtered <- enrichGO(
  gene          = genes_filtered,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# Convert to data frame
ego_filtered_df <- as.data.frame(ego_filtered)
head(ego_filtered_df)


# Extract stress-related pathways
library(dplyr)

stress_filtered <- ego_filtered_df %>%
  filter(grepl("stress", Description, ignore.case = TRUE))

stress_filtered

head(df_filtered)
length(genes_filtered)
head(ego_filtered_df)

# Inspect top stress-related pathways
head(stress_filtered, 10)  # Top 10 stress-related GO terms

# Visualise stress-related GO terms
# Barplot of top stress GO terms
library(ggplot2)

# Pick top 10 stress GO terms by FoldEnrichment
top_stress <- stress_filtered %>%
  arrange(desc(FoldEnrichment)) %>%
  head(10)

ggplot(top_stress, aes(x = reorder(Description, FoldEnrichment), y = FoldEnrichment)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top Stress-Related GO Terms (Filtered Peaks)",
       x = "GO Term",
       y = "Fold Enrichment")

# Count number of genes per stress GO term
top_stress %>%
  select(Description, Count)

# Maternal stress pathways from both broad and filtered peaks
# For broad peaks
ego_broad_df <- as.data.frame(ego_broad)

# For filtered peaks
ego_filtered_df <- as.data.frame(ego_filtered)

stress_broad <- ego_broad_df %>%
  filter(grepl("stress", Description, ignore.case = TRUE))

ls()

# Filter stress-related GO terms
library(dplyr)

stress_broad <- ego_broad_df %>%
  filter(grepl("stress", Description, ignore.case = TRUE))

stress_filtered <- ego_filtered_df %>%
  filter(grepl("stress", Description, ignore.case = TRUE))

# Visualise results
head(stress_broad)
head(stress_filtered)

# Load required libraries
library(dplyr)
library(ggplot2)

# Combine stress-related GO pathways from broad and filtered peaks
# Step 1: Convert GO enrichment results to data frames
ego_broad_df <- as.data.frame(ego_broad)
ego_filtered_df <- as.data.frame(ego_filtered)

# Step 2: Filter for stress-related pathways
stress_broad <- ego_broad_df %>%
  filter(grepl("stress", Description, ignore.case = TRUE))

stress_filtered <- ego_filtered_df %>%
  filter(grepl("stress", Description, ignore.case = TRUE))

# Step 3: Add a column indicating the peak type
stress_broad$PeakType <- "Broad"
stress_filtered$PeakType <- "Filtered"

# Step 4: Combine both data frames
stress_combined <- rbind(stress_broad, stress_filtered)

# Step 5: Optional – keep top 10 pathways per peak type (by Count or FoldEnrichment)
stress_combined <- stress_combined %>%
  group_by(PeakType) %>%
  slice_max(order_by = Count, n = 10) %>%
  ungroup()

# Step 6: Visualize using ggplot2
ggplot(stress_combined, aes(x = reorder(Description, Count), y = Count, fill = PeakType)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip()

library(dplyr)
library(ggplot2)

str(stress_combined)

stress_combined <- stress_combined %>%
  mutate(
    Description = sapply(Description, `[`, 1),  # take first element if list
    Count = as.numeric(Count)                  # ensure Count is numeric
  )

library(dplyr)

# Convert columns to proper types
stress_combined <- stress_combined %>%
  mutate(
    Description = as.character(Description),
    Count = as.numeric(unlist(Count))  # unlist if it's a list
  )

# Check structure
str(stress_combined)

# Now compute top 5 pathways
top5_pathways <- stress_combined %>%
  group_by(Description) %>%
  summarise(TotalCount = sum(Count, na.rm = TRUE)) %>%
  arrange(desc(TotalCount)) %>%
  slice_head(n = 5)

top5_pathways


library(dplyr)

# Convert the relevant columns to proper types
stress_combined <- stress_combined %>%
  mutate(
    Description = as.character(Description),
    Count = as.numeric(unlist(Count))  # Convert Count from list/Rle to numeric
  )

# Check structure to confirm
str(stress_combined)

# Compute top 5 pathways
top5_pathways <- stress_combined %>%
  group_by(Description) %>%
  summarise(TotalCount = sum(Count, na.rm = TRUE)) %>%
  arrange(desc(TotalCount)) %>%
  slice_head(n = 5)  # safer than slice(1:5)

# View results
top5_pathways

# Create stress_top5
stress_top5 <- stress_combined %>%
  inner_join(top5_pathways, by = "Description")

# Check
head(stress_top5)
str(stress_top5)

# Plot
ggplot(stress_top5, aes(x = reorder(Description, TotalCount), y = Count, fill = PeakType)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Top 5 Maternal Stress-related Pathways",
       x = "GO Term",
       y = "Gene Count") +
  theme_minimal()


ggplot(stress_top5, aes(x = reorder(Description, TotalCount), y = Count, fill = PeakType)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.9), hjust = -0.1) +
  coord_flip() +
  labs(title = "Top 5 Maternal Stress-related Pathways",
       x = "GO Term",
       y = "Gene Count") +
  theme_minimal()

library(ggplot2)

# Combine stress-related GO terms from broad and filtered peaks
library(dplyr)

# Check ego_broad_df and ego_filtered_df exist
stress_broad <- ego_broad_df %>%
  filter(grepl("stress", Description, ignore.case = TRUE)) %>%
  mutate(PeakType = "Broad")

stress_filtered <- ego_filtered_df %>%
  filter(grepl("stress", Description, ignore.case = TRUE)) %>%
  mutate(PeakType = "Filtered")

# Combine both datasets
stress_combined <- bind_rows(stress_broad, stress_filtered)

# Add a column for total gene count per GO term across peak types
stress_combined <- stress_combined %>%
  group_by(Description) %>%
  mutate(TotalCount = sum(Count)) %>%
  ungroup()

# Select top 5 maternal stress-related pathways
stress_top5 <- stress_combined %>%
  arrange(desc(TotalCount)) %>%
  distinct(Description, .keep_all = TRUE) %>%
  slice_head(n = 5)

# Visualise with a barplot
library(ggplot2)

ggplot(stress_top5, aes(x = reorder(Description, TotalCount), y = Count, fill = PeakType)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "Top 5 Maternal Stress-related Pathways",
       x = "GO Term",
       y = "Gene Count") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

# Save the results
write.csv(stress_top5, "stress_top5_pathways.csv", row.names = FALSE)
