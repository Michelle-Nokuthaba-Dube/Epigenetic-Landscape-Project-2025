# Download mm10 reference genome
curl -L -o mm10.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz

# Unzip
gunzip mm10.fa.gz

# Build Bowtie2 index
bowtie2-build mm10.fa mm10
