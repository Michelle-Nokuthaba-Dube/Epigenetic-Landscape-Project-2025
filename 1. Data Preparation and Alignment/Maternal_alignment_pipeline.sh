Raw data retrieval

prefetch SRR11392448
fasterq-dump SRR11392448 --split-files -O fastq/

Quality trimming
fastp \
  -i fastq/SRR11392448.fastq \
  -o fastq/SRR11392448_trimmed.fastq \
  --html SRR11392448_fastp.html \
  --json SRR11392448_fastp.json \
  --thread 4

  Alignment to mm10 (using Bowtie2)
  bowtie2 \
  -x genome/mm10 \
  -U fastq/SRR11392448_trimmed.fastq \
  -S alignments/SRR11392448.sam \
  --threads 4

  Bam conversion, sorting, indexing (SAMtools)
  samtools view -bS SRR11392448.sam > SRR11392448.bam

  samtools sort SRR11392448.bam -o SRR11392448_sorted.bam

  samtools index SRR11392448_sorted.bam

  samtools flagstat SRR11392448_sorted.bam > SRR11392448_flagstat.txt

