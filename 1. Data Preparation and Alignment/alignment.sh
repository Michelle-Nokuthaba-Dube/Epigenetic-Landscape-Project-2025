# 1. Trim reads
fastp -i SRR11392452_1.fastq -o SRR11392452_trimmed.fastq
fastp -i SRR11392444_1.fastq -o SRR11392444_trimmed.fastq

# 2. Align to mm10
bowtie2 -x ../genome/mm10 -U SRR11392452_trimmed.fastq \
  -S ../alignments/SRR11392452.sam 2> ../logs/SRR11392452.align.log

bowtie2 -x ../genome/mm10 -U SRR11392444_trimmed.fastq \
  -S ../alignments/SRR11392444.sam 2> ../logs/SRR11392444.align.log

# 3. Convert SAM → BAM
samtools view -bS ../alignments/SRR11392452.sam > ../alignments/SRR11392452.bam
samtools view -bS ../alignments/SRR11392444.sam > ../alignments/SRR11392444.bam

# 4. Sort BAM
samtools sort ../alignments/SRR11392452.bam -o ../alignments/SRR11392452.sorted.bam
samtools sort ../alignments/SRR11392444.bam -o ../alignments/SRR11392444.sorted.bam

# 5. Index BAM
samtools index ../alignments/SRR11392452.sorted.bam
samtools index ../alignments/SRR11392444.sorted.bam
