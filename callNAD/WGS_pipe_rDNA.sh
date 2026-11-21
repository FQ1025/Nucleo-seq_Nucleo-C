#!/bin/bash
i=$1
## trim and fastqc
trim_galore --paired -j 20 --gzip --fastqc -q 20 -o cleandata rawdata/${i}_R1.fq.gz rawdata/${i}_R2.fq.gz

## bowtie2 alignment and mapq30 filteration
bowtie2 -p 20 -x /lustre1/lch3000_pkuhpc/fangq/ref/hg19_UCSC/Sequence/WholeGenomeFasta_rDNA/genome_rDNA \
                                -1 cleandata/${i}_R1_val_1.fq.gz -2 cleandata/${i}_R2_val_2.fq.gz | samtools view -@ 20 -h | samtools sort -@ 20 \
                                -T alignment/${i}_rDNA -O bam -o alignment/${i}_rDNA.sorted.bam

samtools view -@ 20 -h -S -q 30 alignment/${i}_rDNA.sorted.bam | grep -v chrM | samtools view -@ 20 -bS -o alignment/${i}_rDNA.rmM.MAPQ30.sorted.bam

cd alignment
samtools index -b ${i}_rDNA.rmM.MAPQ30.sorted.bam
cd ..
/lustre1/lch3000_pkuhpc/fangq/softwares/jdk-21.0.3/bin/java -jar /lustre1/lch3000_pkuhpc/fangq/softwares/picard.jar MarkDuplicates \
			I=alignment/${i}_rDNA.rmM.MAPQ30.sorted.bam \
			O=alignment/${i}_rDNA.rmdup.rmM.MAPQ30.sorted.bam \
			M=alignment/${i}_rDNA.dup.matrics.txt REMOVE_DUPLICATES=true

### quality input of bam file
##samtools flagstat
#samtools flagstat -@ 20 alignment/${i}_rDNA.rmdup.rmM.MAPQ30.sorted.bam > alignment/${i}_rDNA.rmdup.rmM.MAPQ30.sorted.stat

cd alignment
samtools index -b ${i}_rDNA.rmdup.rmM.MAPQ30.sorted.bam
cd ..

samtools idxstats -@ 10 alignment/${i}_rDNA.rmdup.rmM.MAPQ30.sorted.bam > ${i}_counts.stats
samtools idxstats -@ 10 alignment/${i}_rDNA.sorted.bam > ${i}_mtcounts.stats
cat ${i}_counts.stats | awk '{chrom[$1]=$1; reads[$1]=$3; total+=$3} END {for (chr in reads) print chr, (reads[chr] / total) * 100 "%"}' | sort -t 'r' -k2,2n -k1,1 > ${i}.cov 
##bam to bed
bamCoverage -p 20 --normalizeUsing RPKM --binSize 50 -b alignment/${i}_rDNA.rmdup.rmM.MAPQ30.sorted.bam -o bigwig/${i}_rDNA.bw -of bigwig
#bigWigToBedGraph bigwig/${i}.bw bedgraph/${i}.bdg
### get 1k coverage via samtools
#samtools bedcov -Q 30 /lustre1/lch3000_pkuhpc/fangq/projects/NAD/hg19_1k_rDNA_rmCent.bed alignment/${i}_rDNA.rmdup.rmM.MAPQ30.sorted.bam > bedgraph/${i}_rDNA_samcov_1k.bdg
#bedtools coverage -a hg19_1k.bed -b alignment/${i}.rmdup.rmM.MAPQ30.sorted.bam > bedgraph/${i}_bedcov_1k.bdg
