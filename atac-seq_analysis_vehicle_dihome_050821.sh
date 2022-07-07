# Performing the ATAC-seq analysis for the 12-13 DiHOME study interactively

module load java-jdk/1.10.0_1
module load fastqc/0.11.7
module load gcc/6.2.0
module load multiqc/1.8.0
module load java-jdk/1.8.0_92
module load trimmomatic/0.36
module load intel/2017
module load bowtie2/2.3.4.3

# sam to bam and remove mtdna reads

module load gcc/6.2.0
module load samtools/1.10

# MarkDuplicates

module load java-jdk/1.8.0_92
module load picard/2.18.29

# call the peaks

module load gcc/6.2.0
module load macs2/2.1.0

# remove blacklisted regions with bedtools

module load gcc/6.2.0
module load bedtools/2.29.0

# read counts for each annotation using featurecounts

module load gcc/6.2.0
module load subread/1.5.3

# annotate the peaks using homer

module load gcc/6.2.0
module load homer/4.10.0

# 1) For each set of paired fastq files, ran fastqc, concatenated the fastq files, and reran fastqc, then performed multiqc to summarize all .html files.

fastqc [file]
multiqc [directory]

# 2) Use cutaadpt to trim the Nextera adapters that are present on a significant proportion of reads (sample-specific pbs scripts)

java -jar /apps/software/java-jdk-1.8.0_92/trimmomatic/0.36/trimmomatic-0.36.jar PE -threads 6 -trimlog SampleTA1.txt -phred33 Sample-1_S1_L001_R1_001.fastq.gz Sample-1_S1_L001_R2_001.fastq.gz SampleTA1_R1_paired.fastq.gz SampleTA1_R1_unpaired.fastq.gz SampleTA1_R2_paired.fastq.gz SampleTA1_R2_unpaired.fastq.gz ILLUMINACLIP:/apps/software/java-jdk-1.8.0_92/trimmomatic/0.36/adapters/NexteraPE-PE.fa:2:30:10:2:keepBothReads

## Checked a couple files with fastqc and indeed the adapters are trimmed

# 3) Align the reads using Bowtie2 (sample-specific pbs scripts)

bowtie2 --threads 8 -x /scratch/kmagnaye/12-13/hg19/hg19 --very-sensitive -X 2000 -1 /scratch/kmagnaye/12-13/Raw_Data/SampleTC4/SampleTA1_R1_paired.fastq.gz -2 /scratch/kmagnaye/12-13/Raw_Data/SampleTC4/SampleTA1_R2_paired.fastq.gz -S /scratch/kmagnaye/12-13/Mapped_Paired_Data/SampleTA1.sam

# 4) Convert sam to bam and then sort

samtools view --threads 7 -bS SampleTA1.sam | samtools sort --threads 7 - > SampleTA1.sorted.bam

# 5) Index the BAM file (a few minutes)

samtools index SampleTA1.sorted.bam

# 6) Remove reads aligning to the mitochondria

samtools idxstats SampleTA1.sorted.bam | cut -f1 | grep -v chrM | xargs samtools view --threads 7 -b SampleTA1.sorted.bam > SampleTA1.sorted.nomt.bam

# 7) Remove PCR duplicates using the MarkDuplicates function in Picard

## mark the pcr duplicates

java -jar /apps/software/java-jdk-1.8.0_92/picard/2.18.29/picard.jar MarkDuplicates \
  QUIET=true INPUT=SampleTA1.sorted.nomt.bam OUTPUT=SampleTA1.sorted.nomt.marked.bam METRICS_FILE=SampleTA1.sorted.metrics \
  REMOVE_DUPLICATES=false CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=/tmp

## remove them 

samtools view --threads 7 -h -b -F 1024 SampleTA1.sorted.nomt.marked.bam > SampleTA1.sorted.nomt.marked.nodup.bam

## reindex

samtools index SampleTA1.sorted.nomt.marked.nodup.bam

# 8) Remove non-unique alignments (aka multi-mapped reads with MAPQ < 10)

samtools view --threads 7 -h -q 10 SampleTA1.sorted.nomt.marked.nodup.bam > SampleTA1.sorted.nomt.marked.nodup.rmmulti10.bam

sleep 1

samtools index SampleTA1.sorted.nomt.marked.nodup.rmmulti10.bam

# 9) Extract fragment insert sizes

samtools view SampleTA1.sorted.nomt.marked.nodup.rmmulti10.bam | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > SampleTA1_fragment_length_count.txt

# 10) Call peaks for all samples using MACS2

# setting the PYTHON_EGG_CACHE environment
export PYTHON_EGG_CACHE=/scratch/kmagnaye/Python-Eggs

# calling macs2 (narrow)
macs2 callpeak -t SampleTA1.sorted.nomt.marked.nodup.rmmulti10.bam --nomodel --shift -100 --extsize 200 -q 0.05 -f BAMPE -g hs --keep-dup all -n SampleTA1 -B --outdir .
# calling macs2 (broad)
macs2 callpeak -t SampleTA1.sorted.nomt.marked.nodup.rmmulti10.bam --nomodel --shift -100 --extsize 200 -q 0.05 --broad -f BAMPE -g hs --keep-dup all -n SampleTA1 -B --outdir .

# 11) Use bedtools intersect function to merge the bed files of the vehicle and dihome treated samples and merge peaks within 10 bp (keeps peaks if found in at least one individual)

cat SampleTA1_narrow_mapq10_peaks.narrowPeak.filt SampleTA2_narrow_mapq10_peaks.narrowPeak.filt SampleTA4_narrow_mapq10_peaks.narrowPeak.filt SampleTB1_narrow_mapq10_peaks.narrowPeak.filt SampleTB3_narrow_mapq10_peaks.narrowPeak.filt SampleTB4_narrow_mapq10_peaks.narrowPeak.filt SampleTC2_narrow_mapq10_peaks.narrowPeak.filt SampleTC3_narrow_mapq10_peaks.narrowPeak.filt > SampleT_vehicle_dihome_narrow_mapq10_filt
sort -k1,1 -k2,2n -k3,3n SampleT_vehicle_dihome_narrow_mapq10_filt > SampleT_vehicle_dihome_narrow_mapq10_filt_sort
cat SampleT_vehicle_dihome_narrow_mapq10_filt_sort | awk '{print $1"\t"$2"\t"$3"\t"$7}' > SampleT_vehicle_dihome_narrow_mapq10_filt_sort.bed
bedtools merge -d 10 -c 4 -o mean -i SampleT_vehicle_dihome_narrow_mapq10_filt_sort.bed > SampleT_vehicle_dihome_narrow_mapq10_filt_mergedpeaks.bed

sort -k1,1 -k2,2n -k3,3n SampleT_vehicle_dihome_narrow_mapq10_filt_mergedpeaks.bed > SampleT_vehicle_dihome_narrow_mapq10_filt_mergedpeaks_sorted.bed
cat SampleT_vehicle_dihome_narrow_mapq10_filt_mergedpeaks_sorted.bed | awk '{print $1"\t"$2"\t"$3}' > SampleT_vehicle_dihome_narrow_mapq10_filt_mergedpeaks_sorted_final.bed
awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' SampleT_vehicle_dihome_narrow_mapq10_filt_mergedpeaks_sorted_final.bed > SampleT_vehicle_dihome_narrow_mapq10_filt_mergedpeaks_sorted_final.saf

# 12) For each sample, count the number of reads for each annotation using featureCounts (subread package)

featureCounts --read2pos 5 -a SampleT_vehicle_dihome_narrow_mapq10_filt_mergedpeaks_sorted_final.saf -F SAF -T 7 -p -o All_peaks_countMatrix_vehicle_dihome.txt *rmmulti10.bam

# 13) Use HOMER to perform pathway enrichments and motif analyis of each of the six sets of regions (example of one)

/apps/software/gcc-6.2.0/homer/4.10.0/bin/annotatePeaks.pl gained_sites_strict_6_051021.bed hg19 > gained_sites_strict_6_051021.annotated

/apps/software/gcc-6.2.0/homer/4.10.0/bin/annotatePeaks.pl gained_sites_strict_6_051021.bed -go /scratch/kmagnaye/12-13/Mapped_Paired_Data/gained_sites_strict

/apps/software/gcc-6.2.0/homer/4.10.0/bin/findMotifsGenome.pl gained_sites_strict_6_051021.bed hg19 /scratch/kmagnaye/12-13/Mapped_Paired_Data/gained_sites_strict -bg /scratch/kmagnaye/12-13/Mapped_Paired_Data/SampleT_vehicle_dihome_narrow_mapq10_filt_mergedpeaks_sorted_final.bed





