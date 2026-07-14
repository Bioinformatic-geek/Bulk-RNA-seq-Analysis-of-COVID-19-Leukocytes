#!/bin/bash

# ============================================================
# RNA-Seq Pipeline: FASTQ to Count Matrix
# Dataset: GSE157103 (COVID-19 vs Healthy, leukocytes)
# Tools: FastQC, Trimmomatic, HISAT2, SAMtools, featureCounts
# ============================================================

# --- Paths ---
PROJECT=~/rnaseq_project
RAW=$PROJECT/raw_data
TRIMMED=$PROJECT/trimmed
ALIGNED=$PROJECT/aligned
COUNTS=$PROJECT/counts
REF=$PROJECT/reference
QC=$PROJECT/fastqc_reports
THREADS=4   # change this to match your nproc

# --- Sample list ---
SAMPLES=(SRR12544419 SRR12544420 SRR12544421
         SRR12544527 SRR12544528 SRR12544529)
# ============================================================
# STEP 1: Quality Control on raw reads
# ============================================================
echo "=== STEP 1: FastQC on raw reads ==="
fastqc $RAW/*.fastq.gz -o $QC/raw/ -t $THREADS
echo "FastQC complete. Check $QC/raw/ for reports."

# ============================================================
# STEP 2: Trimming with Trimmomatic
# ============================================================
echo "=== STEP 2: Trimming ==="

for sample in "${SAMPLES[@]}"; do
    echo "Trimming $sample..."
    trimmomatic PE \
      $RAW/${sample}_1.fastq.gz $RAW/${sample}_2.fastq.gz \
      $TRIMMED/${sample}_1_paired.fastq.gz $TRIMMED/${sample}_1_unpaired.fastq.gz \
      $TRIMMED/${sample}_2_paired.fastq.gz $TRIMMED/${sample}_2_unpaired.fastq.gz \
      ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
      LEADING:3 TRAILING:3 \
      SLIDINGWINDOW:4:15 \
      MINLEN:36 \
      -threads $THREADS
    echo "Done: $sample"
done

# FastQC on trimmed reads
fastqc $TRIMMED/*_paired.fastq.gz -o $QC/trimmed/ -t $THREADS

# ============================================================
# STEP 3: Alignment with HISAT2
# ============================================================
echo "=== STEP 3: Alignment ==="

for sample in "${SAMPLES[@]}"; do
    echo "Aligning $sample..."
    hisat2 -x $REF/hg38_index \
      -1 $TRIMMED/${sample}_1_paired.fastq.gz \
      -2 $TRIMMED/${sample}_2_paired.fastq.gz \
      -S $ALIGNED/${sample}.sam \
      --dta \
      -p $THREADS \
      2> $ALIGNED/${sample}_alignment_stats.txt  # saves alignment rate to file

    # Convert to BAM, sort, index
    samtools sort $ALIGNED/${sample}.sam -o $ALIGNED/${sample}.bam -@ $THREADS
    samtools index $ALIGNED/${sample}.bam

    # Delete SAM immediately to save space
    rm $ALIGNED/${sample}.sam

    echo "Done: $sample"
done

# ============================================================
# STEP 4: Quantification with featureCounts
# ============================================================
echo "=== STEP 4: featureCounts ==="

featureCounts \
  -T $THREADS \
  -p \
  -a $REF/Homo_sapiens.GRCh38.109.gtf \
  -o $COUNTS/raw_counts.txt \
  $ALIGNED/*.bam

echo "=== Pipeline complete! ==="
