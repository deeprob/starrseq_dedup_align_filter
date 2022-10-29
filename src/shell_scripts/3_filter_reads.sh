#!/bin/bash
set -ue
# filter and merge the aligned but unfiltered reads from replicates obtained after aligning fastq to reference
# filtration steps
# 1) Mark duplicates using picard under lenient setting
# 2) Get rid of unmapped, secondary alignments, mapping quality score < 30, reads using samtools
# 3) Merge aligned and filtered biological replicates

# example command
# $ bash filter_reads.sh -i "/data5/deepro/starrseq/miseq_test/aligned_reads/input" -r "S1" -o "/data5/deepro/starrseq/miseq_test/filtered_libraries/input" -d false

# get all arguments
while getopts i:r:o:f: flag
do
    case "${flag}" in
        i) INPUT_PREFIX=${OPTARG};;
        r) BIOL_REPS=(${OPTARG});;
        o) OUTPUT_PREFIX=${OPTARG};;
        f) ROI_FILE=${OPTARG};;
    esac
done

for rep in "${BIOL_REPS[@]}"
do

    if [ -z $ROI_FILE]
    then
        samtools index ${INPUT_PREFIX}_${rep}.bam
        samtools view -F 2828 -f 2 -q 30 -@ 64 -o ${OUTPUT_PREFIX}_${rep}.bam ${INPUT_PREFIX}_${rep}.bam chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr2 chr20 chr21 chr22 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrX chrY
        samtools sort -@ 64 -o ${OUTPUT_PREFIX}_${rep}_tmp.bam ${OUTPUT_PREFIX}_${rep}.bam
        mv ${OUTPUT_PREFIX}_${rep}_tmp.bam ${OUTPUT_PREFIX}_${rep}.bam
    else
        samtools index ${INPUT_PREFIX}_${rep}.bam
        samtools view -F 2828 -f 2 -q 30 -@ 64 -L ${ROI_FILE} -o ${OUTPUT_PREFIX}_${rep}.bam ${INPUT_PREFIX}_${rep}.bam chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr2 chr20 chr21 chr22 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chrX chrY
        samtools sort -@ 64 -o ${OUTPUT_PREFIX}_${rep}_tmp.bam ${OUTPUT_PREFIX}_${rep}.bam
        mv ${OUTPUT_PREFIX}_${rep}_tmp.bam ${OUTPUT_PREFIX}_${rep}.bam
    fi

done

# merge all replicates into a file and transfer it to the output library: step 3
TMP_PRE=(${BIOL_REPS[@]/#/${OUTPUT_PREFIX}_})
TMP=${TMP_PRE[@]/%/.bam}
samtools merge -f ${OUTPUT_PREFIX}.bam $TMP
