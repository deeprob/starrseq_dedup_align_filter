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
while getopts i:r:o:d:f: flag
do
    case "${flag}" in
        i) INPUT_PREFIX=${OPTARG};;
        r) BIOL_REPS=(${OPTARG});;
        o) OUTPUT_PREFIX=${OPTARG};;
        f) ROI_FILE=${OPTARG};;
        d) FLAG_DUP=${OPTARG};;
    esac
done

for rep in "${BIOL_REPS[@]}"
do
    # mark duplicates using picard if this option is true: step 1
    if [ $FLAG_DUP == true ]
    then
        picard MarkDuplicates -I ${INPUT_PREFIX}_${rep}.bam -O ${OUTPUT_PREFIX}_${rep}_dup.bam -M ${OUTPUT_PREFIX}_${rep}_dup_metrics.txt --VALIDATION_STRINGENCY LENIENT --REMOVE_DUPLICATES true --VERBOSITY WARNING
        # use samtools to filter, get only the regions of interests, filter options taken directly from starrpeaker github page: step 2
        samtools view -F 3852 -f 2 -q 30 -L ${ROI_FILE} -o ${OUTPUT_PREFIX}_${rep}.bam ${OUTPUT_PREFIX}_${rep}_dup.bam 
    else
        echo "Duplicates removed by the user"
        samtools view -F 2828 -f 2 -q 30 -L ${ROI_FILE} -o ${OUTPUT_PREFIX}_${rep}.bam ${INPUT_PREFIX}_${rep}.bam 
    fi 

done

# merge all replicates into a file and transfer it to the output library: step 3
TMP_PRE=(${BIOL_REPS[@]/#/${OUTPUT_PREFIX}_})
TMP=${TMP_PRE[@]/%/.bam}
samtools merge -f ${OUTPUT_PREFIX}.bam $TMP
