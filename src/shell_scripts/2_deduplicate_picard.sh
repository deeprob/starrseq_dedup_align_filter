#!/bin/bash
set -ue

# 1) Mark duplicates using picard under lenient setting

# get all arguments
while getopts i:r:o: flag
do
    case "${flag}" in
        i) INPUTFILE_PREFIX=${OPTARG};;
        r) BIOL_REPS=(${OPTARG});;
        o) OUTPUTFILE_PREFIX=${OPTARG};;
    esac
done

picard_tmp=${OUTPUTFILE_PREFIX}_tmp

mkdir -p $picard_tmp 

for rep in "${BIOL_REPS[@]}"
do
    # mark duplicates using picard if this option is true: step 1
    picard MarkDuplicates -I ${INPUTFILE_PREFIX}_${rep}.bam -O ${OUTPUTFILE_PREFIX}_${rep}.bam -M ${OUTPUTFILE_PREFIX}_${rep}_metrics.txt --VALIDATION_STRINGENCY LENIENT --REMOVE_DUPLICATES true --VERBOSITY WARNING --TMP_DIR $picard_tmp

done

rm -r $picard_tmp
