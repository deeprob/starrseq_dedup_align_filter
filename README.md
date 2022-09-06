# STARRSeq deduplicate, align and filter pipeline

STARRSeq sequencing data was cleaned to get a robust set of reads per region tested.

# STEPS

1. Remove PCR duplicates using *STARRDust*, an in-house tool.

2. Align deduplicated reads to the reference genome using *BWA MEM*.

3. Only keep reads that align to the reference regions and filter reads which do not pass quality control criterion using *SAMTOOLS*.

# Output

Filtered reads aligned to regions of interest

# Script descriptions

1. root/src/0_dedup_align_filter.py: Runs the deduplicate, align and filter pipeline.
