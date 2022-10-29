# STARRSeq analysis - deduplicate, align and filter pipeline

# Pipeline steps

## 1. Deduplication
**Definition:** Removal of PCR duplicates in the library

**Tools supported:** *starrdust* or *picard*

**Description:** 

## 2. Alignment
**Definition:** Align reads to the reference genome

**Tools supported:** *BWA MEM*

**Description:** 

## 3. Filter

**Definition:** Filter reads that do not pass certain quality control criterion using *SAMTOOLS*

**Tools supported:** *samtools*

**Description:** 

# Required tools

# Conda environment creation

# How to run the pipeline

1. Create metadata file
2. Create conda environment
3. Run 0_dedup_align_filter.py
