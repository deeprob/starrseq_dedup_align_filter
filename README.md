# STARRSeq analysis - deduplicate, align and filter pipeline

# Pipeline steps

## 1. Deduplication
**Definition:** Removal of PCR duplicates in the library

**Tools supported:** *starrdust* or *picard*

## 2. Alignment
**Definition:** Align reads to the reference genome

**Tools supported:** *BWA MEM*

## 3. Filter

**Definition:** Filter reads that do not pass certain quality control criterion

**Tools supported:** *samtools*

# How to run the pipeline

- Create metadata file
- Create conda environment
- Run 0_dedup_align_filter.py

## Create metadata file
Metadata file is a json file containing library specific information along with path to genome fasta file and region of interest bed file (if available). An example metadata file is present in *example_metadata* folder.

## Conda environment creation
```bash
foo@bar:~$ conda create -n starrseq -c conda-forge -c bioconda -c anaconda python=3.9 pyfastx bwa samtools picard
```

## Run 0_dedup_align_filter.py
*0_dedup_align_filter.py* script requires the following arguments. 
1. The path to the metadata file.
2. The starrseq library name as provided in the metadata file which needs to be analyzed.
3. The path to the raw directory where raw fastq files obtained from the sequencer is stored. 
4. The path to the directory where analyzed files will be stored. Approprite sub directories will be created automatically.

Example slurm scripts to run dedup align filter is given in *slurm* folder.  
