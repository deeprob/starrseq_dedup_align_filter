import os
import subprocess
import json
from argparse import Namespace
from itertools import starmap
import starrdust.starrdust as sd


#### GLOBALS ####
CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))

###############################
# read meta file; create args #
###############################

def create_args(meta_file, lib_name):
    with open(meta_file, "r") as f: 
        meta_dict = json.load(f)
        
    args = Namespace(
        # from metadata file
        library_prefix = meta_dict[lib_name]["prefix"],
        library_reps = meta_dict[lib_name]["replicates"],
        library_pair= meta_dict[lib_name]["read_pairs"],
        library_umi = meta_dict[lib_name]["umi"],
        library_suffix = meta_dict[lib_name]["suffix"],
        library_short = meta_dict[lib_name]["shortform"],
        reference_genome = meta_dict["genome"]["ref_fasta"],
        reference_genome_twobit = meta_dict["genome"]["ref_twobit"],
        roi_file = meta_dict["roi"]["sorted"]
    )

    return args

####################
# format filenames #
####################

def get_analyzed_filename_prefix_suffix(store_dir, lib_short, lib_prefix, analysis):
    lib_prefix_formatted = os.path.join(store_dir, analysis, lib_short, lib_prefix)
    suffix_dict = {"deduped": ".fastq", "aligned": ".bam", "filtered": ".bam"}
    lib_suffix_formatted = suffix_dict[analysis]
    return lib_prefix_formatted, lib_suffix_formatted

#########################
# deduplication helpers #
#########################

def get_deduped_files_helper(read1, read2, umi, read1_out, read2_out, mapq_r1=30, mapq_r2=30, min_length_r1=100, min_length_r2=100):
    sd.starrdust(read1, read2, umi, mapq_r1, mapq_r2, min_length_r1, min_length_r2, None, read1_out, read2_out)
    return

def get_deduped_files(lib_pre, lib_rep, lib_pairs, lib_umi_idx, lib_suff, outfile_pre, outfile_suf):
    lib_rep_list = lib_rep.split()
    lib_pair_list = lib_pairs.split()

    def get_read_file(pre, rep, pair, suff):
        fpath = "_".join([f for f in [pre, rep, pair] if f])
        return fpath + suff
    
    for rep in lib_rep_list:
        lib_read1_path, lib_read2_path, lib_umi_path = starmap(get_read_file, [(lib_pre, rep, i, lib_suff) for i in [lib_pair_list[0], lib_pair_list[1], lib_umi_idx]])
        lib_read1_outpath, lib_read2_outpath = starmap(get_read_file, [(outfile_pre, rep, i, outfile_suf) for i in [lib_pair_list[0], lib_pair_list[1]]])
        os.makedirs(os.path.dirname(outfile_pre), exist_ok=True)
        get_deduped_files_helper(lib_read1_path, lib_read2_path, lib_umi_path, lib_read1_outpath, lib_read2_outpath)
    return

def remove_dups(
    library_prefix, library_replicates, library_read_pairs, library_umi_index, library_suffix, library_deduped_prefix, library_deduped_suffix
    ):
    """
    Remove duplicates using starrdust
    """
    get_deduped_files(library_prefix, library_replicates, library_read_pairs, library_umi_index, library_suffix, library_deduped_prefix, library_deduped_suffix)
    return
 
#####################
# alignment helpers #
#####################

def align_reads_helper(reference_genome, library_prefix, library_replicates, 
                       library_read_pairs, library_suffix, library_aligned_prefix, 
                       cores=64):
    """
    Align reads using bwa - bash script under the ./shell_scripts dir
    """
    os.makedirs(os.path.dirname(library_aligned_prefix), exist_ok=True)
    cmd = [
        "bash", f"{CURRENT_DIR_PATH}/shell_scripts/1_align_fastq_to_ref.sh", 
        "-g", f"{reference_genome}", "-i", f"{library_prefix}", 
        "-r", f"{library_replicates}", "-p", f"{library_read_pairs}", 
        "-s", f"{library_suffix}", "-o", f"{library_aligned_prefix}",
         "-t", f"{cores}"]

    subprocess.run(cmd)
    return

def align_reads(
    library_prefix, library_replicates, library_read_pairs, library_suffix, library_aligned_prefix,
    reference_genome, cores=64
    ):
    """Align reads to the reference genome"""
    # align input reads to the reference genome
    align_reads_helper(reference_genome, library_prefix, library_replicates, library_read_pairs, library_suffix, library_aligned_prefix, cores)
    return

##################
# filter helpers #
##################

def filter_reads_helper(library_prefix, library_replicates, library_filtered_prefix, 
                        roi_file, dedup_flag):
    """
    Filter reads using samtools - bash script under the ./shell_scripts dir
    """
    os.makedirs(os.path.dirname(library_filtered_prefix), exist_ok=True)
    cmd = [
        "bash", f"{CURRENT_DIR_PATH}/shell_scripts/2_filter_reads.sh", 
        "-i", f"{library_prefix}", "-r", f"{library_replicates}", 
        "-o", f"{library_filtered_prefix}", "-f", f"{roi_file}", 
        "-d", f"{dedup_flag}"]

    subprocess.run(cmd)
    return

def filter_reads(library_aligned_prefix, library_replicates, library_filtered_prefix,
                roi_file, umi):
    """filter bad reads"""
    dedup = "true"
    if umi:
        dedup = "false"
    # filter input reads samtools -F 2828 or 2852 (depending on the use of starrdust or picard), -f 2 -q 30
    filter_reads_helper(library_aligned_prefix, library_replicates, library_filtered_prefix, roi_file, dedup)
    return
