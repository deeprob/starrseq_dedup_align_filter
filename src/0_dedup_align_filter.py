import os
import argparse
import utils as ut


def main(
    library_prefix, library_replicates, library_read_pairs, library_umi_index, library_suffix,
    reference_genome, roi_file,
    dedup_flag=True, align_flag=True, filter_flag=True
    ):
    """
    Deduplicate, align and filter starrseq reads
    """
    # initialize running prefix and suffix
    library_running_prefix, library_running_suffix = library_prefix, library_suffix
    
    # Deduplication
    # deduplicate only if there is an umi flag
    umi_flag = False
    if library_umi_index:
        umi_flag = True
        # get deduped prefix and suffix
        library_deduped_prefix, library_deduped_suffix = ut.get_analyzed_filename_prefix_suffix(library_prefix, "deduped")
        # deduplicate
        if dedup_flag:
            # use starrdust to remove duplicates
            ut.remove_dups(
                library_running_prefix, library_replicates, library_read_pairs, library_umi_index, library_running_suffix, library_deduped_prefix, library_deduped_suffix)
            # modify running prefix and suffix
            library_running_prefix, library_running_suffix = library_deduped_prefix, library_deduped_suffix
    
    # Alignment
    # get aligned prefix and suffix
    library_aligned_prefix, library_aligned_suffix = ut.get_analyzed_filename_prefix_suffix(library_prefix, "aligned")
    # align
    if align_flag:
        ut.align_reads(
            library_running_prefix, library_replicates, library_read_pairs, library_running_suffix, library_aligned_prefix,
            reference_genome)
    # modify running prefix and suffix
    library_running_prefix, library_running_suffix = library_aligned_prefix, library_aligned_suffix
    
    # Filtration
    # get filtered prefix and suffix
    library_filtered_prefix, library_filtered_suffix = ut.get_analyzed_filename_prefix_suffix(library_prefix, "filtered")
    # filter
    if filter_flag:
        ut.filter_reads(
            library_running_prefix, library_replicates, library_filtered_prefix,
            roi_file, 
            umi=umi_flag)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq analysis")
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored, look at $prepare meta$ package")
    parser.add_argument("lib", type=str, help="The library name as given in the meta file, this will run dedup, align and filter")
    parser.add_argument("-d", "--dedup", action="store_false", help="Do not run the deduplication pipeline")
    parser.add_argument("-a", "--align", action="store_false", help="Do not run the alignment pipeline")
    parser.add_argument("-f", "--filter", action="store_false", help="Do not run the filter pipeline")
    parser.add_argument("-t", "--threads", type=int, help="number of threads to use", default=64) # TODO: include threads in main argument

    cli_args = parser.parse_args()
    lib_args = ut.create_args(cli_args.meta_file, cli_args.lib)
    main(
        lib_args.library_prefix, 
        lib_args.library_reps, 
        lib_args.library_pair,
        lib_args.library_umi,
        lib_args.library_suffix,
        lib_args.reference_genome,
        lib_args.roi_file,
        cli_args.dedup,
        cli_args.align,
        cli_args.filter
        )
