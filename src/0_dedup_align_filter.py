import os
import argparse
import utils as ut


def main(
    raw_dir, library_short,
    library_prefix, library_replicates, library_read_pairs, library_umi_index, library_suffix,
    reference_genome, roi_file, store_dir,
    deduplication_tool,
    dedup_flag=True, align_flag=True, filter_flag=True
    ):
    """
    Deduplicate, align and filter starrseq reads
    """
    ### initialize running prefix and suffix:: running prefix has dir struture along with filename prefix ###
    library_running_prefix, library_running_suffix = os.path.join(raw_dir, library_short, library_prefix), library_suffix
    
    ### Pre alignment deduplication with starrdust ###
    if deduplication_tool == "starrdust":
        # check to make sure umi index is given
        assert library_umi_index != ""
        # get deduped prefix and suffix
        library_deduped_prefix, library_deduped_suffix = ut.get_analyzed_filename_prefix_suffix(store_dir, library_short, library_prefix, "deduped_starrdust")
        # deduplicate
        if dedup_flag:
            # use starrdust to remove duplicates
            ut.remove_dups_with_starrdust(
                library_running_prefix, 
                library_replicates, 
                library_read_pairs, 
                library_umi_index, 
                library_running_suffix, 
                library_deduped_prefix, 
                library_deduped_suffix
                )
        # modify running prefix and suffix
        library_running_prefix, library_running_suffix = library_deduped_prefix, library_deduped_suffix
    
    ### Alignment ###
    # get aligned prefix and suffix
    library_aligned_prefix, library_aligned_suffix = ut.get_analyzed_filename_prefix_suffix(store_dir, library_short, library_prefix, "aligned")
    # align
    if align_flag:
        ut.align_reads(
            library_running_prefix, 
            library_replicates, 
            library_read_pairs, 
            library_running_suffix, 
            library_aligned_prefix,
            reference_genome)
    # modify running prefix and suffix
    library_running_prefix, library_running_suffix = library_aligned_prefix, library_aligned_suffix

    ### Post alignment deduplication with picard ###
    if deduplication_tool == "picard":
        # get deduped prefix and suffix
        library_deduped_prefix, library_deduped_suffix = ut.get_analyzed_filename_prefix_suffix(store_dir, library_short, library_prefix, "deduped_picard")
        # deduplicate
        if dedup_flag:
            # use starrdust to remove duplicates
            ut.remove_dups_with_picard(
                library_running_prefix, 
                library_replicates, 
                library_deduped_prefix
                )
        # modify running prefix and suffix
        library_running_prefix, library_running_suffix = library_deduped_prefix, library_deduped_suffix    

    ### Filtration ###
    # get filtered prefix and suffix
    library_filtered_prefix, library_filtered_suffix = ut.get_analyzed_filename_prefix_suffix(store_dir, library_short, library_prefix, "filtered")
    # filter
    if filter_flag:
        ut.filter_reads(
            library_running_prefix, 
            library_replicates, 
            library_filtered_prefix,
            roi_file
            )
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq analysis")
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored, look at $prepare meta$ package")
    parser.add_argument("lib", type=str, help="The library name as given in the meta file, this will run dedup, align and filter")
    parser.add_argument("raw_dir", type=str, help="Dir where raw data is stored")
    parser.add_argument("store_dir", type=str, help="Dir where analyzed data will be stored")
    parser.add_argument("-d", "--deduptool", type=str, help="deduplication tool to use; options starrdust/picard/None", default="")
    parser.add_argument("-p", "--dedup", action="store_false", help="Do not run the deduplication pipeline")
    parser.add_argument("-a", "--align", action="store_false", help="Do not run the alignment pipeline")
    parser.add_argument("-f", "--filter", action="store_false", help="Do not run the filter pipeline")
    parser.add_argument("-t", "--threads", type=int, help="number of threads to use", default=64) # TODO: include threads in main argument

    cli_args = parser.parse_args()
    lib_args = ut.create_args(cli_args.meta_file, cli_args.lib)
    main(
        cli_args.raw_dir,
        lib_args.library_short,
        lib_args.library_prefix, 
        lib_args.library_reps, 
        lib_args.library_pair,
        lib_args.library_umi,
        lib_args.library_suffix,
        lib_args.reference_genome,
        lib_args.roi_file,
        cli_args.store_dir,
        cli_args.deduptool,
        cli_args.dedup,
        cli_args.align,
        cli_args.filter,
        )
