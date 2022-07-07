import starrdust as sd

def main():
    base  = '../Reads/'
    base += 'F1_SCRT1/'

    sd.starrdust(
        r1_input_file=base + 'F1_SCRT1_R2_S4_R1_001.fastq.gz',
        r2_input_file=base + 'F1_SCRT1_R2_S4_R3_001.fastq.gz',
        umi_input_file=base + 'F1_SCRT1_R2_S4_R2_001.fastq.gz',
        r1_qscore_min=30,
        r2_qscore_min=30,
        r1_readlen_min=140,
        r2_readlen_min=140,
        dedup_file=None)

if __name__ == '__main__':
    main()