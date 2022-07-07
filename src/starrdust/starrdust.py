
import time    as tt
import numpy   as np
from . import utils   as ut
import hashlib as hl


def starrdust(
    r1_input_file,
    r2_input_file,
    umi_input_file,
    r1_qscore_min,
    r2_qscore_min,
    r1_readlen_min,
    r2_readlen_min,
    dedup_file,
    r1_outfname="",
    r2_outfname=""
    ):

    liner = ut.liner_engine(online=True)
    liner.send('\n[STARRdust - STARR-Seq Read Deduplicator]\n')
    liner.send('\n [Experiment Information]\n')
    liner.send('   R1 Input File  : {}\n'.format(r1_input_file))
    liner.send('   R2 Input File  : {}\n'.format(r2_input_file))
    liner.send('  UMI Input File  : {}\n'.format(umi_input_file))
    liner.send('   R1   Min QScore: {}\n'.format(r1_qscore_min))
    liner.send('   R2   Min QScore: {}\n'.format(r2_qscore_min))
    liner.send('   R1   Min Length: {}\n'.format(r1_readlen_min))
    liner.send('   R2   Min Length: {}\n'.format(r2_readlen_min))
    liner.send('  Dedup Dictionary: {}\n'.format(dedup_file))
    total_reads = 0
    bad_reads   = 0
    low_reads   = 0
    short_reads = 0
    final_reads = 0
    dedup_reads = 0
    r1_reader  = ut.stream_fastq_engine(
        filepath=r1_input_file)
    r2_reader  = ut.stream_fastq_engine(
        filepath=r2_input_file)
    umi_reader = ut.stream_fastq_engine(
        filepath=umi_input_file)
    r1_outfname = ut.get_outfq_fname(infile=r1_input_file) if not r1_outfname else r1_outfname
    r2_outfname = ut.get_outfq_fname(infile=r2_input_file) if not r2_outfname else r2_outfname

    ut.remove_file(filepath=r1_outfname)
    ut.remove_file(filepath=r2_outfname)

    r1_writer = open(r1_outfname, 'w')
    r2_writer = open(r2_outfname, 'w')
    q1status = r1_qscore_min > 0
    q2status = r2_qscore_min > 0
    mcharset = set('ATGC')
    rstatus  = np.random.randint(0.5*10**4, 10**4)
    if dedup_file is None:
        dedup_data = {}
    else:
        dedup_data = ut.pickleload(filepath=dedup_file)
    t0 = tt.time()
    liner.send('\n [Starting Deduplication]\n')

    while True:
        h1,r1,q1,qv1  = next(r1_reader)
        h2,r2,q2,qv2  = next(r2_reader)
        _,umi,__,qvu = next(umi_reader)
        reader_status = int(r1 is None) + int(r2 is None) + int(umi is None)
        if 0 < reader_status < 3:
            raise RuntimeError('Truncated Split')
        if reader_status == 3:
            break
        total_reads += 1
        # if (total_reads) % rstatus == 0:
        #     liner.send('   Processed: {:,} Reads'.format(total_reads))
        if ('N' in r1) or ('N' in r2) or ('N' in umi):
            bad_reads += 1
            continue
        charset = set(umi)
        if not (charset <= mcharset):
            bad_reads += 1
            continue
        q1s = np.mean(qv1) >= r1_qscore_min if q1status else True
        q2s = np.mean(qv2) >= r2_qscore_min if q2status else True
        ums = np.mean(qvu) >= 30

        if not all((q1s, q2s, ums)):
            low_reads += 1
            continue
        if (len(r1) < r1_readlen_min) or \
           (len(r2) < r2_readlen_min):
            short_reads += 1
            continue
        final_reads += 1
        deduped = False
        if not umi in dedup_data:
            dedup_data[umi] = set()
        read = r1[:r1_readlen_min] + r2[:r2_readlen_min]
        read = hl.md5(read.encode()).hexdigest()
        if read in dedup_data[umi]:
            dedup_reads += 1
            deduped = True
        else:
            dedup_data[umi].add(read)
        if not deduped:
            r1_writer.write('@{}\n{}\n+\n{}\n'.format(
                h1, r1, q1))
            r2_writer.write('@{}\n{}\n+\n{}\n'.format(
                h2, r2, q2))

    liner.send('   Processed: {:,} Reads'.format(total_reads))
    liner.send('\*   Status: Completed\n')
    plen = ut.get_printlen(
        value=total_reads)

    liner.send('\n [Deduplication Stats]\n')
    liner.send(
        '   Total Reads: {:{},d}\n'.format(
            total_reads,
            plen))
    liner.send(
        '     Bad Reads: {:{},d} ({:6.2f} %)\n'.format(
            bad_reads,
            plen,
            ut.safediv(
                A=(100. * bad_reads),
                B=total_reads)))
    liner.send(
        '   Short Reads: {:{},d} ({:6.2f} %)\n'.format(
            short_reads,
            plen,
            ut.safediv(
                A=(100. * short_reads),
                B=total_reads)))
    liner.send(
        '   LowQS Reads: {:{},d} ({:6.2f} %)\n'.format(
            low_reads,
            plen,
            ut.safediv(
                A=(100. * low_reads),
                B=total_reads)))
    liner.send(
        '   Final Reads: {:{},d} ({:6.2f} %)\n'.format(
            final_reads,
            plen,
            ut.safediv(
                A=(100. * final_reads),
                B=total_reads)))
    liner.send(
        '   Dedup Reads: {:{},d} ({:6.2f} %)\n'.format(
            dedup_reads,
            plen,
            ut.safediv(
                A=(100. * dedup_reads),
                B=total_reads)))
    liner.send(
        '  Time Elapsed: {:.2f} sec\n'.format(
            tt.time() - t0))
    r1_writer.close()
    r2_writer.close()
    liner.close()
    return total_reads, bad_reads, short_reads,  low_reads, dedup_reads, final_reads
