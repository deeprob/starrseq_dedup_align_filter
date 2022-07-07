import sys
import os

import contextlib as cl
import pickle     as pk

import pyfastx as pf
import numpy   as np


def coroutine(func):
    def online(*args, **kwargs):
        _coroutine = func(*args, **kwargs)
        next(_coroutine)
        return _coroutine
    return online

@cl.contextmanager
def ignored(*exceptions):
    try:
        yield
    except exceptions:
        pass

@coroutine
def liner_engine(online=True):
    clrlen = 0
    try:
        while True:
            printstr = (yield)
            if online and printstr.startswith('\*'):
                sys.stdout.write('\n')
                clrlen   = 0
                printstr = printstr.removeprefix('\*')
            if online:
                sys.stdout.write('\r' + ' '*clrlen)
                sys.stdout.write('\r' + printstr)
                clrlen = len(printstr)
                sys.stdout.flush()
    except GeneratorExit:
        if online:
            sys.stdout.write('\n')

def get_outfq_fname(infile):
    outfname = infile.removesuffix('.gz')
    outfname = outfname.removesuffix('.fastq')
    outfname = outfname.removesuffix('.fq')
    return outfname + '.deduped.fastq'

def get_adjusted_path(path, suffix):
    if not isinstance(path, str):
        return None
    path = path.strip()
    path = path.removesuffix('/')
    if not suffix is None and \
       not path.endswith(suffix):
        path += str(suffix)
    return path

def remove_file(filepath):
    if not filepath is None:
        filepath = get_adjusted_path(
            path=filepath,
            suffix=None)
        with ignored(OSError):
            os.remove(filepath)

def picklesave(obj, filepath):
    if os.path.exists(filepath):
        raise RuntimeError(
            '{} exists!'.format(filepath))
    with open(filepath, 'wb') as outfile:
        pk.dump(
            obj=obj,
            file=outfile,
            protocol=5,
            fix_imports=True)

def pickleload(filepath):
    with open(filepath, 'rb') as infile:
        obj = pk.load(
            file=infile,
            fix_imports=True)
    return obj

def safediv(A, B):
    return 0. if B == 0. else float(A) / B

def safelog10(A):
    return np.log10(A) if A > 0. else 0.

def get_printlen(value):
    svalue = str(value)
    if len(svalue) <= 3:
        return len(svalue)
    else:
        return len(str(value)) + \
            int(safelog10(value) / 3)

def stream_fastq_engine(
    filepath):

    try:
        ptr = iter(pf.Fastq(
            file_name=filepath,
            build_index=False,
            full_index=False))

        rcount = 0
        tcount = 0

        while True:

            try:
                if rcount < tcount:
                    next(ptr)
                    rcount += 1
                    continue
                fqe  = next(ptr)
                head = fqe[0]
                read = fqe[1]
                qual = fqe[2]
                qvec = np.frombuffer(
                    qual.encode(),
                    dtype=np.int8) - 33
                yield (head, read, qual, qvec)
                rcount += 1
                tcount += 1
            except Exception as E:
                break

    except Exception as E:
        raise E

    yield None, None, None, None