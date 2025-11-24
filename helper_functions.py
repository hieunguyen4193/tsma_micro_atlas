import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pathlib
from tqdm import tqdm
import re
pattern = re.compile("^[1-9][0-9]M")
import warnings
warnings.filterwarnings("ignore")
import pyfaidx

def get_refseq(path_to_all_fa, chrom, start, end):
    
    refseq = pyfaidx.Fasta(os.path.join(path_to_all_fa, "{}.fa".format(chrom)))
    return(str.upper(refseq.get_seq(name = "{}".format(chrom), start = start, end = end).seq))
            
def get_CpG_status(read_start, read_end, read, cpg_pos, mode = "string"):    
    if (read_start <= cpg_pos) and (read_end >= cpg_pos):
        seq = read[cpg_pos - read_start: cpg_pos + 1 - read_start + 1]
    else:
        seq = "not_covered"
    if mode == "string":
        return seq
    elif mode == "num":
        if seq == "not_covered":
            seq = -1
        elif seq == "CG":
            seq = 1
        else:
            seq = 0
        return seq

