import pandas as pd
import numpy as np
from scipy import optimize
import pathlib
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import os
import argparse
import sys
from helper_functions import *

def _parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract reads from TSMA panel regions from a BAM file and save per-region CSVs."
    )
    parser.add_argument("-b", "--bamfile", required=True,
                        help="Path to input BAM file (sorted/indexed).")
    parser.add_argument("-o", "--outputdir", default="./output",
                        help="Base output directory (will create per-sample subdir).")
    return parser.parse_args(argv)

if __name__ == "__main__":
    args = _parse_args()
else:
    args = _parse_args([])

bamfile = args.bamfile
outputdir = args.outputdir

sampleid = bamfile.split("/")[-1].split(".")[0]
outputdir = os.path.join(outputdir, sampleid)
os.system(f"mkdir -p {outputdir}")

tsma_regiondf = pd.read_excel("./assets/TSMA panel.xlsx", skiprows=2)
tsma_regiondf.columns = ["chrom", "start", "end", "annotation", "panel"]
tsma_regiondf["panel"] = tsma_regiondf["panel"].apply(lambda x: x.replace(" ", "_"))
tsma_regiondf["regionname"] = tsma_regiondf[["chrom", "start", "end"]].apply(lambda x: f"{x[0]}_{x[1]}_{x[2]}", axis = 1)

all_panels = tsma_regiondf.panel.unique()

for panel in all_panels:
    all_regions = tsma_regiondf[tsma_regiondf.panel == panel].regionname.unique()
    os.system(f"mkdir -p {os.path.join(outputdir, panel)}")

    for regionname in tqdm(all_regions):
        if os.path.isfile(os.path.join(outputdir, panel, "readdf_{}.csv".format(regionname))) == False:
            region = f"chr{regionname.split('_')[0]}:{regionname.split('_')[1]}-{regionname.split('_')[2]}"

            sampleid = bamfile.split("/")[-1].split(".")[0]

            # ***** collect reads from regions from the TSMA
            bamfile_obj = pysam.AlignmentFile(bamfile).fetch(region = region)
            reads = []
            for read in bamfile_obj:
                reads.append(read)
            readdf = pd.DataFrame()
            readdf["chrom"] = [read.to_dict()["ref_name"] for read in reads]
            readdf["start"] = [read.to_dict()["ref_pos"] for read in reads]
            readdf["cigar"] = [read.to_dict()["cigar"] for read in reads]
            readdf["flen"] = [read.to_dict()["length"] for read in reads]
            readdf["seq"] = [read.to_dict()["seq"] for read in reads]
            readdf["methyl_string"] = [read.to_dict()["tags"][2] for read in reads]
            readdf["XR"] = [read.to_dict()["tags"][3] for read in reads]
            readdf["XG"] = [read.to_dict()["tags"][4] for read in reads]
            readdf["sample"] = sampleid
            readdf["region"] = region

            readdf.to_csv(os.path.join(outputdir, panel, "readdf_{}.csv".format(regionname)))