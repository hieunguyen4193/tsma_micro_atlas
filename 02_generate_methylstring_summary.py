import argparse    
import sys
import pandas as pd 
import numpy as np 
import seaborn as sns 
import matplotlib.pyplot as plt 
from tqdm import tqdm 
import warnings
warnings.filterwarnings("ignore")
import os
import pathlib

def get_parser(argv):
    parser = argparse.ArgumentParser(
    description=(
        "Generate methylation string summary dataframe from read based methylation model output."
    )
    )
    parser.add_argument("--region",
                        type=str,
                        help="Name of the region to process.",
                        )
    return parser.parse_args(argv)

def main(args):
    region = args.region
    inputdir = "/Volumes/HNWD02/outdir/output_methyl_ReadBased_Models/01_output/TSMA_panel"
    path_to_all_fa = "/Volumes/HNSD01/storage/ref/hg19"
    outdir = "/Volumes/HNSD01/outdir"
    outputdir = os.path.join(outdir, "output_methyl_ReadBased_Models")
    path_to_01_output = os.path.join(outputdir, "01_output")
    path_to_02_output = os.path.join(outputdir, "02_output")
    os.system(f"mkdir -p {path_to_02_output}")

    all_regions = [item.name for item in pathlib.Path(inputdir).glob("*") if ".DS_Store" not in item.name]

    # for region in tqdm(all_regions):
    if os.path.isfile(os.path.join(path_to_02_output, f"{region}_methylString.csv")) == False:
        print(f"working on region: {region}")
        all_readdf = [item for item in pathlib.Path(inputdir).glob(f"{region}/*.tmpdf.csv")]
        df = pd.DataFrame()
        for input_file in tqdm(all_readdf):
                tmpdf = pd.read_csv(input_file, index_col = [0])
                cpg_coords = [item for item in tmpdf.columns if item not in ['chrom', 'start', 'cigar', 'flen', 'seq', 'methyl_string', 'XR', 'XG',
                        'sample', 'region', 'check_cigar', 'end']]
                tmpdf = tmpdf[cpg_coords]
                # tmpdf = tmpdf[(tmpdf != -1).all(axis=1)]
                tmpdf["methyl_string"] = tmpdf[cpg_coords].apply(lambda x: ''.join(x.astype(str).tolist()), axis = 1)
                # tmpdf = tmpdf.drop_duplicates(subset = ["methyl_string"])
                tmpdf["SampleID"] = input_file.name.replace(".tmpdf.csv", "")
                tmpdf["cover"] = tmpdf[cpg_coords].apply(lambda x: len([i for i in x if i != -1]), axis = 1)
                df = pd.concat([df, tmpdf], axis = 0)
        metadata = pd.read_excel("./metadata/metadata_TSMA_R8281_R8282.xlsx")[["LABCODE", "TYPE"]]
        df = df.merge(metadata, right_on = "LABCODE", left_on = "SampleID")
        df.to_csv(os.path.join(path_to_02_output, f"{region}_methylString.csv"))
    else:
        print(f"region {region} already done!")
if __name__ == '__main__':
    main(get_parser(sys.argv[1:]))


# example cmds: 
# parallel -j 20 python 02_generate_methylstring_summary.py --region {} ::: ${regions}