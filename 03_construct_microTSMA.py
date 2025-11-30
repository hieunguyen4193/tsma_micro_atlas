import pandas as pd 
import numpy as np 
import seaborn as sns 
import matplotlib.pyplot as plt 
from tqdm import tqdm 
import warnings
warnings.filterwarnings("ignore")
import os
import pathlib
import random 
import argparse
from scipy import optimize


def _parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="construct micro TSMA at each region, perform LOOCV"
    )
    parser.add_argument("-r", "--region", required=True,
                        help="Region identifier (e.g., 10_104170538_104170638).")
    parser.add_argument("-l", "--loocv_val_sample", required=True,
                        help="Sample ID for LOOCV validation.")
    return parser.parse_args(argv)

if __name__ == "__main__":
    args = _parse_args()
else:
    args = _parse_args([])
    
inputdir = "/Volumes/HNWD02/outdir/output_methyl_ReadBased_Models/01_output/TSMA_panel"
path_to_all_fa = "/Volumes/HNSD01/storage/ref/hg19"
outdir = "/Volumes/HNSD01/outdir"
outputdir = os.path.join(outdir, "output_methyl_ReadBased_Models")
path_to_01_output = os.path.join(outputdir, "01_output")
path_to_02_output = os.path.join(outputdir, "02_output")

metadata = pd.read_excel("/Users/hieunguyen/src/tsma_micro_atlas/metadata/metadata_TSMA_R8281_R8282.xlsx")

# loocv_val_sample = random.sample(metadata.LABCODE.tolist(), 1)[0]

# region = "10_104170538_104170638"
# loocv_val_sample = "LAAL31TS"

region = args.region
loocv_val_sample = args.loocv_val_sample

paneldf = pd.read_excel("./assets/12967_2024_5416_MOESM1_ESM.xlsx")
all_regions = paneldf.Region_name.unique()

predictiondf = pd.DataFrame()

# for region in tqdm(all_regions):
path_to_03_output = os.path.join(outputdir, "03_output", 
                                    f"val_on_{loocv_val_sample}", 
                                    "microAtlas_data",
                                    f"region_{region}")
os.system(f"mkdir -p {path_to_03_output}")
df = pd.read_csv(os.path.join(path_to_02_output, f"{region}_methylString.csv"), index_col = [0])
if loocv_val_sample in df["LABCODE"].unique():   
    df_test = df[df["LABCODE"] == loocv_val_sample]
    df = df[df["LABCODE"] != loocv_val_sample]

    cpg_coords = [item for item in df.columns if item not in 
    ["methyl_string", "SampleID","cover","LABCODE","TYPE"] ]

    df = df[df["cover"] == df["cover"].max()]

    # Filter rows where all cpg_coord columns are either 0 or 1
    maskdf = pd.DataFrame()
    for c in cpg_coords:
        maskdf[c] = [df[df[c].isin([0,1])].shape[0]]
    maskdf = maskdf.T.reset_index()
    full_cover_cpgs = maskdf[maskdf[0] == df.shape[0]]["index"].tolist()
    df = df[full_cover_cpgs + ["methyl_string", "SampleID","cover","LABCODE","TYPE"]]
    df["full_cover_methyl_string"] = df[full_cover_cpgs].astype(str).agg(''.join, axis=1)

    # save a list of full covered CpG sites in the region
    pd.DataFrame({
        "CpG" : full_cover_cpgs
    }).to_csv(os.path.join(path_to_03_output, f"{region}_full_covered_CpGs.csv"), index=False)

    wbc_patterns = df[df["TYPE"] == "Control"].full_cover_methyl_string.unique()

    df_nowbc = df[df["full_cover_methyl_string"].isin(wbc_patterns) == False]
    if df_nowbc.shape[0] == 0:
        tmp_prediction = pd.DataFrame({
            "region": region,
            "loocv_val_sample": loocv_val_sample,
            "prediction": "WBC"
        }, index=[0])
        # predictiondf = pd.concat([predictiondf, tmp_prediction], axis=0)
        print(f"The region {region} only has WBC patterns after filtering, skipping...")
        # continue
    else:
        # ***** group by samples
        dfcount = df_nowbc.groupby(["SampleID", "full_cover_methyl_string"])["TYPE"].count().reset_index()
        dfcount_wide = dfcount.pivot(index='full_cover_methyl_string', columns='SampleID', values='TYPE').fillna(0)
        for n in dfcount_wide.columns:
            dfcount_wide[n] = dfcount_wide[n]/dfcount_wide[n].sum()

        plt.figure(figsize=(25,25))
        sns.heatmap(dfcount_wide)# ***** filter out WBC patterns
        plt.savefig(os.path.join(path_to_03_output, f"{region}_methyl_pattern_by_sample_heatmap.pdf"))
        plt.close()

        # ***** group by class
        toodf = df_nowbc.groupby(["TYPE", "full_cover_methyl_string"])["SampleID"].count().reset_index().\
            pivot(index='full_cover_methyl_string', columns='TYPE', values='SampleID').fillna(0)
        for c in toodf.columns:
            toodf[c] = toodf[c]/toodf[c].sum()

        toodf = toodf[(toodf != 0).any(axis=1)]
        toodf.to_csv(os.path.join(path_to_03_output, f"{region}_methyl_pattern_by_tissue_type.csv"), index=True, header=True)

        plt.figure(figsize=(12, 12))
        sns.heatmap(toodf, square=True, linewidths=0.5, linecolor='gray', cbar_kws={"shrink": 0.8})
        plt.tight_layout()
        plt.savefig(os.path.join(path_to_03_output, f"{region}_methyl_pattern_by_tissue_type_heatmap.pdf"))
        # plt.show()
        plt.close()

        # ***** calculate pairwise cosine similarity between samples
        from sklearn.metrics.pairwise import cosine_similarity

        # Calculate cosine similarity between all columns of toodf
        cosine_sim = cosine_similarity(toodf.T)

        # Create a DataFrame with proper labels
        cosine_sim_df = pd.DataFrame(cosine_sim, 
                                    index=toodf.columns, 
                                    columns=toodf.columns)

        # Visualize the cosine similarity matrix
        plt.figure(figsize=(8, 6))
        sns.heatmap(cosine_sim_df, annot=True, fmt='.3f', cmap='coolwarm', 
                    square=True, linewidths=0.5, cbar_kws={"shrink": 0.8})
        plt.title('Cosine Similarity between Tissue Types')
        plt.tight_layout()
        plt.savefig(os.path.join(path_to_03_output, f"{region}_cosine_similarity_heatmap.pdf"))
        plt.close()

        sum_cos_sim = np.sum(cosine_sim_df.sum().sort_values(ascending=False) - 1)
        pd.DataFrame({"region": region, "sum_cosine_similarity": sum_cos_sim}, index=[0]).to_csv(os.path.join(path_to_03_output, f"{region}_sum_cosine_similarity.csv"), index=False)

        # micro atlas prediction 
        df_test = df_test[full_cover_cpgs + ["methyl_string", "SampleID","cover","LABCODE","TYPE"]]
        df_test = df_test[df_test["cover"] == df_test["cover"].max()]
        df_test["full_cover_methyl_string"] = df_test[full_cover_cpgs].astype(str).agg(''.join, axis=1)
        df_test_nowbc = df_test[df_test["full_cover_methyl_string"].isin(wbc_patterns) == False]

        test_countdf = df_test_nowbc.groupby("full_cover_methyl_string")["SampleID"].count().reset_index()
        test_countdf[loocv_val_sample] = test_countdf["SampleID"]/test_countdf["SampleID"].sum()
        test_countdf = test_countdf.drop("SampleID", axis = 1)

        if test_countdf.shape[0] == 0:
            test_countdf = pd.DataFrame({
                "full_cover_methyl_string": toodf.index.tolist(),
                loocv_val_sample: 0
            })
            prediction = "WBC"
        else:
            atlasdf = toodf.merge(test_countdf, on="full_cover_methyl_string", how = "left").set_index("full_cover_methyl_string").fillna(0)
            if atlasdf[loocv_val_sample].sum() == 0:
                prediction = "WBC"
            else:
                mixture, residual = optimize.nnls(atlasdf[[item for item in atlasdf.columns if item != loocv_val_sample]].to_numpy(), 
                                            atlasdf[loocv_val_sample].to_numpy())
                mixture /= np.sum(mixture)
                mixturedf = pd.DataFrame({
                    "Tissue": [item for item in atlasdf.columns if item != loocv_val_sample],
                    "Proportion": mixture
                })
                prediction = mixturedf.loc[mixturedf['Proportion'].idxmax(), 'Tissue']
        tmp_prediction = pd.DataFrame({
            "region": region,
            "loocv_val_sample": loocv_val_sample,
            "prediction": prediction
        }, index=[0])
else:
    tmp_prediction = pd.DataFrame({
        "region": region,
        "loocv_val_sample": loocv_val_sample,
        "prediction": "no data"
    }, index=[0])
    print(f"The sample {loocv_val_sample} not found in region {region}, skipping...")
# predictiondf = pd.concat([predictiondf, tmp_prediction], axis=0)
tmp_prediction.to_excel(os.path.join(path_to_03_output, f"{region}_prediction_on_{loocv_val_sample}.xlsx"))
# predictiondf.to_excel(os.path.join(outputdir, "03_output", f"prediction_on_LOOCV_{loocv_val_sample}.xlsx"))

# parallel -j 20 python ./03_construct_microTSMA.py -r {} -l ${sample} ::: $regions