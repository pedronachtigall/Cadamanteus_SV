#!/usr/bin/env python3
'''
Script designed to calculate the coverage average for specific regions from a bed/bedGraph file.
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

inputs:
coverage.bed/bedGraph
regions.bed

output:
region_coverage.tsv

PS: This script was designed for the manuscript: "A Segregating Structural Variant Defines Novel Venom Phenotypes in the Eastern Diamondback Rattlesnake"
'''

import os
import pandas as pd

#parse regions.bed and generate a dictionary - name:chr:st-end
def _GetRegions_(region):
    final = {}
    df = pd.read_csv(region,sep='\t',names=['chr','st','end','name'])
    df.drop_duplicates(inplace=True)
    df.sort_values(by=['chr','st'], ascending=True, inplace=True)
    df.reset_index(drop=True,inplace=True)
    for index, row in df.iterrows():
        final[row['name']] = [row['chr'],row['st'],row['end']]
    return final

def _Run_(cov, region):
    R = _GetRegions_(region)

    #normalize bedgraph -> Zi = (Xi – min(X)) / (max(X) – min(X)) * 100
    COV = pd.read_csv(cov, sep='\t', names=['chr','st','end','cov','pvalue','peak'])
    COV["cov_norm"] = ( COV["cov"] - COV['cov'].min() ) / ( COV['cov'].max() - COV['cov'].min() ) * 100
    COV["peak_norm"] = ( COV["peak"] - COV['peak'].min() ) / ( COV['peak'].max() - COV['peak'].min() ) * 100
    COV['cov_norm'] = COV['cov_norm'].replace(0.0, 0.00001)
    COV['peak_norm'] = COV['peak_norm'].replace(0.0, 0.00001)

    FINAL = {}
    for k in list(R.keys()):#[1:11]:
        FINAL[k] = COV[COV['st'].between(R[k][1], R[k][2])]['cov_norm'].mean()

    return FINAL

path = "/path/to/working/dir/"
region = "SVMP_region.bed"
TSV = {}
for i in os.listdir(path):#[:4]:
    sample = i.split("_")[2]
    cov = path+i
    F = _Run_(cov, region)
    TSV["gene"] = []
    gene = []
    TSV.setdefault(sample, [])
    for k in F.keys():
        gene.append(k)
        TSV[sample].append(F[k])
    TSV["gene"] = gene

out = pd.DataFrame.from_dict(TSV).T

out.to_csv("/path/to/some/dir/SVMP_resolving.tsv", sep="\t")

#END
