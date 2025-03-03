#!/usr/bin/env python3
'''
Script designed to estimate CNV of MYO genes for the manuscript "A Segregating Structural Variant Defines Novel Venom Phenotypes in the Eastern Diamondback Rattlesnake"
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''

import os
import pandas as pd
import matplotlib.pyplot as plt

def _ParseCNV_(depthIN, bedIN, REFgene):
    depth = pd.read_csv(depthIN,sep='\t',names=['chr','position','depth'])
    Mdepth = int(depth['depth'].mean()) # get mean of depth of all regions
    STDdepth = depth['depth'].std() # get standard deviation of depth of all regions

    bed = pd.read_csv(bedIN,sep='\t',names=['chr','start','end','gene','score','strand'])

    Tdepthtemp = {}
    Tdepth = {}
    for index, row in bed.iterrows():
        subdepth = depth[depth['chr'] == row['chr']]
        subdepth = subdepth[subdepth['position'].isin(range(row['start'],row['end']))]
        norm = subdepth['depth'].mean()-Mdepth/STDdepth
        if row['gene'] in Tdepthtemp.keys():
            Tdepthtemp[row['gene']] = (Tdepthtemp[row['gene']] + norm ) / 2
        if not row['gene'] in Tdepthtemp.keys():
            Tdepthtemp[row['gene']] = norm
    AverageSL = sum([Tdepthtemp[x] for x in REFgene])/len(REFgene)
    minSL = (min([Tdepthtemp[x] for x in REFgene])/AverageSL)-0.1
    maxSL = (max([Tdepthtemp[x] for x in REFgene])/AverageSL)+0.1
    for k in Tdepthtemp.keys():
        Tdepth[k] = Tdepthtemp[k]/AverageSL
        if Tdepthtemp[k]/AverageSL < 0:
            Tdepth[k] = 0.01
        if str(Tdepthtemp[k]/AverageSL) == "nan":
            Tdepth[k] = 0.01

    COLOR = {"MYO-1":"orange","MYO-2":"orange","MYO-3":"orange","MYO-4":"orange"}
    for g in REFgene:
        COLOR[g] = "grey"
    FAM = [COLOR[x] for x in list(Tdepth.keys())]

    plt.scatter(list(Tdepth.keys()),list(Tdepth.values()), c=FAM)
    plt.axhline(y = 1, color = 'r', linestyle = '-')
    plt.axhline(y = minSL, color = 'r', linestyle = '--')
    plt.axhline(y = maxSL, color = 'r', linestyle = '--')
    plt.xticks(rotation=90)
    plt.ylim((0,12))
    plt.savefig(depthIN.replace(".txt",".adjusted.png"), bbox_inches='tight')
    plt.close()
    return Tdepth

path = "/path/to/working/dir/"
bedIN = path+"MYO_region.bed"
output = path+"depth_summary"
REFgene = ["ATPSynLipid-1", "ATPase-lys70", "CD63", "Calreticulin_", "DAZ-2", "GADD45", "Glutaredoxin1", "Leptin1", "PDI_e427_-_Exon", "Nexin2"]

FINAL = {}
for i in os.listdir(path):
    if i.endswith(".depth.txt"):
        s = i.replace(".depth.txt","")
        depthIN = path+i
        Tdepth = _ParseCNV_(depthIN, bedIN, REFgene)
        FINAL["depth"] = []
        for k in Tdepth.keys():
            FINAL["depth"].append(k)
            FINAL.setdefault(s,[])
            FINAL[s].append(Tdepth[k])

OUT = open(output,"w")
for k in FINAL.keys():
    OUT.write(k+"\t"+"\t".join([str(x) for x in FINAL[k]])+"\n")
OUT.close()

#END
