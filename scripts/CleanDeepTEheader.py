#!/usr/bin/env python3
'''
Clean up headers output from deepTE to run RepatMasker and annotate TEs and plot RepeatLandscape
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

usage:
CleanDeepTEheader.py input_lib.fasta output_lib_clean.fasta

Putative Header from DeepTE to clean:
ClassI
ClassII
ClassIII_Helitron
ClassII_DNA_CACTA_MITE
ClassII_DNA_CACTA_nMITE
ClassII_DNA_Harbinger_MITE
ClassII_DNA_Harbinger_nMITE
ClassII_DNA_Harbinger_unknown
ClassII_DNA_Mutator_MITE
ClassII_DNA_Mutator_nMITE
ClassII_DNA_Mutator_unknown
ClassII_DNA_PiggyBac_MITE
ClassII_DNA_PiggyBac_nMITE
ClassII_DNA_TcMar_MITE
ClassII_DNA_TcMar_nMITE
ClassII_DNA_TcMar_unknown
ClassII_DNA_hAT_MITE
ClassII_DNA_hAT_nMITE
ClassII_DNA_hAT_unknown
ClassII_MITE
ClassII_nMITE
ClassI_LTR
ClassI_LTR_BEL
ClassI_LTR_Copia
ClassI_LTR_ERV
ClassI_LTR_Gypsy
ClassI_nLTR
ClassI_nLTR_DIRS
ClassI_nLTR_LINE
ClassI_nLTR_LINE_I
ClassI_nLTR_LINE_Jockey
ClassI_nLTR_LINE_L1
ClassI_nLTR_LINE_R2
ClassI_nLTR_LINE_RTE
ClassI_nLTR_PLE
ClassI_nLTR_SINE_tRNA
unknown
'''

import sys
from Bio import SeqIO

def _Run_(TEin, TEout):
    TElib = SeqIO.to_dict(SeqIO.parse(TEin, 'fasta'))
    TEOUT = open(TEout,"w")
    written = []
    for k in TElib:
        if k.split("__")[-1].startswith("ClassI_"):
            if "_LTR" in k:
                id = k.split("#")[0]+"#LTR/"
                if len(k.split("__")[-1].split("_")) == 2:
                    id = id+"Unknown"
                    TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
                if len(k.split("__")[-1].split("_")) > 2:
                    id = id+k.split("__")[-1].split("_")[2]
                    TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
            if "_nLTR" in k:
                if "_LINE" in k:
                    id = k.split("#")[0]+"#LINE/"
                    if len(k.split("__")[-1].split("_")) == 3:
                        id = id+"Unknown"
                        TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
                    if len(k.split("__")[-1].split("_")) > 3:
                        id = id+k.split("__")[-1].split("_")[3]
                        TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
                if "_SINE" in k:
                    id = k.split("#")[0]+"#SINE/"
                    if len(k.split("__")[-1].split("_")) == 3:
                        id = id+"Unknown"
                        TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
                    if len(k.split("__")[-1].split("_")) > 3:
                        id = id+k.split("__")[-1].split("_")[3]
                        TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
                if "_DIRS" in k:
                    id = k.split("#")[0]+"#LTR/DIRS"
                    TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
                if "_PLE" in k:
                    id = k.split("#")[0]+"#Penelope/Penelope"
                    TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
        if k.split("__")[-1] == "ClassI":
            id = k.split("#")[0]+"#LTR/Unknown"
            TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
        if k.split("__")[-1].startswith("ClassII_"):
            if "_MITE" in k:
                id = k.split("#")[0]+"#MITE/"
                if len(k.split("__")[-1].split("_")) == 2:
                    id = id+"Unknown"
                    TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
                if len(k.split("__")[-1].split("_")) > 2:
                    id = id+k.split("__")[-1].split("_")[2]
                    TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
            if "_nMITE" in k:
                id = k.split("#")[0]+"#DNA/"
                if len(k.split("__")[-1].split("_")) == 2:
                    id = id+"Unknown"
                    TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
                if len(k.split("__")[-1].split("_")) > 2:
                    id = id+k.split("__")[-1].split("_")[2]
                    TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
        if k.split("__")[-1] == "ClassII":
            id = k.split("#")[0]+"#DNA/Unknown"
            TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
        if k.split("__")[-1].startswith("ClassIII_"):
            id = k.split("#")[0]+"#DNA/"+k.split("ClassIII_")[-1]
            TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
        if k.split("__")[-1].startswith("unknown"):
            id = k.split("#")[0]+"#Unknown"
            TEOUT.write(">"+id+"\n"+str(TElib[k].seq)+"\n")
    TEOUT.close()


def _main_():
    if len (sys.argv) != 3:
        print("Script designed to clean header from DeepTE output")
        print("Basic usage: CleanDeepTEheader.py input_lib.fasta clean_lib.fasta")
        quit()

    TEin = sys.argv[1]
    TEout = sys.argv[2]
    _Run_(TEin, TEout)

_main_()

#END
