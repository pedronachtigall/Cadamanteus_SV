#!/usr/bin/env python3
'''
Script designed to adjust the GFF output by RepeatMasker (using the parameter -gff) based on TElib used as input.
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

inputs:
TElib.fa
RMout.gff

output:
RMout_adjusted.gff
'''

import sys

def _ParseTElib_(TElib):
    final = {}
    a = open(TElib,"r")
    for line in a:
        if line.startswith(">"):
            idseq = line.strip().replace(">","").split("#")[0]
            TEfam = line.strip().split("#")[1].split("/")[0]
            idfull = line.strip().replace(">","")
            final[idseq] = [TEfam, idfull]
    a.close()
    return final

def _Run_(GFFrm, TElib, GFFout):
    TE = _ParseTElib_(TElib)
    gff = open(GFFrm, "r")
    OUT = open(GFFout,"w")
    count = 1
    for line in gff:
        if not line.startswith("#") and not line.startswith("\n"):
            line1 = line.strip().split("\t")
            if line1[-1].startswith("Target "):
                idseq = line1[-1].replace("Target ","").split()[0].replace("\"","").replace("Motif:","").replace("-int","")
            if line1[-1].startswith("Target="):
                idseq = line1[-1].replace("Target=","").split()[0].replace("-int","")
            chr = line1[0]
            source = line1[1]
            st = line1[3]
            end = line1[4]
            score = line1[5]
            strand = line1[6]
            frame = line1[7]
            if idseq in TE.keys(): #present in the TElib
                feature = TE[idseq][0]
                desc = "ID=TE_"+str(count)+";Name="+TE[idseq][1].replace("#","__")+";Classification="+TE[idseq][1].split("#")[1]+";"
            if idseq not in TE.keys(): #not in TElib but annotated by RepeatMasker (simple, low_comp, satellite)
                if ")n" in idseq: #Simple_repeat
                    feature = "Simple_repeat"
                    desc = "ID=TE_"+str(count)+";Name="+idseq+"__Simple_repeat;Classification=Simple_repeat;"
                    if idseq == "(GAATG)n" or idseq == "(CATTC)n": #Satellite
                        feature = "Satellite"
                        desc = "ID=TE_"+str(count)+";Name="+idseq+"__Satellite;Classification=Satellite;"
                if "-rich" in idseq: #Low_complexity
                    feature = "Low_complexity"
                    desc = "ID=TE_"+str(count)+";Name="+idseq+"__Low_complexity;Classification=Low_complexity;"
            newline = [chr,source,feature,st,end,score,strand,frame,desc]
            OUT.write("\t".join(newline)+"\n")
            count += 1
    gff.close()
    OUT.close()

def _main_():
    if len (sys.argv) != 4:
        print("Script designed to adjust RepeatMasker GFF using the TElib header.")
        print("Basic usage: AdjustingGFF_RM.py RMout.gff TElib.fasta RMout_adjusted.gff")
        quit()

    GFFrm = sys.argv[1]
    TElib = sys.argv[2]
    GFFout = sys.argv[3]
    _Run_(GFFrm, TElib, GFFout)

_main_()

#END
