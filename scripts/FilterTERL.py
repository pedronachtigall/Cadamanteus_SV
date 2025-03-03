#!/usr/bin/env python3
'''
Script designed to filter out "NonTE" sequences as classified by TERL to generate the final TElib.
 #PS: it is based on the TERL adjusted to generate a "better" header in the output

Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)

inputs:
TERL_classified.fasta
Input_used_for_TERL.fasta

output:
TERL_filtered.fasta

'''

import sys
from Bio import SeqIO

def _ParseTERL_(TERL):
    final = []
    a = open(TERL,"r")
    for line in a:
        if line.startswith(">"):
            if "NonTE" in line:
                id = line.split("__")[0].replace(">","")
                final.append(id)
    a.close()
    return set(final)

def _ParseTElib_(TEin):
    final = {}
    for record in SeqIO.parse(TEin,"fasta"):
        id = str(record.id)
        seq = str(record.seq)
        final[id] = seq
    return final

def _Run_(TERL, TEin, TEout):
    terl = _ParseTERL_(TERL)
    TElib = SeqIO.to_dict(SeqIO.parse(TEin, 'fasta'))
    #TElib = _ParseTElib_(TEin)
    TEOUT = open(TEout, "w")
    written = set([])
    for k in TElib:
    #for k in TElib.keys():
        if k not in terl and k not in written:
            TEOUT.write(">"+k+"\n"+str(TElib[k].seq)+"\n")
            #TEOUT.write(">"+k+"\n"+TElib[k]+"\n")
            written.add(k)
    TEOUT.close()
    print(len(terl),"sequences removed based on NonTE classification")
    print("   >>> final set has",len(written),"sequences")

def _Main_():
    if len (sys.argv) != 4:
        print("Script designed to filter out \"NonTE\" sequences as classified by TERL to generate the final TElib.")
        print("Basic usage: FilterTERL.py TERLoutput input_TElib.fasta filtered.fasta")
        quit()

    TERL = sys.argv[1]
    TEin = sys.argv[2]
    TEout = sys.argv[3]
    _Run_(TERL, TEin, TEout)

_Main_()

#END
