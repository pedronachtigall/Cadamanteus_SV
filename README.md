

# Crotalus adamanteus structural variant
[![Published in MBE](https://img.shields.io/badge/published%20in-MBE-blue)](https://doi.org/10.1093/gigascience/giad116)
[![Data available in the Fisghare](https://img.shields.io/badge/data%20available%20in%20the-figshare-red)](https://figshare.com/projects/Eastern_diamondback_rattlesnake_Crotalus_adamanteus_-_haplotype-resolved_genome_assembly/200614)

This repository contains commands and scripts used in the manuscript "A Segregating Structural Variant Defines Novel Venom
Phenotypes in the Eastern Diamondback Rattlesnake" published in *Molecular Biology and Evolution*.

All datasets used in the present study are detailed in the Supplementary file of the published manuscript.

## Genome Assembly
The genome assembly is described in the "[GenomeAssembly.md](https://github.com/pedronachtigall/Cadamanteus_SV/blob/main/GenomeAssembly.md)" file.

## Genome annotation
### Repeat annotation
The repeat annotation was performed using [RepeatModeler2](https://github.com/Dfam-consortium/RepeatModeler), to generate a *de novo* species-specific library, and [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker), to perform the annotation. We complement the species-sepcific TE library using a curated TE library of snakes previosuly published <sup>[Castoe et al., 2013](https://doi.org/10.1073/pnas.1314475110)</sup>.

The pipeline with commands and scripts used to perform the repeat annotation is decribed in the following tutorial: https://github.com/pedronachtigall/Repeat-annotation-pipeline

For this step, we used the primary assembly to perform the repeat annotation and the soft-masked primary assembly as the source for gene annotation.

### Gene annotation
The gene annotation was performed using the [funannotate](https://github.com/nextgenusfs/funannotate) pipeline. We followed the commands decribed in the "[NonToxin annotation](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Guide#nontoxin-annotation)" section of the ToxCodAn-Genome's guide.

For this step, we used the soft-maked primary assembly to perform the gene annotation using funannotate and set the RNA-seq of several tissues to be used as transcript evidence.

Then, we lifted annotations from primary to both haplotypes using [LiftOff](https://github.com/agshumate/Liftoff).

```
liftoff -g Cadam_primary_chromosomes.funannotate.gff Cadam_hap1_chromosomes.fasta Cadam_primary_chromosomes.fasta -o Cadam_hap1_chromosomes.funannotate.liftoff.gff
liftoff -g Cadam_primary_chromosomes.funannotate.gff Cadam_hap2_chromosomes.fasta Cadam_primary_chromosomes.fasta -o Cadam_hap2_chromosomes.funannotate.liftoff.gff
```

The mitochondrial genome was annotated using [MitoZ](https://github.com/linzhi2013/MitoZ).
```
MitoZ.py annotate --genetic_code auto --clade Chordata --outprefix mitogenome_annotation --thread_number 20 --fastafile Cadam_mitogenome.fasta
```

### Toxin annotation
We used [ToxCodAn-Genome](https://github.com/pedronachtigall/ToxCodAn-Genome) to annotate toxins.

For this step, we used the primary assembly, the venom-gland transcriptome data annotated using [ToxCodAn](https://github.com/pedronachtigall/ToxCodAn) and the [Viperidae](https://raw.githubusercontent.com/pedronachtigall/ToxCodAn-Genome/main/Databases/Viperidae_db.fasta) database. We used the venom-gland transcriptome data of a male and female individuals (SRR21096543 and SRR21096547).

First, we assembled the venom-gland transcriptome. For the transcriptome assembly, we used the "[TRassembly.py](https://raw.githubusercontent.com/pedronachtigall/ToxCodAn-Genome/refs/heads/main/bin/TRassembly.py)" script from ToxCodAn-Genome, which uses the genome-guided mode of Trinity and StringTie and the *de novo* mode of Trinity and rnaSPAdes.

```
#trim reads
trim_galore --paired --phred33 --length 75 -q 25 --stringency 1 -e 0.1 -o SRR21096543_tg SRR21096543_R1.fastq.gz SRR21096543_R2.fastq.gz
trim_galore --paired --phred33 --length 75 -q 25 --stringency 1 -e 0.1 -o SRR21096547_tg SRR21096547_R1.fastq.gz SRR21096547_R2.fastq.gz

#assembly transcripts
TRassembly.py -g Cadam_primary_chromosomes.fasta -r SRR21096543_R1_val_1.fastq.gz,SRR21096543_R2_val_2.fastq.gz -c 20 -M 20G --output SRR21096543_TRassembly
TRassembly.py -g Cadam_primary_chromosomes.fasta -r SRR21096547_R1_val_1.fastq.gz,SRR21096547_R2_val_2.fastq.gz -c 20 -M 20G --output SRR21096547_TRassembly
```

Annotate the transcriptome using ToxCodAn.
```
wget https://github.com/pedronachtigall/ToxCodAn/raw/refs/heads/master/models.zip
unzip models.zip
toxcodan.py -c 10 -t SRR21096543_TRassembly/transcripts.fasta -m models -o SRR21096543_ToxCodAn -s SRR21096543
toxcodan.py -c 10 -t SRR21096547_TRassembly/transcripts.fasta -m models -o SRR21096547_ToxCodAn -s SRR21096547

cat SRR21096543_ToxCodAn/SRR21096543_Toxins_cds_RedundancyFiltered.fasta SRR21096543_ToxCodAn/SRR21096543_PutativeToxins_cds_SPfiltered.fasta \
SRR21096547_ToxCodAn/SRR21096547_Toxins_cds_RedundancyFiltered.fasta SRR21096547_ToxCodAn/SRR21096547_PutativeToxins_cds_SPfiltered.fasta > Cadam_VG_toxins.toxcodan.fasta
```

Then, run ToxCodAn-Genome setting the annotated toxins from the vneomg-land transcriptome and the Viperidae database.
```
wget https://raw.githubusercontent.com/pedronachtigall/ToxCodAn-Genome/main/Databases/Viperidae_db.fasta

toxcodan-genome.py -c 10 -g Cadam_primary_chromosomes.fasta -C Cadam_VG_toxins.toxcodan.fasta -d Viperidae_db.fasta

#the annotation file to be reviewed is: ToxCodAnGenome_output/toxin_annotation.gtf
```

The toxin genes were manually reviewed and curated following the "[Checking annotations](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Guide#checking-annotations)" section of the ToxCodAn-Genome's guide to ensure a comprehensive toxin annotation.

The toxin annotation was lifted from primary to both haplotypes using LiftOff.
```
liftoff -g Cadam_primary_chromosomes.toxin.gtf Cadam_hap1_chromosomes.fasta Cadam_primary_chromosomes.fasta -o Cadam_hap1_chromosomes.toxin.liftoff.gtf
liftoff -g Cadam_primary_chromosomes.toxin.gtf Cadam_hap2_chromosomes.fasta Cadam_primary_chromosomes.fasta -o Cadam_hap2_chromosomes.toxin.liftoff.gtf
```

## Estimate expression level
We used Bowtie2 and RSEM to estimate the expression level of venom-gland from several individuals as described in the "[Estimating expression level](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Guide#estimating-expression-level)" section of the ToxCodAn-Genome's guide. We mapped the reads into the CDSs retrieved from annotation files and set the ```mismatch-rate``` parameter to 0.02.

We used the merged annotation and the primary assembly to retrieve the toxin and nontoxin set of CDSs. The CDSs were retrieved using [GffRead](https://github.com/gpertea/gffread). The datasets used to estimate the expression level of venom-gland is in the file "vg_rna.list". We also trimmed adapters and removed low-quality reads using trim_galore.

```
#retrieve CDSs
gffread -x Cadam_primary_annotation_cds.fasta -g Cadam_primary_chromosomes.fasta Cadam_primary_chromosomes.MERGE.gtf

#estimate expression level
rsem-prepare-reference --bowtie2 Cadam_primary_annotation_cds.fasta CADAM
for i in ; do
	rsem-calculate-expression -p 20 --paired-end --bowtie2 --bowtie2-mismatch-rate 0.02 rna/${i}_tg/${i}_R1_val_1.fq.gz rna/${i}_tg/${i}_R2_val_2.fq.gz CADAM ${i}_rsem
done
```

## Gene tree
When needed (e.g., to analyze SVMPs and SVSPs), we inferred the gene tree using their CDSs. We used [MAFFT](https://github.com/GSLBiotech/mafft) and [IQ-TREE](https://github.com/Cibiv/IQ-TREE) for alignment of sequences and infer gene tree, respectively. The input files, which consist of a fasta file with sequences from the target gene family, were set manually.

```
mafft --auto GENE.fasta > GENE.ALIGNED.fasta
iqtree -s GENE.ALIGNED.fasta -m TEST -bb 1000 -alrt 1000
```

The trees were inspected and adjusted using [FigTree](https://github.com/rambaut/figtree/).

## Cite
If you follow the pipelines and/or scripts in this repository, please cite:
```
ADD bibtext here
```
