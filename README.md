

# Crotalus adamanteus structural variant
[![Published in MBE](https://img.shields.io/badge/published%20in-MBE-blue)](https://doi.org/10.1093/gigascience/giad116)
[![Data available in the Fisghare](https://img.shields.io/badge/data%20available%20in%20the-figshare-red)](https://figshare.com/projects/Eastern_diamondback_rattlesnake_Crotalus_adamanteus_-_haplotype-resolved_genome_assembly/200614)

This repository contains commands and scripts used in the manuscript "A Segregating Structural Variant Defines Novel Venom
Phenotypes in the Eastern Diamondback Rattlesnake" published in *Molecular Biology and Evolution*.

## Genome Assembly

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

For this step, we used the primary assembly, the venom-gland transcriptome data annotated using [ToxCodAn](https://github.com/pedronachtigall/ToxCodAn) and the [Viperidae](https://raw.githubusercontent.com/pedronachtigall/ToxCodAn-Genome/main/Databases/Viperidae_db.fasta) database.

First, assembly the venomg-land transcriptome and annotate it using ToxCodAn.
```

```

Then, run ToxCodAn-Genome setting the annotated toxins from the vneomg-land transcriptome and the Viperidae database.
```
wget https://raw.githubusercontent.com/pedronachtigall/ToxCodAn-Genome/main/Databases/Viperidae_db.fasta



```

The toxin genes were manually reviewed and curated following the "[Checking annotations](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Guide#checking-annotations)" section of the ToxCodAn-Genome's guide to ensure a comprehensive toxin annotation.

The toxin annotation was lifted from primary to both haplotypes using LiftOff.
```
liftoff -g Cadam_primary_chromosomes.toxin.gtf Cadam_hap1_chromosomes.fasta Cadam_primary_chromosomes.fasta -o Cadam_hap1_chromosomes.toxin.liftoff.gtf
liftoff -g Cadam_primary_chromosomes.toxin.gtf Cadam_hap2_chromosomes.fasta Cadam_primary_chromosomes.fasta -o Cadam_hap2_chromosomes.toxin.liftoff.gtf
```

## Estimate expression level
We used Bowtie2 and RSEM to estimate the expression level of venom-gland from several individuals as described in the "[Estimating expression level](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Guide#estimating-expression-level)" section of the ToxCodAn-Genome's guide. We mapped the reads into the CDSs retrieved from annotation files and set the ```mismatch-rate``` parameter to 0.02.

```
rsem-prepare-reference --bowtie2 Cadam_primary_annotation_cds.fasta CADAM
for i in DRR0105 DRR0044 DRR0106 DRR0107 DRR0108 KW0944 KW1264 KW1942 KW2161 KW2170 KW2171 KW2184 MM0114 MM0127 MM0143 MM0198 TJC1661 TJC1665; do
	rsem-calculate-expression -p 20 --paired-end --bowtie2 --bowtie2-mismatch-rate 0.02 ../rna/${i}_tg/${i}_R1_val_1.fq.gz ../rna/${i}_tg/${i}_R2_val_2.fq.gz CADAM ${i}_rsem
done
```

## Cite

```
ADD bibtext here
```
