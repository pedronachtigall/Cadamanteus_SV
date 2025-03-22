

# Crotalus adamanteus structural variant
[![Published in MBE](https://img.shields.io/badge/published%20in-MBE-blue)](https://doi.org/10.1093/molbev/msaf058)
[![Data available in the Fisghare](https://img.shields.io/badge/data%20available%20in%20the-figshare-red)](https://figshare.com/projects/Eastern_diamondback_rattlesnake_Crotalus_adamanteus_-_haplotype-resolved_genome_assembly/200614)

This repository contains commands and scripts used in the manuscript "A Segregating Structural Variant Defines Novel Venom
Phenotypes in the Eastern Diamondback Rattlesnake" published in *Molecular Biology and Evolution*.

All datasets used in the present study are detailed in the Supplementary file of the published manuscript.

## Genome assembly
The genome assembly pipeline is described in the "[GenomeAssembly.md](https://github.com/pedronachtigall/Cadamanteus_SV/blob/main/GenomeAssembly.md)" file.

## Genome annotation
### Repeat annotation
The repeat annotation was performed using [RepeatModeler2](https://github.com/Dfam-consortium/RepeatModeler), to generate a *de novo* species-specific library, and [RepeatMasker](https://github.com/Dfam-consortium/RepeatMasker), to perform the annotation. We complement the species-sepcific TE library using a curated TE library of snakes previosuly published <sup>[Castoe et al., 2013](https://doi.org/10.1073/pnas.1314475110)</sup>.

The pipeline with commands and scripts used to perform the repeat annotation is decribed in the following tutorial: https://github.com/pedronachtigall/Repeat-annotation-pipeline

For this step, we used the primary genome assembly to perform the repeat annotation and the soft-masked primary genome assembly as the source for gene annotation.

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

```

The toxin genes were manually reviewed and curated following the "[Checking annotations](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Guide#checking-annotations)" section of the ToxCodAn-Genome's guide to ensure a comprehensive toxin annotation.

The toxin annotation was lifted from primary to both haplotypes using LiftOff.
```
liftoff -g Cadam_primary_chromosomes.toxin.gtf Cadam_hap1_chromosomes.fasta Cadam_primary_chromosomes.fasta -o Cadam_hap1_chromosomes.toxin.liftoff.gtf
liftoff -g Cadam_primary_chromosomes.toxin.gtf Cadam_hap2_chromosomes.fasta Cadam_primary_chromosomes.fasta -o Cadam_hap2_chromosomes.toxin.liftoff.gtf
```

## Estimate expression level
We used Bowtie2 and RSEM to estimate the expression level of venom-gland from 18 individuals <sup>[Hogan et al., 2024](https://doi.org/10.1073/pnas.2313440121)</sup> as described in the "[Estimating expression level](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Guide#estimating-expression-level)" section of the ToxCodAn-Genome's guide. We mapped the reads into the CDSs retrieved from annotation files and set the ```mismatch-rate``` parameter to 0.02.

We used the merged annotation and the primary assembly to retrieve the toxin and nontoxin set of CDSs. The CDSs were retrieved using [GffRead](https://github.com/gpertea/gffread). The datasets used to estimate the expression level of venom-gland is in the file "vg_rna.list". We also trimmed adapters and removed low-quality reads using trim_galore.

```
#retrieve CDSs
gffread -x Cadam_primary_annotation_cds.fasta -g Cadam_primary_chromosomes.fasta Cadam_primary_chromosomes.MERGE.gtf

#estimate expression level
rsem-prepare-reference --bowtie2 Cadam_primary_annotation_cds.fasta CADAM
for i in vg_rna.list; do
	rsem-calculate-expression -p 20 --paired-end --bowtie2 --bowtie2-mismatch-rate 0.02 ${i}_tg/${i}_R1_val_1.fq.gz ${i}_tg/${i}_R2_val_2.fq.gz CADAM ${i}_rsem
done
```

## Gene tree
When needed (e.g., to analyze SVMPs and SVSPs), we inferred the gene tree using their CDSs. We used [MAFFT](https://github.com/GSLBiotech/mafft) and [IQ-TREE](https://github.com/Cibiv/IQ-TREE) for alignment of sequences and infer gene tree, respectively. The input files, which consist of a fasta file with sequences from the target gene family, were set manually.

```
mafft --auto GENE.fasta > GENE.ALIGNED.fasta
iqtree -s GENE.ALIGNED.fasta -m TEST -bb 1000 -alrt 1000
```

The trees were inspected and adjusted using [FigTree](https://github.com/rambaut/figtree/).

## Comparing haplotypes

### Sequence-level
We used [syri](https://github.com/schneebergerlab/syri) and followed its guide to compare haplotypes at sequence-level.
```
minimap2 -ax asm5 --eqx Cadam_hap1_chromosomes.fasta Cadam_hap2_chromosomes.fasta > out.sam
syri -c out.sam -r Cadam_hap1_chromosomes.fasta -q Cadam_hap2_chromosomes.fasta -k -F S
plotsr --sr syri.out --genomes genomes.txt -H 8 -W 5 -o ALL_hap1Xhap2.pdf
```

### Gene-level
We used [GENESPACE](https://github.com/jtlovell/GENESPACE) and followed its guide to compare haplotypes at gene-level.

### Microsynteny of toxin genes
We used BLAST to perform a microsynteny of the toxin gene regions between haplotypes and primary assemblies.

Regions screened:
- **SVMP**
  - Primary - Cadam_mi-1:26000000-26900000
  - Haplotype_1 - Cadam_hap1_mi-1:20900000-21700000
  - Haplotype_2 - Cadam_hap2_mi-1:22200000-23100000

- **SVSP**
  - Primary - Cadam_mi-2:8900000-10200000
  - Haplotype_1 - Cadam_hap1_mi-2:8800000-9600000
  - Haplotype_2 - Cadam_hap2_mi-2:8800000-9900000

- **MYO**
  - Primary - Cadam_ma-2:286300000-288700000
  - Haplotype_1 - Cadam_hap1_ma-2:282500000-283000000
  - Haplotype_2 - Cadam_hap2_ma-2:284250000-286100000

- **PLA2**
  - Primary - Cadam_mi-8:295000-330000
  - Haplotype_1 - Cadam_hap1_mi-8:300000-335000
  - Haplotype_2 - Cadam_hap2_mi-8:294000-330000

- **CTL**
  - Primary - Cadam_mi-2:172000-1490000
  - Haplotype_1 - Cadam_hap1_mi-2:74000-1561000
  - Haplotype_2 - Cadam_hap2_mi-2:24000-1530000

```
#retrieve regions from assembly
samtools faidx Cadam_primary_chromosomes.fasta "${chr}:${start}-${end}" | awk '/^>/{print ">pri_${GENE}"; next}{print}' > ${GENE}_pri.fasta
samtools faidx Cadam_hap1_chromosomes.fasta "${chr}:${start}-${end}" | awk '/^>/{print ">hap1_${GENE}"; next}{print}' > ${GENE}_hap1.fasta
samtools faidx Cadam_hap2_chromosomes.fasta "${chr}:${start}-${end}" | awk '/^>/{print ">hap2_${GENE}"; next}{print}' > ${GENE}_hap2.fasta

#run BLAST
makeblastdb -dbtype nucl -in ${GENE}_pri.fasta -out blastDB/PRI
makeblastdb -dbtype nucl -in ${GENE}_hap1.fasta -out blastDB/HAP1
makeblastdb -dbtype nucl -in ${GENE}_hap2.fasta -out blastDB/HAP2
blastn -num_threads 6 -perc_identity 90 -query ${GENE}_hap1.fasta -db blastDB/PRI -out ${GENE}_priXhap1_blast.out -outfmt 6
blastn -num_threads 6 -perc_identity 90 -query ${GENE}_hap2.fasta -db blastDB/PRI -out ${GENE}_priXhap2_blast.out -outfmt 6
blastn -num_threads 6 -perc_identity 90 -query ${GENE}_hap1.fasta -db blastDB/HAP2 -out ${GENE}_hap1Xhap2_blast.out -outfmt 6
```

We used the BLAST output to plot alignments using [ggplot2](https://ggplot2.tidyverse.org/) in R.

## Exon-capture data analysis
We used a set of exon-capture data available for 139 individuals of *C. adamanteus* sampled throughout the species distribution <sup>[Margres et al., 2017](https://doi.org/10.1534/genetics.117.202655);[Margres et al., 2019](https://doi.org/10.1093/molbev/msy207)</sup>.

We trimmed adapters and removed low-quality reads.
```
for i in capture_data.list; do
	trim_galore --paired --phred33 --length 75 -q 25 --stringency 1 -e 0.1 -o ${i}_tg Sample_${i}/${i}*R1*.fastq.gz Sample_${i}/${i}*R2*.fastq.gz
done
```

We mapped reads using [Bowtie2](https://github.com/BenLangmead/bowtie2) with the primary genome assembly as reference.
```
mkdir B2_index
bowtie2-build Cadam_primary_chromosomes.fasta B2_index/Cadam

for i in capture_data.list; do
	bowtie2 -p 20 --very-sensitive --no-unal --no-mixed --no-discordant -k 10 -X 700 -x B2_index/Cadam -1 ${i}_tg/${i}*_1.fq.gz -2 ${i}_tg/${i}*_2.fq.gz | samtools view -u - | samtools sort -o ${i}.mapped.bam -
	samtools view -@ 40 -q 30 -b ${i}.bam -o ${i}.mapped.q30.bam
	samtools index ${i}.mapped.q30.bam
done
```

### Genotyping SVMP
We used the mapped reads (filtered for multi-mapped reads) to retrieve coverage for SVMP genes.
```
for i in capture_data.list; do
	samtools depth -b SVMP_region.bed ${i}.mapped.q30.bam > ${i}.SVMP.depth.txt
done
```
We calculated the average of coverage for all SVMP genes and comapred to each SVMP gene for each individual using the [```SVMP_counter.py```](https://raw.githubusercontent.com/pedronachtigall/Cadamanteus_SV/refs/heads/main/scripts/SVMP_counter.py) script. Then, we manually checked the output and the read coverage using IGV to infer the final genotype.

### Estimating CNV of MYO
We used the mapped reads (not filtered for multi-mapped reads) to retrieve coverage for MYO and 10 nontoxin genes (i.e., ATPSynLipid-1, ATPase-lys70, CD63, Calreticulin, DAZ-2, GADD45, Glutaredoxin-1, Leptin-1, PDI, and Nexin-2) located at the same chromosome (i.e., ma-2).
```
for i in capture_data.list; do
	samtools depth -b MYO_region.bed ${i}.bam > ${i}.MYO.depth.txt
done
```
We calculated the average of coverage of MYO to the average of coverage of the nontoxin genes to estimate the CNV for each individual using the [```MYO_counter.py```](https://raw.githubusercontent.com/pedronachtigall/Cadamanteus_SV/refs/heads/main/scripts/MYO_counter.py) script.

## Cite
If you follow the pipelines and/or use any of the scripts in this repository, please cite:

Nachtigall et al. (2025) A Segregating Structural Variant Defines Novel Venom Phenotypes in the Eastern Diamondback Rattlesnake. Molecular Biology and Evolution. [https://doi.org/10.1093/molbev/msaf058](https://doi.org/10.1093/molbev/msaf058)

```
@article{Nachtigall:2025,
    author = {Nachtigall, Pedro G and Nystrom, Gunnar S and Broussard, Emilie M and Wray, Kenneth P and Junqueira-de-Azevedo, Inácio L M and Parkinson, Christopher L and Margres, Mark J and Rokyta, Darin R},
    title = {A Segregating Structural Variant Defines Novel Venom Phenotypes in the Eastern Diamondback Rattlesnake},
    journal = {Molecular Biology and Evolution},
    pages = {msaf058},
    year = {2025},
    month = {03},
    issn = {1537-1719},
    doi = {10.1093/molbev/msaf058},
    url = {https://doi.org/10.1093/molbev/msaf058},
    eprint = {https://academic.oup.com/mbe/advance-article-pdf/doi/10.1093/molbev/msaf058/62445184/msaf058.pdf},
}
```
