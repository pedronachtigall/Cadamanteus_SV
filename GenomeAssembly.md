## Genome Assembly

Pipeline to perform the haplotype-resolved genome assembly of *Crotalus adamanteus*.

The raw data is listed below:
| Sample ID | Data type | NCBI accession |
| :-------- | :-------: | :------------: | 
| DRR0105   | HiFi | SRR21092035, SRR21092036 |
| DRR0105   | HiC | SRR28357486 |

### Trim adapters
We used [cutadapt](https://github.com/marcelm/cutadapt) to remove reads with adapters from HiFi data and [Trim_Galore](https://github.com/FelixKrueger/TrimGalore) to trim adapters and low-quality reads from HiC data.
```
#trim and merge hifi reads
cutadapt -j 20 -n 3 -O 35 --revcomp --discard-trimmed --anywhere="ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT" --anywhere="ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT" SRR21092035.hifi_reads.fastq.gz --output SRR21092035_HiFi_cutadapt.fastq
cutadapt -j 20 -n 3 -O 35 --revcomp --discard-trimmed --anywhere="ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT" --anywhere="ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT" SRR21092036.hifi_reads.fastq.gz --output SRR21092036_HiFi_cutadapt.fastq
cat SRR21092035_HiFi_cutadapt.fastq SRR21092036_HiFi_cutadapt.fastq > Cadam.hifi.trimmed.fastq

#trim hic reads
trim_galore --paired --phred33 --length 75 -q 25 --stringency 1 -e 0.1 -o SRR28357486_hic_tg SRR28357486_R1.fastq.gz SRR28357486_R2.fastq.gz
```

### Remove of contaminants
We checked and removed any bacterial and/or human reads (if needed) in the hifi data using [kraken2](https://github.com/DerrickWood/kraken2). TaxonomyID of Sauria 32561 was used to retrieve reads not matching to bacteria or human. We used a custom kraken2 database comprising the standard libraries and other squamata genomes.
```
mkdir kraken2 && cd kraken2
kraken2 --threads 20 --db krakendb ../Cadam.hifi.trimmed.fastq --report krakendb_report.txt --output krakendb_output.txt
extract_kraken_reads.py -k krakendb_output.txt --report krakendb_report.txt -s ../Cadam.hifi.trimmed.fastq -t 32561 -o ../Cadam.hifi.fastq --include-children
```

### Mitochondrial genome assembly
We used [MITGARD](https://github.com/pedronachtigall/MITGARD) to perform the mitogenome assembly. The available mitochondrial genome of *C. adamanteus* (NC_041524.1) was used as the reference.
```
mkdir MITGARD_output && cd MITGARD_output
MITGARD-LR.py -s Cadam_mitogenome -m pacbio_hifi -r ../Cadam.hifi.fastq -R NC_041524.1.fasta
```

### Draft genome assembly
We used [hifiasm](https://github.com/chhylp123/hifiasm) to perform the draft genome assembly using HiFi and HiC reads to acquire the primary and both resolved haplotypes.
```
hifiasm -o Cadam_DRR0105 -t32 --h1 SRR28357486_hic_tg/SRR28357486_R1_val_1.fq --h2 SRR28357486_hic_tg/SRR28357486_R2_val_2.fq Cadam.hifi.fastq
awk '/^S/{print ">"$2;print $3}' Cadam_DRR0105.hic.p_ctg.gfa > Cadam_DRR0105.hic.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' Cadam_DRR0105.hic.hap1.p_ctg.gfa > Cadam_DRR0105.hic.hap1.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' Cadam_DRR0105.hic.hap2.p_ctg.gfa > Cadam_DRR0105.hic.hap2.p_ctg.fasta
```

### Scaffold draft genome
We used [YaHS](https://github.com/c-zhou/yahs) to scaffold the primary genome assembled by hifiasm. To map reads against the draft genome, we used [chromap](https://github.com/haowenz/chromap).
```
mkdir YAHS_primary && cd YAHS_primary
ln -s ../Cadam_DRR0105.hic.p_ctg.fasta .
ln -s ../SRR28357486_hic_tg/SRR28357486_R1_val_1.fq hic_R1.fastq
ln -s ../SRR28357486_hic_tg/SRR28357486_R2_val_2.fq hic_R2.fastq

chromap -i -r Cadam_DRR0105.hic.p_ctg.fasta -o Cadam_DRR0105.hic.p_ctg.fasta.index
samtools faidx Cadam_DRR0105.hic.p_ctg.fasta

chromap --preset hic -r Cadam_DRR0105.hic.p_ctg.fasta -x Cadam_DRR0105.hic.p_ctg.fasta.index --remove-pcr-duplicates -1 hic_R1.fastq -2 hic_R2.fastq --SAM -o aligned.sam -t 32
samtools view -@ 32 -bh aligned.sam | samtools sort -@ 32 -n > aligned.bam
rm aligned.sam

yahs Cadam_DRR0105.hic.p_ctg.fasta aligned.bam
```

### Review of scaffolded genome
We used [Juicer](https://github.com/aidenlab/juicer) and the [3D-DNA](https://github.com/aidenlab/3d-dna) to manually review the scaffolded genome following the standard [DNA Genome Assembly Cookbook](https://aidenlab.org/assembly/manual 180322.pdf) instructions.
```
#prepare folder to run juicer
mkdir Juicer && cd Juicer
mkdir fastq references restriction_sites splits

#softlink and index the genome sequence in the references/ folder
cd references
ln -s ../../yahs.out_scaffolds_final.fa .
bwa index yahs.out_scaffolds_final.fa

#create a sizes.genome file containing all the chromosome/scaffold sizes
samtools faidx yahs.out_scaffolds_final.fa
cut -f1,2 yahs.out_scaffolds_final.fa.fai > ../sizes.genome

#create restriction enzyme proposed cutting sites
cd ../restriction_sites
python /path/to/juicer-1.6/misc/generate_site_positions.py Arima yahs.out_scaffolds_final.fa ../references/yahs.out_scaffolds_final.fa

#get your HiC reads ready
cd ../fastq
ln -s ../../../SRR28357486_hic_tg/SRR28357486_R1_val_1.fq hic_R1.fastq
ln -s ../../../SRR28357486_hic_tg/SRR28357486_R2_val_2.fq hic_R2.fastq

#splits reads
cd ../splits/
split -a 3 -l 90000000 -d --additional-suffix=_R2.fastq ../fastq/hic_R2.fastq &
split -a 3 -l 90000000 -d --additional-suffix=_R1.fastq ../fastq/hic_R1.fastq &
cd ..

mkdir -p scripts/common
cp /path/to/juicer-1.6/CPU/common/* scripts/common
cp /path/to/juicer_tools.2.20.00.jar scripts/common/juicer_tools.jar

#run juicer
cd ../..
/path/to/juicer-1.6/CPU/juicer.sh -d Juicer -p Juicer/sizes.genome -y Juicer/restriction_sites/Cadam_chr.primary.fasta_Arima.txt -z Juicer/references/yahs.out_scaffolds_final.fa -g yahs.out_scaffolds_final.fa -t 32 -D Juicer

#visualize candidate assembly
cd Juicer/aligned/
awk -f /path/to/3d-dna/utils/generate-assembly-file-from-fasta.awk ../references/yahs.out_scaffolds_final.fa > yahs.out_scaffolds_final.fa.assembly
/path/to/3d-dna/visualize/run-assembly-visualizer.sh yahs.out_scaffolds_final.fa.assembly merged_nodups.txt
```
The .hic and .assembly files must be reviewed using [JuiceBox](https://github.com/aidenlab/Juicebox).

To help in the review process, we checked for telomeric sequences within scaffolds using [tidk](https://github.com/tolkit/telomeric-identifier) and the conserved telomeric motif of vertebrates (TTAGGG). We also mapped the assembled mitochondrial genome to check for mitochondrial contigs/scaffolds that can be removed.
```
tidk search -s TTAGGG --dir tidk_out --output Cadam_pri_prereview yahs.out_scaffolds_final.fa
tidk search -s TTAGGG --dir tidk_out --output Cadam_pri_prereview --extension bedgraph yahs.out_scaffolds_final.fa
tidk plot --csv tidk_out/Cadam_pri_prereview_telomeric_repeat_windows.csv --output tidk_out/Cadam_pri_prereview_telomeric_repeat_windows

mkdir mitocheck && cd mitocheck
makeblastdb -in yahs.out_scaffolds_final.fa -out blastDB/genome -dbtype nucl
blastn -query ../../mitogenome/Cadam_mitogenome.fasta -db blastDB/genome -out blast.out -evalue 1E-6 -qcov_hsp_perc 40 -num_threads 20 -outfmt 6
```

After reviewing, we generate the final scaffolded primary assembly.
```
#the working directory must be Juicer/aligned/
/path/to/run-asm-pipeline-post-review.sh –r draft.final.review.assembly ../references/yahs.out_scaffolds_final.fa merged_nodups.txt
```

### Assign chromosomes
To identify chromosomes, we used a set of chromosome-specific markers (NCBI accessions SAMN00177542 and SAMN00152474) of snakes<sup>[Matsubara
et al., 2006](https://doi.org/10.1073/pnas.0605274103)</sup> and the chromosome-level genome available of the closely related species *C. viridis*<sup>[Schield et al., 2019](http://www.genome.org/cgi/doi/10.1101/gr.240952.118)</sup>. We used [BLAST](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) and [minimap2](https://github.com/lh3/minimap2) to perform this analysis and manually checked the outputs.
```
makeblastdb -in yahs.out_scaffolds_final.FINAL.fasta -out blastDB/ALL -dbtype nucl
blastn -query /path/to/markers_sequence.fasta -out blast_markers.out -db blastDB/ALL -num_threads 20 -max_target_seqs 1 -outfmt 6
minimap2 -ax splice yahs.out_scaffolds_final.FINAL.fasta /path/to/markers_sequence.fasta > mm2_markers.sam
minimap2 -k19 -w19 -t32 /path/to/Crovir_chr_only.fasta yahs.out_scaffolds_final.FINAL.fasta > mm2_cvir.paf
```

After characterizing chromosomes, we renamed the scaffolds and split the file into chromsomes and unplaced for downstream analysis.

### Confirm sex chromosomes
To double-check the sex chromosomes, we mapped whole genome seqeuncing data available for a male and female individuals. We used [bwa](https://github.com/lh3/bwa) and checked the read coverage manually using [IGV](https://igv.org/).

The raw data is listed below:
| Sample ID | Sex | NCBI accession |
| :-------- | :-------: | :------------: | 
| KW0944 | Male | SRR12802469 |
| KW1264 | Female | SRR12802470 |

```
mkdir wgs
for i in KW0944 KW1264; do
	bwa mem -t 40 Cadam_primary_chromosomes.fasta wgs/${i}/${i}_tg/${i}_R1_val_1.fq.gz wgs/${i}/${i}_tg/${i}_R2_val_2.fq.gz > wgs/${i}_align.sam
	samtools view -@ 40 -b -S -F 4 -o wgs/${i}_mapped.bam wgs/${i}_align.sam
	rm wgs/${i}_align.sam
	samtools sort -@ 40 wgs/${i}_mapped.bam -o wgs/${i}_mapped_sorted.bam
	rm wgs/${i}_mapped.bam
	samtools view -@ 40 -q 30 -b wgs/${i}_mapped_sorted.bam -o wgs/${i}_mapped_sorted.q30.bam
	samtools index wgs/${i}_mapped_sorted.q30.bam
	bedtools genomecov -ibam wgs/${i}_mapped_sorted.q30.bam -bg > wgs/${i}_q30_cov.bedGraph
done
```

### Haplotype-resolved
The haplotype-resolved assemblies were generated using [RagTag](https://github.com/malonge/RagTag) with each haplotype contigs as a query and the primary chromosome-level assembly as a reference.
```
mkdir haplotype_resolved && cd haplotype_resolved
ragtag.py scaffold -t 10 Cadam_primary_chromosomes.fasta ../Cadam.hic.hap1.p_ctg.fasta -o ragtag_hifiasm_hap1
ragtag.py scaffold -t 10 Cadam_primary_chromosomes.fasta ../Cadam.hic.hap2.p_ctg.fasta -o ragtag_hifiasm_hap2
```

### Checking genome quality
#### BUSCO
We used [BUSCO](https://busco.ezlab.org/) with the Tetrapoda database to assess the completeness.
```
busco -i Cadam_primary_chromosomes.fasta -m genome -l tetrapoda_odb10 -c 20 -o busco_primary
busco -i Cadam_hap1_chromosomes.fasta -m genome -l tetrapoda_odb10 -c 20 -o busco_hap1
busco -i Cadam_hap2_chromosomes.fasta -m genome -l tetrapoda_odb10 -c 20 -o busco_hap2
```

#### Inspector
We used [Inspector](https://github.com/Maggi-Chen/Inspector) to assess genomic statistics and assembly quality by calculating the QV score.
```
inspector.py -c Cadam_primary_chromosomes.fasta -r Cadam.hifi.fastq -o inspector_pri_out/ --datatype hifi --thread 10
inspector.py -c Cadam_hap1_chromosomes.fasta -r Cadam.hifi.fastq -o inspector_hap1_out/ --datatype hifi --thread 10
inspector.py -c Cadam_hap2_chromosomes.fasta -r Cadam.hifi.fastq -o inspector_hap2_out/ --datatype hifi --thread 10
```

#### VerityMap
We used [VerityMap](https://github.com/ablab/VerityMap) to detect error-prone regions in the assemblies.
```
python /path/to/VerityMap/veritymap/main.py -t 20 -d hifi-diploid --reads Cadam.hifi.fastq -o Cadam_pri_VM Cadam_primary_chromosomes.fasta
python /path/to/VerityMap/veritymap/main.py -t 20 -d hifi-diploid --reads Cadam.hifi.fastq -o Cadam_hap1_VM Cadam_hap1_chromosomes.fasta
python /path/to/VerityMap/veritymap/main.py -t 20 -d hifi-diploid --reads Cadam.hifi.fastq -o Cadam_hap2_VM Cadam_hap2_chromosomes.fasta
```
#### Minimap2
We used [minimap2](https://github.com/lh3/minimap2) to check the read coverage in specific regions. Additionally, we filtered multi-mapped and low-quality reads using [samtools](https://github.com/samtools/samtools) to check the coverage of unique and high-quality mapped reads im those specific regions.
```
minimap2 -t 40 -ax map-hifi Cadam_primary_chromosomes.fasta Cadam.hifi.fastq | samtools view -@ 40 -bS -o Cadam_primary_hifi.bam -
samtools sort -@ 40 Cadam_primary_hifi.bam -o Cadam_primary_hifi.mapped.bam
samtools index Cadam_primary_hifi.bam
samtools view -@ 40 -q 30 -b Cadam_primary_hifi.mapped.bam -o Cadam_primary_hifi.mapped.q30.bam

minimap2 -t 40 -ax map-hifi Cadam_hap1_chromosomes.fasta Cadam.hifi.fastq | samtools view -@ 40 -bS -o Cadam_hap1_hifi.bam -
samtools sort -@ 40 Cadam_hap1_hifi.bam -o Cadam_hap1_hifi.mapped.bam
samtools index Cadam_hap1_hifi.bam
samtools view -@ 40 -q 30 -b Cadam_hap1_hifi.mapped.bam -o Cadam_hap1_hifi.mapped.q30.bam

minimap2 -t 40 -ax map-hifi Cadam_hap2_chromosomes.fasta Cadam.hifi.fastq | samtools view -@ 40 -bS -o Cadam_hap2_hifi.bam -
samtools sort -@ 40 Cadam_hap2_hifi.bam -o Cadam_hap2_hifi.mapped.bam
samtools index Cadam_hap2_hifi.bam
samtools view -@ 40 -q 30 -b Cadam_hap2_hifi.mapped.bam -o Cadam_hap2_hifi.mapped.q30.bam
```

#### NucFreq
We used [NucFreq](https://github.com/mrvollger/NucFreq) to check for collapsed regions, which may reflect assembly artifacts occurring in error-prone regions.
```
#retrieve a bam file for the specific region
samtools view -h -o ${chr}_${begin}-${end}.bam mapped.bam ${chr}:${begin}-${end}
samtools index ${chr}_${begin}-${end}.bam
#run nucplot
NucPlot.py ${chr}_${begin}-${end}.bam ${chr}_${begin}-${end}.png -t $thread
```

#### Plotting coverage of specific regions
We used [pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks) to plot coverage of mapped hifi reads using minimap2, Inspector, and VerityMap. We followed the well-detailed documentation of pyGenomeTracks.

The bam files were converted into bedGraph to plot charts: ```bedtools genomecov -ibam mapped.bam -bg > mapped.bedGraph```

Regions screened:
- **SVMP**
  - Primary - Cadam_mi-1:26000000-26900000
  - Haplotype_1 - Cadam_hap1_mi-1:20900000-21700000
  - Haplotype_2 - Cadam_hap2_mi-1:22200000-23100000

- **SVSP**
  - Primary - Cadam_mi-2:8900000-10200000
  - Haplotype_1 - Cadam_hap1_mi-2:8800000-9600000
  - Haplotype_2 - Cadam_hap2_mi-2:8800000-9900000
