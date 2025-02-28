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
We used [Juicer](https://github.com/aidenlab/juicer) to manually review the scaffolded genome following the standard [DNA Genome Assembly Cookbook](https://aidenlab.org/assembly/manual 180322.pdf) instructions.

```

```
To help in the review process, we also checked for telomeric sequences within scaffolds using [tidk](https://github.com/tolkit/telomeric-identifier)
```
```

After reviewing, we generate the final scaffolded primary assembly.
```
```

### Haplotype-resolved
The haplotype-resolved assemblies were generated using [RagTag](https://github.com/malonge/RagTag) with each haplotype contigs as a query and the primary chromosome-level assembly as a reference.


### Checking genome quality
#### BUSCO

```
busco -i Cadam_primary_chromosomes.fasta -m genome -l tetrapoda_odb10 -c 20 -o busco_primary
busco -i Cadam_hap1_chromosomes.fasta -m genome -l tetrapoda_odb10 -c 20 -o busco_hap1
busco -i Cadam_hap2_chromosomes.fasta -m genome -l tetrapoda_odb10 -c 20 -o busco_hap2
```

#### Inspector

#### VerityMap
