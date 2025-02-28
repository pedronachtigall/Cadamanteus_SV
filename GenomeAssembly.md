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
cat SRR21092035_HiFi_cutadapt.fastq SRR21092036_HiFi_cutadapt.fastq > Cadam.hifi.fastq

#trim hic reads
trim_galore --paired --phred33 --length 75 -q 25 --stringency 1 -e 0.1 -o SRR28357486_hic_tg SRR28357486_R1.fastq.gz SRR28357486_R2.fastq.gz
```

### Mitochondrial genome assembly
We used [MITGARD](https://github.com/pedronachtigall/MITGARD) to perform the mitogenome assembly. The available mitochondrial genome of *C. adamanteus* () was used as the reference.
```

```

### Draft genome assembly
We used [hifiasm](https://github.com/chhylp123/hifiasm) to perform the draft genome assembly using HiFi and HiC reads to acquire both resolved haplotypes.

### Scaffold genome
We used [YaHS](https://github.com/c-zhou/yahs) to scaffold the primary genome assembled by hifiasm.
