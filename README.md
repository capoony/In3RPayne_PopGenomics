## Bioinformatic pipeline for population genomic inference of _In(3R)Payne_

see Material and Methods in [Kapun _et al._ (2023)]() for more details

### 1) Mapping sequencing data

For each of the newly sequenced library, test the quality with FASTQC, trim the raw reads with cutadapt, map the reads with bbmap, sort and deduplicate the BAM file with picard and realign around InDels with GATK

```bash
## loop through raw reads in folder /data and store BAM in /data/mapping
for i in <samplenames>
do

sh	shell/mapping.sh \
/data/${i}_R1.fq.gz	\
/data/${i}_R2.fq.gz	\	
${i} \
/data/mapping/${i} \
bbmap

done 

```