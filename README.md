# Bioinformatic pipeline for population genomic inference of _In(3R)Payne_

see Material and Methods in [Kapun _et al._ (2023)]() for more details

## 1) Map and Generae Phased Data

### 1.1) USA

For each of the newly sequenced library in the USA, test the quality with FASTQC, trim the raw reads with cutadapt, map the reads with bbmap, sort and deduplicate the BAM file with picard and realign around InDels with GATK

```bash
## loop through raw reads in folder /data and store BAM in /data/mapping
for file in /data/USA/*_R1.fq.gz
do

    tmp=${file##*/}
    i=${tmp%%_*}

    sh shell/mapping.sh \
        /data/USA/${i}_R1.fq.gz	\
        /data/USA/${i}_R2.fq.gz	\	
        ${i} \
        /data/USA/mapping/${i} \
        bbmap

done 
```

Following the approach in [Kapun _et al._ (2014)]() bioninformatically obtain haploid genomes from hemiclones

```bash

## synchronize as MPILEUP based on hologenome from Kapun et al. (2020) only including 3L and 3R
samtools mpileup \
    -B \
    -f reference/Dmel_6.04_hologenome_v2.fasta \
    -b /data/USA_BAM.txt \
    -l data /regions.bed.txt > /data/USA/USA.mpileup

## use PoPoolation2 (Kofler et al. 2011) to covert MPILEUP to SYNC format
java -jar /scripts/popoolation2_1201/mpileup2sync.jar \
    --input /data/USA/USA.mpileup \
    --threads 20 \
    --output /data/USA/USA.sync

## calculate max coverage threshold (see Kapun et al. 2014)
python2.7 /scripts/max_coverage.py \
    --input /data/USA/USA.sync \
    --max-coverage 0.05 \
    --output /data/USA/USA.cov

## use GNU parallel to parallelize the phasing (see Kapun et al. 2014 for details)
parallel -a /data/USA/USA.sync \
    -k \
    --pipepart \
    --cat python2.7 /scripts/extract_consensus.py \
        --input {} \
        --min-coverage 10 \
        --min-count 10 \
        --max-coverage /data/USA/USA.cov \
        --output /data/consensus/usa_min10_max005_mc10 \
        | gzip > /data/consensus/usa_min10_max005_mc10.consensus.gz

 ## convert consensus to VCF and retain site if > 50% of all samples with non-N's
python2.7 /scripts/cons2vcf.py \
    --input /data/consensus/usa_min10_max005_mc10.consensus.gz \
    --output /data/consensus/usa_min10_max005_mc10 \
    --ind 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58 \
    --N-cutoff 0.5 \
    --names Florida_S_1142,Florida_S_1145,Florida_S_1153,Florida_S_1155,Florida_S_1156,Florida_S_1157,Florida_S_1163,Florida_S_1164,Florida_S_1167,Florida_S_1218,Florida_S_1189,Florida_S_1170,Florida_S_1178,Florida_S_1203,Florida_S_1204,Florida_S_1158,Florida_S_1149,Florida_S_1174,Florida_S_1160,Florida_I_1153,Florida_I_1165,Florida_I_1169,Florida_I_1203,Florida_I_1218,Florida_I_1142,Florida_I_1146,Florida_I_1147,Florida_I_1149,Florida_I_1150,Florida_I_1178,Florida_I_1143,Florida_I_1156,Florida_I_1160,Florida_I_1161,Florida_I_1162,Florida_I_1164,Florida_I_1174,Florida_I_1152,Florida_I_1158,Maine_S_10-96,Maine_S_10-95,Maine_S_10-82,Maine_S_10-53,Maine_S_10-73,Maine_S_10-24,Maine_S_10-72,Maine_S_10-12,Maine_S_10-77,Maine_S_10-89,Maine_S_10-76,Maine_S_10-69,Maine_S_10-93,Maine_S_10-57,Maine_S_10-58,Maine_S_10-60,Maine_S_10-67,Maine_S_10-84,Maine_S_10-79,Maine_S_10-81
```

### 1.2) Zambia from DGN dataset (see Lack, _et al._ [2016]())

```bash
## download raw data from SRA and map with same pipeline as above
while IFS=',' read -r name SRA
do 
    sh /shell/obtain-n-map-africa.sh \
        /data/Zambia/
        ${name} \
        ${SRA}
done < /data/Zambia_SRA.txt

## now merge all Zambian strains to a big mpileup
samtools mpileup \
    -B \
    -f /reference/Dmel_6.04_hologenome_v2.fasta \
    -b /data/Zambia_BAM.txt 
    -l /data/regions.bed \
    | gzip > /data/Zambia/afr.mpileup.gz

## convert to cons fileformat given that these libraries should be from haploid embryos (see Kapopoulou et al. 2020 for more details)
gunzip -c /data/Zambia/afr.mpileup.gz \
    | parallel \
    --pipe  \
    -k \
    --cat python2.7 /scripts/mpileup2cons.py \
        -c 10,200 \
        -m {} \
        -u 0.9 \
        -b 20 | gzip > /data/consensus/afr.consensus.gz

```

### 1.3) Portugal (see Kapun, _et al._ [2014]() and Franssen, _et al._[2016]())

```bash

## get and map data 
while IFS=',' read -r name FWD REV
do 
    sh /shell/obtain-n-map-portugal.sh \
        /data/Portugal/
        ${name} \
        ${FWD} \
        ${REV}
done < /data/Portugal_SRA.txt

## now merge all Portugese strains to a big mpileup
samtools mpileup \
    -B \
    -f /reference/Dmel_6.04_hologenome_v2.fasta \
    -b /data/Portugal_BAM.txt 
    -l /data/regions.bed \
    | gzip > /data/Portugal/Portugal.mpileup.gz


## use PoPoolation2 (Kofler et al. 2011) to covert MPILEUP to SYNC format
java -jar /scripts/popoolation2_1201/mpileup2sync.jar \
    --input /data/Portugal/Portugal.mpileup \
    --threads 20 \
    --output /data/Portugal/Portugal.sync

## calculate max coverage threshold (see Kapun et al. 2014)
python2.7 /scripts/max_coverage.py \
    --input /data/Portugal/Portugal.sync \
    --max-coverage 0.05 \
    --output /data/Portugal/Portugal.cov

## use GNU parallel to parallelize the phasing (see Kapun et al. 2014 for details)
parallel -a /data/Portugal/Portugal.sync \
    -k \
    --pipepart \
    --cat python2.7 /scripts/extract_consensus.py \
        --input {} \
        --min-coverage 10 \
        --min-count 10 \
        --max-coverage /data/Portugal/Portugal.cov \
        --output /data/consensus/portugal_min10_max005_mc10 \
        | gzip > /data/consensus/portugal_min10_max005_mc10.consensus.gz

```
### 1.4) Australia 

We downloaded the raw sequencing data of [Rane _et al._ 2015]() from SRA (accession: [PRJNA221876](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA221876/)) and then classified the karyotpyes as follows:

```bash

## loop through raw reads in folder /data and store BAM in /data/mapping
for file in /data/Australia/*_R1.fq.gz
do

    tmp=${file##*/}
    i=${tmp%%_*}

    sh shell/mapping.sh \
        /data/Australia/${i}_1.fq.gz	\
        /data/Australia/${i}_2.fq.gz	\	
        ${i} \
        /data/Australia/mapping/${i} \
        bbmap

done 

## synchronize as MPILEUP based on hologenome from Kapun et al. (2020) only including 3L and 3R
samtools mpileup \
    -B \
    -f reference/Dmel_6.04_hologenome_v2.fasta \
    -b /data/Australia_BAM.txt \
    -l data /regions.bed.txt > /data/Australia/Australia.mpileup

## use PoPoolation2 (Kofler et al. 2011) to covert MPILEUP to SYNC format
java -jar /scripts/popoolation2_1201/mpileup2sync.jar \
    --input /data/Australia/Australia.mpileup \
    --threads 20 \
    --output /data/Australia/Australia.sync

```
Then we used the dataset of In(3R)Payne-specific marker SNPs from [Kapun _et al._ 2014]() to test for the frequency of inversion-specific alleles in the individuals



Then we further downloaded SNP data from from DataDryad [doi:10.5061/dryad.5q0m8.](https://datadryad.org/stash/dataset/doi:10.5061/dryad.5q0m8)

## References

Kapun, M., Schalkwyk, H. van, McAllister, B., Flatt, T. & Schlötterer, C. 2014. Inference of chromosomal inversion dynamics from Pool-Seq data in natural and laboratory populations of Drosophila melanogaster. Molecular Ecology 23: 1813–1827.

Kapun, M., Barrón, M.G., Staubach, F., Obbard, D.J., Wiberg, R.A.W., Vieira, J., et al. 2020. Genomic Analysis of European Drosophila melanogaster Populations Reveals Longitudinal Structure, Continent-Wide Selection, and Previously Unknown DNA Viruses. Mol Biol Evol 37: 2661–2678. Oxford Academic.

Kofler, R., Pandey, R.V. & Schlotterer, C. 2011. PoPoolation2: identifying differentiation between populations using sequencing of pooled DNA samples (Pool-Seq). Bioinformatics 27: 3435–3436.
