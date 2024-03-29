# Bioinformatic pipeline for population genomic inference of _In(3R)Payne_

see Material and Methods in [Kapun _et al._ (2023)](https://www.biorxiv.org/content/10.1101/2023.01.31.526462v1) for more details. Most scripts were originally written in Python v.2. They were lifted over to Python v.3 using the [`2to3`](https://docs.python.org/3/library/2to3.html) tool and storded in the /scripts/v3 folder. The original files can be found in the /scripts folder. 

## 1) Map and generate phased sequencing data

### 1.1) USA

For each of the newly sequenced library in the USA, test the quality with FASTQC, trim the raw reads with cutadapt, map the reads with bbmap, sort and deduplicate the BAM file with picard and realign around InDels with GATK

```bash
mkdir /data/mapping

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

Following the approach in [Kapun _et al._ (2014)](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12594) we bioinformatically obtained haploid genomes from hemiclones. This requires [PoPoolation2](https:/sourceforge.net/p/popoolation2/wiki/Main/) from ([Kofler et al. 2011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3232374/))

```bash
mkdir /data/USA
mkdir /data/consensus

## synchronize as MPILEUP based on hologenome from Kapun et al. (2020) only including 3L and 3R
samtools mpileup \
    -B \
    -f /data/Dmel_6.04_hologenome_v2.fasta \
    -b /data/USA_BAM.txt \
    -l /data/regions.bed.txt > /data/USA/USA.mpileup

## use PoPoolation2 (Kofler et al. 2011) to covert MPILEUP to SYNC format
java -jar /scripts/popoolation2_1201/mpileup2sync.jar \
    -input /data/USA/USA.mpileup \
    -threads 20 \
    -output /data/USA/USA.sync

## calculate max coverage threshold (see Kapun et al. 2014)
python3 /scripts/v3/max_coverage.py \
    -input /data/USA/USA.sync \
    -max-coverage 0.05 \
    -output /data/USA/USA.cov

## use GNU parallel to parallelize the phasing (see Kapun et al. 2014 for details)
parallel -a /data/USA/USA.sync \
    -k \
    -pipepart \
    -cat python3 /scripts/v3/extract_consensus.py \
        -input {} \
        -min-coverage 10 \
        -min-count 10 \
        -max-coverage /data/USA/USA.cov \
        -output /data/consensus/USA \
        | gzip > /data/consensus/USA.consensus.gz

 ## convert consensus to VCF and retain site if > 50% of all samples with non-N's
python3 /scripts/v3/cons2vcf.py \
    -input /data/consensus/USA.consensus.gz \
    -output /data/consensus/America \
    -ind 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58 \
    -N-cutoff 0.5 \
    -names Florida_S_1142,Florida_S_1145,Florida_S_1153,Florida_S_1155,Florida_S_1156,Florida_S_1157,Florida_S_1163,Florida_S_1164,Florida_S_1167,Florida_S_1218,Florida_S_1189,Florida_S_1170,Florida_S_1178,Florida_S_1203,Florida_S_1204,Florida_S_1158,Florida_S_1149,Florida_S_1174,Florida_S_1160,Florida_I_1153,Florida_I_1165,Florida_I_1169,Florida_I_1203,Florida_I_1218,Florida_I_1142,Florida_I_1146,Florida_I_1147,Florida_I_1149,Florida_I_1150,Florida_I_1178,Florida_I_1143,Florida_I_1156,Florida_I_1160,Florida_I_1161,Florida_I_1162,Florida_I_1164,Florida_I_1174,Florida_I_1152,Florida_I_1158,Maine_S_10-96,Maine_S_10-95,Maine_S_10-82,Maine_S_10-53,Maine_S_10-73,Maine_S_10-24,Maine_S_10-72,Maine_S_10-12,Maine_S_10-77,Maine_S_10-89,Maine_S_10-76,Maine_S_10-69,Maine_S_10-93,Maine_S_10-57,Maine_S_10-58,Maine_S_10-60,Maine_S_10-67,Maine_S_10-84,Maine_S_10-79,Maine_S_10-81
```

### 1.2) Worldwide samples from DGN dataset (see Lack, _et al._ [2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5100052/))

```bash

mkdir /data/DGNs

## download raw data from SRA and map with same pipeline as above
while IFS=',' read -r name SRA
do 
    sh /shell/obtain-n-map-africa.sh \
        /data/DGN/
        ${name} \
        ${SRA}
done < /data/DGN_SRA.txt

## now merge all strains to a big mpileup
samtools mpileup \
    -B \
    -f /data/Dmel_6.04_hologenome_v2.fasta \
    -b /data/DGN_BAM.txt 
    -l /data/regions.bed \
    | gzip > /data/DGN/DGN.mpileup.gz

## convert major allele to cons fileformat given that these libraries should be from haploid embryos (see Kapopoulou et al. 2020 for more details)
gunzip -c /data/DGN/DGN.mpileup.gz \
    | parallel \
    -pipe  \
    -k \
    -cat python3 /scripts/v3/mpileup2cons.py \
        -c 10,200 \
        -m {} \
        -u 0.9 \
        -b 20 | gzip > /data/consensus/DGN.consensus.gz

```

### 1.3) Zambia from DGN dataset (see Lack, _et al._ [2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5100052/))

```bash

mkdir /data/Zambia

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
    -f /data/Dmel_6.04_hologenome_v2.fasta \
    -b /data/Zambia_BAM.txt 
    -l /data/regions.bed \
    | gzip > /data/Zambia/Zambia.mpileup.gz

## convert to cons fileformat given that these libraries should be from haploid embryos (see Kapopoulou et al. 2020 for more details)
gunzip -c /data/Zambia/Zambia.mpileup.gz \
    | parallel \
    -pipe  \
    -k \
    -cat python3 /scripts/v3/mpileup2cons.py \
        -c 10,200 \
        -m {} \
        -u 0.9 \
        -b 20 | gzip > /data/consensus/Zambia.consensus.gz

## convert consensus to VCF and retain site if > 50% of all samples with non-N's
python3 /scripts/v3/cons2vcf.py \
    -input /data/consensus/Zambia.consensus.gz \
    -output /data/consensus/Africa \
    -ind 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162 \
    -N-cutoff 0.5 \
    -names ZI81,ZI508,ZI505,ZI486,ZI448,ZI446,ZI405,ZI397N,ZI384,ZI362,ZI341,ZI319,ZI293,ZI292,ZI291,ZI264,ZI237,ZI233,ZI228,ZI226,ZI221,ZI220,ZI99,ZI91,ZI86,ZI85,ZI76,ZI68,ZI61,ZI530,ZI527,ZI524,ZI517,ZI514N,ZI50N,ZI504,ZI491,ZI490,ZI488,ZI477,ZI472,ZI471,ZI468,ZI467,ZI466,ZI460,ZI458,ZI457,ZI456,ZI455N,ZI453,ZI447,ZI445,ZI444,ZI433,ZI431,ZI421,ZI418N,ZI413,ZI402,ZI398,ZI395,ZI394N,ZI392,ZI386,ZI378,ZI377,ZI374,ZI373,ZI370,ZI365,ZI364,ZI359,ZI358,ZI357N,ZI348,ZI344,ZI342,ZI339,ZI336,ZI335,ZI332,ZI33,ZI329,ZI324,ZI321,ZI320,ZI31N,ZI316,ZI314,ZI313,ZI311N,ZI303,ZI295,ZI286,ZI281,ZI28,ZI279,ZI276,ZI271,ZI27,ZI269,ZI268,ZI265,ZI263,ZI261,ZI26,ZI255,ZI254N,ZI253,ZI252,ZI251N,ZI250,ZI240,ZI239,ZI235,ZI232,ZI231,ZI225,ZI219,ZI216N,ZI213,ZI212,ZI211,ZI210,ZI207,ZI206,ZI200,ZI198,ZI196,ZI194,ZI193,ZI191,ZI185,ZI184,ZI182,ZI181,ZI179,ZI178,ZI176,ZI173,ZI172,ZI170,ZI167,ZI165,ZI164,ZI161,ZI157,ZI152,ZI138,ZI136,ZI134N,ZI129,ZI126,ZI118N,ZI117,ZI114N,ZI104,ZI10,ZI56,ZI420,ZI388
```

### 1.4) Portugal (see Kapun, _et al._ [2014](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12594) and Franssen, _et al._ [2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5100052/))

```bash

mkdir /data/Portugal

## get and map data from Franssen et al. 2014
while IFS=',' read -r name FWD REV
do 
    sh /shell/obtain-n-map-portugal.sh \
        /data/Portugal /
        ${name} \
        ${FWD} \
        ${REV}
done < /data/Portugal_SRA.txt

## get and map data from Kapun et al. 2016
while IFS=',' read -r name SRA
do 
    sh /shell/obtain-n-map-africa.sh \
        /data/Portugal /
        ${name} \
        ${SRA}
done < /data/Portugal_SRA2.txt

## now merge all Portugese strains to a big mpileup
samtools mpileup \
    -B \
    -f /data/Dmel_6.04_hologenome_v2.fasta \
    -b /data/Portugal_BAM.txt 
    -l /data/regions.bed \
    | gzip > /data/Portugal/Portugal.mpileup.gz

## use PoPoolation2 (Kofler et al. 2011) to covert MPILEUP to SYNC format
java -jar /scripts/popoolation2_1201/mpileup2sync.jar \
    -input /data/Portugal/Portugal.mpileup \
    -threads 20 \
    -output /data/Portugal/Portugal.sync

## calculate max coverage threshold (see Kapun et al. 2014)
python3 /scripts/v3/max_coverage.py \
    -input /data/Portugal/Portugal.sync \
    -max-coverage 0.05 \
    -output /data/Portugal/Portugal.cov

## use GNU parallel to parallelize the phasing (see Kapun et al. 2014 for details)
parallel -a /data/Portugal/Portugal.sync \
    -k \
    -pipepart \
    -cat python3 /scripts/v3/extract_consensus.py \
        -input {} \
        -min-coverage 10 \
        -min-count 10 \
        -max-coverage /data/Portugal/Portugal.cov \
        -output /data/consensus/Portugal \
        | gzip > /data/consensus/Portugal.consensus.gz

## convert consensus to VCF and retain site if > 50% of all samples with non-N's
python3 /scripts/v3/cons2vcf.py \
    -input /data/consensus/Portugal.consensus.gz \
    -output /data/consensus/Portugal \
    -ind 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 \
    -N-cutoff 0.5 \
    -names F0-8-1,F0-11-1,F0-16-1,F0-21-1,F0-24-1,F0-25-1,hot-168-1,cold-21-0,cold-89-0,cold-91-0,F0-2-0,F0-3-0,F0-4-0,F0-5-0,F0-6-0,F0-7-0,F0-9-0,F0-10-0,F0-12-0,F0-13-0,F0-14-0,F0-15-0,F0-17-0,F0-18-0,F0-27-0,F0-29-0
```

### 1.5) Sweden (see Kapopoulou, _et al._ [2020](https://www.nature.com/articles/s41598-020-79720-1))

Use [Aspera](https:/www.ibm.com/aspera/connect/) and [sratoolkit](https:/www.ncbi.nlm.nih.gov/sra/docs/sradownload/) to download and convert raw data 

```bash

mkdir /data/Sweden

## get and map data 
while IFS=$' \t\n'
read -r A B C D E F SRA H name remainder
do

## used the ASPERA app to download from SRA 
Aspera\ Connect.app/Contents/Resources/ascp \
    -i Aspera\ Connect.app/Contents/Resources/asperaweb_id_dsa.openssh \
    -T anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${SRA:0:6}/$SRA/$SRA.sra \
    /data/Sweden/

## convert to FASTQ
scripts/sratoolkit.2.9.2-mac64/bin/fastq-dump.2.5.2 \
    -gzip \
    -split-3 \
    -outdir /data/Sweden/ \
    -A $name \
    /data/Sweden/$SRA.sra

rm /data/Sweden/$SRA.sra

mkdir -p /data/Sweden/$name

## trim with cutadapt (min base quality >18, minium trimmed readlength 75bp) and map with bbmap (standard parameters) against D. mel v. 6.04
sh shell/mapping.sh \
    /data/Sweden/${name}_1.fastq.gz \
    /data/Sweden/${name}_2.fastq.gz \
    $name \
    /data/Sweden/$name \
    bbmap

done < /data/Sweden/Sweden_SRA.txt

## now merge all Portugese strains to a big mpileup
samtools mpileup \
    -B \
    -f /data/Dmel_6.04_hologenome_v2.fasta \
    -b /data/Sweden_BAM.txt 
    -l /data/regions.bed \
    | gzip > /data/Sweden/Sweden.mpileup.gz

## convert to cons fileformat given that these libraries should be from haploid embryos (see Kapopoulou et al. 2020 for more details)
gunzip -c /data/Sweden/Sweden.mpileup.gz \
    | parallel \
    -pipe  \
    -k \
    -cat python3 /scripts/v3/mpileup2cons.py \
        -c 10,200 \
        -m {} \
        -u 0.9 \
        -b 20 | gzip > /data/consensus/Sweden.consensus.gz

## Now merge data from Portugal and Sweden
python3 /scripts/merge_consensus.py \
    -consensus /data/consensus/Portugal.consensus.gz,/data/consensus/Sweden.consensus.gz \
    | gzip > /data/consensus/Europe.cons.gz

## convert consensus to VCF and retain site if > 50% of all samples with non-N's
python3 /scripts/v3/cons2vcf.py \
    -input /data/consensus/Europe.cons.gz \
    -output /data/consensus/Europe \
    -ind 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40 \
    -N-cutoff 0.5 \
    -names F0-8-1,F0-11-1,F0-16-1,F0-21-1,F0-24-1,F0-25-1,hot-168-1,cold-21-0,cold-89-0,cold-91-0,F0-2-0,F0-3-0,F0-4-0,F0-5-0,F0-6-0,F0-7-0,F0-9-0,F0-10-0,F0-12-0,F0-13-0,F0-14-0,F0-15-0,F0-17-0,F0-18-0,F0-27-0,F0-29-0,SU02n,SU05n,SU07n,SU08,SU21n,SU25n,SU26n,SU29,SU37n,SU58n,SU75n,SU81n,SU93n,SU94
```

### 1.6) Australia 

We downloaded the raw sequencing data of [Rane _et al._ 2015](https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.13161) from SRA (accession: [PRJNA221876](https:/www.ncbi.nlm.nih.gov/bioproject/PRJNA221876/)) and then classified the karyotpyes as follows:

```bash

mkdir /data/Australia

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
    -f /data/Dmel_6.04_hologenome_v2.fasta \
    -b /data/Australia_BAM.txt \
    -l /data/regions.bed.txt > /data/Australia/Australia.mpileup

## use PoPoolation2 (Kofler et al. 2011) to covert MPILEUP to SYNC format
java -jar /scripts/popoolation2_1201/mpileup2sync.jar \
    -input /data/Australia/Australia.mpileup \
    -threads 20 \
    -output /data/Australia/Australia.sync

```
Then we used the dataset of In(3R)Payne-specific marker SNPs `data/fixed095.txt` from [Kapun _et al._ 2014](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12594) to test for the frequency of inversion-specific alleles in the individuals

```bash

## identify the average frequency of inversion-specific alleles in the individual libraries
python3 /scripts/v3/Australia-karyos.py \
    -karyo /data/fixed095.txt \
    -input /data/Australia/Australia.sync \
    -names IiR-1,IiR-2,IiR-3,IiR-4,IiR-5,IiR-6,IiR-7,IiR-8,IiR-9,IiR-10,IiR-11,IiR-12,IiR-13,IiR-14,IiR-15,IiR-16,IiR-17,IiR-18,IiR-19,IsR-1,IsR-2,IsR-3,IsR-4,IsR-5,IsR-6,IsR-7,IsR-8,IsR-9,IsR-10,IsR-11,IsR-12,IsR-13,IsR-14,IsR-15,IsR-16,IsR-17,IsR-18,YSR-1,YSR-2,YSR-3,YSR-4,YSR-5,YSR-6,YSR-7,YSR-8,YSR-9,YSR-10,YSR-11,YSR-12,YSR-13,YSR-14,YSR-15,YSR-16,YSR-17,YSR-18 \
    > /data/Australia/Australia_karyos.txt

```
We only considered samples where the inversion-specific alleles at four SNPs showed an average frequency of >90% (considered "inverted") or <10% (considered "standard").

Then, we further downloaded the imputed SNP data from from DataDryad [doi:10.5061/dryad.5q0m8.](https:/datadryad.org/stash/dataset/doi:10.5061/dryad.5q0m8) and used the karyotypic classification from above to calculate poppulation genetic statistics for each group separately.

## 2) Phylogenetic approach 

Here, we investigated the phylogenetic signal of worldwide samples based on SNPs located within _In(3R)P_ or in distance to the inversion.

### 2.1) Adjust Australian samples to match the other data

At first, we lifted genomic coordinates of the data, originally mapped to the _D. melanogaster_ reference genome v. 5 to v.6 coordinates.

```bash
## convert coordinates
python3 /scripts/v3/Flybase_version_converter.py \
    -input /data/Australia/3L-3R_Imputed_Merged_Variant_calls/3L_IiIsYs.vcf \
    -output /data/Australia/3L-3R_Imputed_Merged_Variant_calls/3L_IiIsYs_v6.vcf

python3 /scripts/v3/Flybase_version_converter.py \
    -input /data/Australia/3L-3R_Imputed_Merged_Variant_calls/3R_IiIsYs.vcf \
    -output /data/Australia/3L-3R_Imputed_Merged_Variant_calls/3R_IiIsYs_v6.vcf

cat /data/Australia/3L-3R_Imputed_Merged_Variant_calls/3L_IiIsYs_v6.vcf \
/data/Australia/3L-3R_Imputed_Merged_Variant_calls/3R_IiIsYs_v6.vcf \
> /data/consensus/Australia.vcf

## then convert from VCF to cons format
python3 /scripts/v3/vcf2cons.py \
    /data/consensus/Australia.vcf \
    | gzip > /data/Australia/Australia.cons.gz

```

### 2.2) Merge all files in *.cons format

Now, we merge the cons files and either draw randomly selected SNPs inside _In(3R)Payne_ or which are in distance to the Inversion boundaries

```bash
## Merge and select Inv SNPs
python3 /scripts/merge_consensus.py \
    - /data/Inv.bed \
    -consensus /data/consensus/Zambia.consensus.gz,/data/consensus/DGN.consensus.gz,/data/consensus/USA.consensus.gz,/data/consensus/Portugal.consensus.gz,/data/consensus/Sweden.consensus.gz,/data/Australia/Australia.cons.gz \
    | gzip > /data/consensus/Inv.cons.gz

## Merge and select Non-Inv SNPs
python3 /scripts/merge_consensus.py \
    - /data/NoInv.bed \
    -consensus /data/consensus/Zambia.consensus.gz,/data/consensus/DGN.consensus.gz,/data/consensus/USA.consensus.gz,/data/consensus/Portugal.consensus.gz,/data/consensus/Sweden.consensus.gz,/data/Australia/Australia.cons.gz \
    | gzip > /data/consensus/NoInv.cons.gz

## convert merged cons to nexus format
python3 /scripts/v3/cons2nexus.py \
    -input /data/consensus/Inv.cons.gz \
    -population 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483 \
    -names ZI81,ZI508,ZI505,ZI486,ZI448,ZI446,ZI405,ZI397N,ZI384,ZI362,ZI341,ZI319,ZI293,ZI292,ZI291,ZI264,ZI237,ZI233,ZI228,ZI226,ZI221,ZI220,ZI99,ZI91,ZI86,ZI85,ZI76,ZI68,ZI61,ZI530,ZI527,ZI524,ZI517,ZI514N,ZI50N,ZI504,ZI491,ZI490,ZI488,ZI477,ZI472,ZI471,ZI468,ZI467,ZI466,ZI460,ZI458,ZI457,ZI456,ZI455N,ZI453,ZI447,ZI445,ZI444,ZI433,ZI431,ZI421,ZI418N,ZI413,ZI402,ZI398,ZI395,ZI394N,ZI392,ZI386,ZI378,ZI377,ZI374,ZI373,ZI370,ZI365,ZI364,ZI359,ZI358,ZI357N,ZI348,ZI344,ZI342,ZI339,ZI336,ZI335,ZI332,ZI33,ZI329,ZI324,ZI321,ZI320,ZI31N,ZI316,ZI314,ZI313,ZI311N,ZI303,ZI295,ZI286,ZI281,ZI28,ZI279,ZI276,ZI271,ZI27,ZI269,ZI268,ZI265,ZI263,ZI261,ZI26,ZI255,ZI254N,ZI253,ZI252,ZI251N,ZI250,ZI240,ZI239,ZI235,ZI232,ZI231,ZI225,ZI219,ZI216N,ZI213,ZI212,ZI211,ZI210,ZI207,ZI206,ZI200,ZI198,ZI196,ZI194,ZI193,ZI191,ZI185,ZI184,ZI182,ZI181,ZI179,ZI178,ZI176,ZI173,ZI172,ZI170,ZI167,ZI165,ZI164,ZI161,ZI157,ZI152,ZI138,ZI136,ZI134N,ZI129,ZI126,ZI118N,ZI117,ZI114N,ZI104,ZI10,ZI56,ZI420,ZI388,FR180,FR217,FR229,FR361,CK1,CO10N,CO13N,CO14,CO15N,CO16,CO1,CO2,CO4N,CO8N,CO9N,EZ25,EZ2,EZ9N,GA125,GA129,GA132,GA141,GA145,GA160,GA185,GU10,GU6,GU7,GU9,KN133N,KN20N,KN35,KR39,KR42,KR4N,KR7,KT6,NG10N,NG1N,NG3N,NG9,SP221,TZ10,TZ14,TZ8,UM118,UM37,UM526,ZO65,ZS11,ZS37,ZS5,CK2,ED2,ED3,ED5N,ED6N,ED10N,EZ5N,GA130,GA191,GU2,KN6,KN34,KT1,NG6N,NG7,RC1,RC5,RG2,RG3,RG4N,RG5,RG6N,RG7,RG8,RG9,RG10,RG11N,RG13N,RG15,RG18N,RG19,RG21N,RG22,RG24,RG25,RG28,RG32N,RG33,RG34,RG35,RG36,RG37N,RG38N,RG39,SP80,SP173,SP188,SP235,SP241,SP254,UG5N,UG7,UG19,UG28N,ZI91a,ZI261a,ZI268,ZI468,ZL130,FR14,FR70,FR151,FR310,RAL-301,RAL-303,RAL-304,RAL-306,RAL-307,RAL-313,RAL-315,RAL-324,RAL-335,RAL-357,RAL-358,RAL-360,RAL-362,RAL-365,RAL-375,RAL-379,RAL-380,RAL-391,RAL-399,RAL-427,RAL-437,RAL-486,RAL-513,RAL-517,RAL-555,RAL-639,RAL-705,RAL-707,RAL-714,RAL-730,RAL-732,RAL-765,RAL-774,RAL-786,RAL-799,RAL-820,RAL-852,Florida_S_1142,Florida_S_1145,Florida_S_1153,Florida_S_1155,Florida_S_1156,Florida_S_1157,Florida_S_1163,Florida_S_1164,Florida_S_1167,Florida_S_1218,Florida_S_1189,Florida_S_1170,Florida_S_1178,Florida_S_1203,Florida_S_1204,Florida_S_1158,Florida_S_1149,Florida_S_1174,Florida_S_1160,Florida_I_1153,Florida_I_1165,Florida_I_1169,Florida_I_1203,Florida_I_1218,Florida_I_1142,Florida_I_1146,Florida_I_1147,Florida_I_1149,Florida_I_1150,Florida_I_1178,Florida_I_1143,Florida_I_1156,Florida_I_1160,Florida_I_1161,Florida_I_1162,Florida_I_1164,Florida_I_1174,Florida_I_1152,Florida_I_1158,Maine_S_10-96,Maine_S_10-95,Maine_S_10-82,Maine_S_10-53,Maine_S_10-73,Maine_S_10-24,Maine_S_10-72,Maine_S_10-12,Maine_S_10-77,Maine_S_10-89,Maine_S_10-76,Maine_S_10-69,Maine_S_10-93,Maine_S_10-57,Maine_S_10-58,Maine_S_10-60,Maine_S_10-67,Maine_S_10-84,Maine_S_10-79,Maine_S_10-81,F0-8-1,F0-11-1,F0-16-1,F0-21-1,F0-24-1,F0-25-1,hot-168-1,cold-21-0,cold-89-0,cold-91-0,F0-2-0,F0-3-0,F0-4-0,F0-5-0,F0-6-0,F0-7-0,F0-9-0,F0-10-0,F0-12-0,F0-13-0,F0-14-0,F0-15-0,F0-17-0,F0-18-0,F0-27-0,F0-29-0,168_hot_R2,150_hot_R2,143_hot_R2,136_hot_R1,129_hot_R1,117_hot_R1,106_hot_R1,100_cold_R1,96_cold_R1,91_cold_R1,89_cold_R2,80_cold_R2,53_cold_R2,52_cold_R2,21_cold_R3,SU02n,SU05n,SU07n,SU08,SU21n,SU25n,SU26n,SU29,SU37n,SU58n,SU75n,SU81n,SU93n,SU94,IiR-1,IiR-2,IiR-3,IiR-4,IiR-5,IiR-6,IiR-7,IiR-8,IiR-9,IiR-10,IiR-11,IiR-12,IiR-13,IiR-14,IiR-15,IiR-16,IiR-17,IiR-18,IiR-19,IsR-1,IsR-2,IsR-3,IsR-4,IsR-5,IsR-6,IsR-7,IsR-8,IsR-9,IsR-10,IsR-11,IsR-12,IsR-13,IsR-14,IsR-15,IsR-16,IsR-17,IsR-18,YSR-1,YSR-2,YSR-3,YSR-4,YSR-5,YSR-6,YSR-7,YSR-8,YSR-9,YSR-10,YSR-11,YSR-12,YSR-13,YSR-14,YSR-15,YSR-16,YSR-17,YSR-18 \
    > /data/consensus/Inv.nexus

python3 /scripts/v3/cons2nexus.py \
    -input /data/consensus/NoInv.cons.gz \
    -population 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401,402,403,404,405,406,407,408,409,410,411,412,413,414,415,416,417,418,419,420,421,422,423,424,425,426,427,428,429,430,431,432,433,434,435,436,437,438,439,440,441,442,443,444,445,446,447,448,449,450,451,452,453,454,455,456,457,458,459,460,461,462,463,464,465,466,467,468,469,470,471,472,473,474,475,476,477,478,479,480,481,482,483 \
    -names ZI81,ZI508,ZI505,ZI486,ZI448,ZI446,ZI405,ZI397N,ZI384,ZI362,ZI341,ZI319,ZI293,ZI292,ZI291,ZI264,ZI237,ZI233,ZI228,ZI226,ZI221,ZI220,ZI99,ZI91,ZI86,ZI85,ZI76,ZI68,ZI61,ZI530,ZI527,ZI524,ZI517,ZI514N,ZI50N,ZI504,ZI491,ZI490,ZI488,ZI477,ZI472,ZI471,ZI468,ZI467,ZI466,ZI460,ZI458,ZI457,ZI456,ZI455N,ZI453,ZI447,ZI445,ZI444,ZI433,ZI431,ZI421,ZI418N,ZI413,ZI402,ZI398,ZI395,ZI394N,ZI392,ZI386,ZI378,ZI377,ZI374,ZI373,ZI370,ZI365,ZI364,ZI359,ZI358,ZI357N,ZI348,ZI344,ZI342,ZI339,ZI336,ZI335,ZI332,ZI33,ZI329,ZI324,ZI321,ZI320,ZI31N,ZI316,ZI314,ZI313,ZI311N,ZI303,ZI295,ZI286,ZI281,ZI28,ZI279,ZI276,ZI271,ZI27,ZI269,ZI268,ZI265,ZI263,ZI261,ZI26,ZI255,ZI254N,ZI253,ZI252,ZI251N,ZI250,ZI240,ZI239,ZI235,ZI232,ZI231,ZI225,ZI219,ZI216N,ZI213,ZI212,ZI211,ZI210,ZI207,ZI206,ZI200,ZI198,ZI196,ZI194,ZI193,ZI191,ZI185,ZI184,ZI182,ZI181,ZI179,ZI178,ZI176,ZI173,ZI172,ZI170,ZI167,ZI165,ZI164,ZI161,ZI157,ZI152,ZI138,ZI136,ZI134N,ZI129,ZI126,ZI118N,ZI117,ZI114N,ZI104,ZI10,ZI56,ZI420,ZI388,FR180,FR217,FR229,FR361,CK1,CO10N,CO13N,CO14,CO15N,CO16,CO1,CO2,CO4N,CO8N,CO9N,EZ25,EZ2,EZ9N,GA125,GA129,GA132,GA141,GA145,GA160,GA185,GU10,GU6,GU7,GU9,KN133N,KN20N,KN35,KR39,KR42,KR4N,KR7,KT6,NG10N,NG1N,NG3N,NG9,SP221,TZ10,TZ14,TZ8,UM118,UM37,UM526,ZO65,ZS11,ZS37,ZS5,CK2,ED2,ED3,ED5N,ED6N,ED10N,EZ5N,GA130,GA191,GU2,KN6,KN34,KT1,NG6N,NG7,RC1,RC5,RG2,RG3,RG4N,RG5,RG6N,RG7,RG8,RG9,RG10,RG11N,RG13N,RG15,RG18N,RG19,RG21N,RG22,RG24,RG25,RG28,RG32N,RG33,RG34,RG35,RG36,RG37N,RG38N,RG39,SP80,SP173,SP188,SP235,SP241,SP254,UG5N,UG7,UG19,UG28N,ZI91a,ZI261a,ZI268,ZI468,ZL130,FR14,FR70,FR151,FR310,RAL-301,RAL-303,RAL-304,RAL-306,RAL-307,RAL-313,RAL-315,RAL-324,RAL-335,RAL-357,RAL-358,RAL-360,RAL-362,RAL-365,RAL-375,RAL-379,RAL-380,RAL-391,RAL-399,RAL-427,RAL-437,RAL-486,RAL-513,RAL-517,RAL-555,RAL-639,RAL-705,RAL-707,RAL-714,RAL-730,RAL-732,RAL-765,RAL-774,RAL-786,RAL-799,RAL-820,RAL-852,Florida_S_1142,Florida_S_1145,Florida_S_1153,Florida_S_1155,Florida_S_1156,Florida_S_1157,Florida_S_1163,Florida_S_1164,Florida_S_1167,Florida_S_1218,Florida_S_1189,Florida_S_1170,Florida_S_1178,Florida_S_1203,Florida_S_1204,Florida_S_1158,Florida_S_1149,Florida_S_1174,Florida_S_1160,Florida_I_1153,Florida_I_1165,Florida_I_1169,Florida_I_1203,Florida_I_1218,Florida_I_1142,Florida_I_1146,Florida_I_1147,Florida_I_1149,Florida_I_1150,Florida_I_1178,Florida_I_1143,Florida_I_1156,Florida_I_1160,Florida_I_1161,Florida_I_1162,Florida_I_1164,Florida_I_1174,Florida_I_1152,Florida_I_1158,Maine_S_10-96,Maine_S_10-95,Maine_S_10-82,Maine_S_10-53,Maine_S_10-73,Maine_S_10-24,Maine_S_10-72,Maine_S_10-12,Maine_S_10-77,Maine_S_10-89,Maine_S_10-76,Maine_S_10-69,Maine_S_10-93,Maine_S_10-57,Maine_S_10-58,Maine_S_10-60,Maine_S_10-67,Maine_S_10-84,Maine_S_10-79,Maine_S_10-81,F0-8-1,F0-11-1,F0-16-1,F0-21-1,F0-24-1,F0-25-1,hot-168-1,cold-21-0,cold-89-0,cold-91-0,F0-2-0,F0-3-0,F0-4-0,F0-5-0,F0-6-0,F0-7-0,F0-9-0,F0-10-0,F0-12-0,F0-13-0,F0-14-0,F0-15-0,F0-17-0,F0-18-0,F0-27-0,F0-29-0,168_hot_R2,150_hot_R2,143_hot_R2,136_hot_R1,129_hot_R1,117_hot_R1,106_hot_R1,100_cold_R1,96_cold_R1,91_cold_R1,89_cold_R2,80_cold_R2,53_cold_R2,52_cold_R2,21_cold_R3,SU02n,SU05n,SU07n,SU08,SU21n,SU25n,SU26n,SU29,SU37n,SU58n,SU75n,SU81n,SU93n,SU94,IiR-1,IiR-2,IiR-3,IiR-4,IiR-5,IiR-6,IiR-7,IiR-8,IiR-9,IiR-10,IiR-11,IiR-12,IiR-13,IiR-14,IiR-15,IiR-16,IiR-17,IiR-18,IiR-19,IsR-1,IsR-2,IsR-3,IsR-4,IsR-5,IsR-6,IsR-7,IsR-8,IsR-9,IsR-10,IsR-11,IsR-12,IsR-13,IsR-14,IsR-15,IsR-16,IsR-17,IsR-18,YSR-1,YSR-2,YSR-3,YSR-4,YSR-5,YSR-6,YSR-7,YSR-8,YSR-9,YSR-10,YSR-11,YSR-12,YSR-13,YSR-14,YSR-15,YSR-16,YSR-17,YSR-18 \
    > /data/consensus/NoInv.nexus
```

Using the NEXUS files, we plotted phylogenetic networks based on the Neighbor-Net inference method with Splits-Tree.

## 3) Population Genetics Analysis 

### 3.1) Genetic variation

Using [VCFtools](https:/vcftools.github.io/index.html) we first calculated nucleotide diveristy (_π_) and Tajima's _D_ as measures of genetic diversity in the different karyotypes

```bash

## make directory 
mkdir /data/PopGen

## loop through continents/countries
for continent in Africa Europe America Australia
do
    ## loop through karyotypes
    for karyotype in Inv Std All
    do
        ## loop through regions
        for file in data/${continent}_*_${karyotype}.samples
        do      
            ## isolate region name
            IFS='_ ' read -r -A array <<< "$file"
            region=${array[2]}

            ## calculate π
            /scripts/vcftools_0.1.13/bin/vcftools \
                -keep data/${continent}_${region}_${karyotype}.samples \
                -window-pi 100000 \
                -window-pi-step 100000 \
                -gzvcf /data/consensus/${continent}.vcf.gz \
                -out /data/PopGen/${continent}_${region}_${karyotype}_100k.pi

            ## calculate Tajima's D
            /scripts/vcftools_0.1.13/bin/vcftools \
                -keep data/${continent}_${region}_${karyotype}.samples \
                -TajimaD 100000 \
                -gzvcf /data/consensus/${continent}.vcf.gz \
                -out /data/PopGen/${continent}_${region}_${karyotype}_100k
        
        done
    done

## then merge the results in a single large file for π and Tajima's D, plot averages by karyotype and Geographic origin and test for significant differences.
python3 /scripts/v3/merge-div.py \
    -input /data/PopGen/Africa_Zambia_Inv_100k.pi,/data/PopGen/Africa_Zambia_Std_100k.pi,/data/PopGen/Africa_Zambia_All_100k.pi,/data/PopGen/Europe_Portugal_Inv_100k.pi,/data/PopGen/Europe_Portugal_Std_100k.pi,/data/PopGen/Europe_Portugal_All_100k.pi,/data/PopGen/Europe_Sweden_Std_100k.pi,/data/PopGen/America_Florida_Inv_100k.pi,/data/PopGen/America_Florida_Std_100k.pi,/data/PopGen/America_Florida_All_100k.pi,/data/PopGen/America_Maine_Std_100k.pi,/data/PopGen/Australia_Queensland_Inv_100k.pi,/data/PopGen/Australia_Queensland_Std_100k.pi,/data/PopGen/Australia_Queensland_All_100k.pi,/data/PopGen/Australia_Victoria_Std_100k.pi \
    -names Zambia-Inv,Zambia-Std,Zambia-All,Portugal-Inv,Portugal-Std,Portugal-All,Sweden-Std,Florida-Inv,Florida-Std,Florida-All,Maine-Std,Queensland-Inv,Queensland-Std,Queensland-All,Victoria-Std \
    > /data/PopGen/AllData.pi

Rscript /scripts/statNplot.R /data/PopGen/AllData.pi

python3 /scripts/v3/merge-div.py \
    -input /data/PopGen/Africa_Zambia_Inv_100k.Tajima.D,/data/PopGen/Africa_Zambia_Std_100k.Tajima.D,/data/PopGen/Africa_Zambia_All_100k.Tajima.D,/data/PopGen/Europe_Portugal_Inv_100k.Tajima.D,/data/PopGen/Europe_Portugal_Std_100k.Tajima.D,/data/PopGen/Europe_Portugal_All_100k.Tajima.D,/data/PopGen/Europe_Sweden_Std_100k.Tajima.D,/data/PopGen/America_Florida_Inv_100k.Tajima.D,/data/PopGen/America_Florida_Std_100k.Tajima.D,/data/PopGen/America_Florida_All_100k.Tajima.D,/data/PopGen/America_Maine_Std_100k.Tajima.D,/data/PopGen/Australia_Queensland_Inv_100k.Tajima.D,/data/PopGen/Australia_Queensland_Std_100k.Tajima.D,/data/PopGen/Australia_Queensland_All_100k.Tajima.D,/data/PopGen/Australia_Victoria_Std_100k.Tajima.D \
    -names Zambia-Inv,Zambia-Std,Zambia-All,Portugal-Inv,Portugal-Std,Portugal-All,Sweden-Std,Florida-Inv,Florida-Std,Florida-All,Maine-Std,Queensland-Inv,Queensland-Std,Queensland-All,Victoria-Std \
    > /data/PopGen/AllData.Tajima.D

Rscript /scripts/statNplot.R /data/PopGen/AllData.Tajima.D
```
In addition, I plotted the values of _π_ and Tajima's _D_ as lineplots in R

```R
### plot π for Inverted and Standard Chromosomes
library(tidyverse)

DATA=read.table("/data/PopGen/AllData.pi",header=T)

DATA.3R=subset(DATA,DATA$C=="3R")
DATA.3R=subset(DATA,DATA$InvStatus!="All")
DATA.3R=subset(DATA.3R,DATA.3R$Origin %in% c("Florida","Portugal","Zambia","Queensland"))
DATA.3R$Origin <- factor(DATA.3R$Origin,levels=c("Zambia","Portugal","Florida","Queensland"))

P=ggplot(DATA.3R,aes(
    x=P,
    y=pi,
    col=InvStatus))+
  geom_line(lwd=1.5)+
  geom_rect(mapping = aes(
        xmin=16432209,
        xmax=24744010,
        ymin=0,
        ymax=(max(pi))),
    color="black",
    alpha=0,
    lwd=0.5,
    lty=2)+
  facet_grid(Origin~ .)+
  theme_classic()

ggsave("/data/PopGen/AllData.pi.pdf",
    P,
    width=12,
    height=8)

### plot Tajima's D for Inverted and Standard Chromosomes
library(tidyverse)

DATA=read.table("/data/PopGen/AllData.Tajima.D",header=T)

DATA.3R=subset(DATA,DATA$C=="3R")
DATA.3R=subset(DATA.3R,DATA.3R$Origin %in% c("Florida","Portugal","Zambia","Queensland"))
DATA.3R$Origin <- factor(DATA.3R$Origin,levels=c("Zambia","Portugal","Florida","Queensland"))

P=ggplot(DATA.3R,aes(
    x=P,
    y=pi,
    col=InvStatus))+
  geom_line(lwd=1.5)+
  geom_rect(mapping = aes(
        xmin=16432209,
        xmax=24744010,
        ymin=0,
        ymax=(max(pi))),
    color="black",
    alpha=0,
    lwd=0.5,
    lty=2)+
  facet_grid(Origin~ .)+
  theme_classic()+
  geom_hline(yintercept=0)

ggsave("/data/PopGen/AllData.Tajima.D.pdf",
    P,
    width=12,
    height=8)
```

### 3.2) Genetic differentiation

Again using [VCFtools](https:/vcftools.github.io/index.html) we calculated pairwise FST among karyotypes in different geographic regions both for single SNPs and in 100kb windows. Then, we merged the FST values, plotted the averages with respect to karyotype and tested for significant differences with respect to the karyotype and origin. Moreover, we tested if window-wise FST is correlated within the inverted genomic region.

```bash

## loop through samples
continent=( "Africa" "Europe" "Europe" "Europe" "America" "America" "America" "Australia" "Australia" "Australia" )
FST=( "ZI-ZS" "PI-PS" "PI-SS" "PS-SS" "FI-FS" "FI-MS" "FS-MS" "II-IS" "II-YS" "IS-YS" )
Pop1=( "Africa_Zambia_Inv" "Europe_Portugal_Inv" "Europe_Portugal_Inv" "Europe_Portugal_Std" "America_Florida_Inv" "America_Florida_Inv" "America_Florida_Std" "Australia_Queensland_Inv" "Australia_Queensland_Inv" "Australia_Queensland_Std" )
Pop2=( "Africa_Zambia_Std" "Europe_Portugal_Std" "Europe_Sweden_Std" "Europe_Sweden_Std" "America_Florida_Std" "America_Maine_Std" "America_Maine_Std" "Australia_Queensland_Std" "Australia_Victoria_Std" "Australia_Victoria_Std" )

for index in ${!FST[@]}
do
    ## window-wise FST
    /scripts/vcftools_0.1.13/cpp/vcftools \
        -gzvcf /data/consensus/${continent[index]}.vcf.gz \
        -weir-fst-pop /data/${Pop1[index]}.samples \
        -weir-fst-pop /data/${Pop2[index]}.samples \
        -fst-window-size 100000 \
        -fst-window-step 100000 \
        -out /data/PopGen/${FST[index]}
    
    ## SNP-wise FST
    /scripts/vcftools_0.1.13/cpp/vcftools \
        -gzvcf /data/consensus/${continent[index]}.vcf.gz \
        -weir-fst-pop /data/${Pop1[index]}.samples \
        -weir-fst-pop /data/${Pop2[index]}.samples \
        -out /data/PopGen/${FST[index]}
done

## now merge all FST values
python3 /scripts/v3/merge-diff.py \
    -input /data/PopGen/ZI-ZS.windowed.weir.fst,/data/PopGen/PI-PS.windowed.weir.fst,/data/PopGen/PI-SS.windowed.weir.fst,/data/PopGen/PS-SS.windowed.weir.fst,/data/PopGen/FI-FS.windowed.weir.fst,/data/PopGen/FI-MS.windowed.weir.fst,/data/PopGen/FS-MS.windowed.weir.fst,/data/PopGen/II-IS.windowed.weir.fst,/data/PopGen/II-YS.windowed.weir.fst,/data/PopGen/IS-YS.windowed.weir.fst \
    -names ZI-ZS,PI-PS,PI-SS,PS-SS,FI-FS,FI-MS,FS-MS,II-IS,II-YS,IS-YS \
 > /data/PopGen/AllData.fst

## plot average FSTs per karyotype and geography, test for significant differences, and investigate if windowed FST within the inverted regions is correlated among geographic regions 
Rscript /scripts/statNcompareFST.R /data/PopGen/AllData.fst

```
Next, we also calculated SNP-wise FST and used [SnpEFF](http:/pcingola.github.io/SnpEff/) to annotate the VCF files.

```bash
## loop through continents/countries
for continent in Africa Europe America Australia
do   
    ## annotate with SNPeff
    java -Xmx4g -jar /scripts/snpEff/snpEff.jar \
        BDGP6.82 \
        /data/consensus/${continent}.vcf.gz \
        | gzip > /data/consensus/${continent}_ann.vcf.gz
done

## Isolate positions with FST >= 0.9 and append annotation
continent=( "Africa" "Europe" "Europe" "Europe" "America" "America" "America" "Australia" "Australia" "Australia" )
FST=( "ZI-ZS" "PI-PS" "PI-SS" "PS-SS" "FI-FS" "FI-MS" "FS-MS" "II-IS" "II-YS" "IS-YS" )

for index in ${!FST[@]}
do  
    ## append to full dataset
    python3 /scripts/v3/AppendAnnotFromVCF.py \
        /data/consensus/${continent[index]}_ann.vcf.gz \
        /data/PopGen/${FST[index]}.weir.fst \
        | sort -k1,1 -k2,2 \
        | > /data/PopGen/${FST[index]}_ann.weir.fst
    
    ## additionally isolate candidates
    awk '$NF>=0.9 && $NF!="nan"' /data/PopGen/${FST[index]}.weir.fst \
    | python3 /scripts/v3/AppendAnnotFromVCF.py \
        /data/consensus/${continent[index]}_ann.vcf.gz \
        - \
        | sort -k1,1 -k2,2 \
        | > /data/PopGen/${FST[index]}_cand_ann.weir.fst
done
```

Now test in R if there is a significant overlap among candidate SNPs and Genes across continents and plot the results as barplots using the [SuperExactTest](https://cran.r-project.org/web/packages/SuperExactTest/index.html) package

```R
## first SNP-wise comparison

library(SuperExactTest)

## read background and idenitfy SNPs inside the inversion
SEu.f=na.omit(read.table(gzfile('/data/PopGen/PI-PS_ann.weir.fst'),header=F,na.string='nan',sep='\t', quote = ''))
SEu.f.in3rp=subset(SEu.f,SEu.f[,1]=='3R' & SEu.f[,2]> 16432209 & SEu.f[,2] < 24744010)
SAf.f=na.omit(read.table(gzfile('/data/PopGen/ZI-ZS_ann.weir.fst'),header=F,na.string='nan',sep='\t', quote = ''))
SAf.f.in3rp=subset(SAf.f,SAf.f[,1]=='3R' & SAf.f[,2]> 16432209 & SAf.f[,2] < 24744010)
SAu.f=na.omit(read.table(gzfile('/data/PopGen/II-IS_ann.weir.fst'),header=F,na.string='nan',sep='\t', quote = ''))
SAu.f.in3rp=subset(SAu.f,SAu.f[,1]=='3R' & SAu.f[,2]> 16432209 & SAu.f[,2] < 24744010)
SNA.f=na.omit(read.table(gzfile('/data/PopGen/FI-FS_ann.weir.fst'),header=F,na.string='nan',sep='\t', quote = ''))
SNA.f.in3rp=subset(SNA.f,SNA.f[,1]=='3R' & SNA.f[,2]> 16432209 & SNA.f[,2] < 24744010)

## now for genes
SAf.f.genes=unique(SAf.f$V5[sapply(SAf.f$V5,FUN=function(x){x!=''})])
SEu.f.genes=unique(SEu.f$V5[sapply(SEu.f$V5,FUN=function(x){x!=''})])
SAu.f.genes=unique(SAu.f$V5[sapply(SAu.f$V5,FUN=function(x){x!=''})])
SNA.f.genes=unique(SNA.f$V5[sapply(SNA.f$V5,FUN=function(x){x!=''})])

SAf.f.in3rp.genes=unique(SAf.f.in3rp$V5[sapply(SAf.f.in3rp$V5,FUN=function(x){x!=''})])
SEu.f.in3rp.genes=unique(SEu.f.in3rp$V5[sapply(SEu.f.in3rp$V5,FUN=function(x){x!=''})])
SAu.f.in3rp.genes=unique(SAu.f.in3rp$V5[sapply(SAu.f.in3rp$V5,FUN=function(x){x!=''})])
SNA.f.in3rp.genes=unique(SNA.f.in3rp$V5[sapply(SNA.f.in3rp$V5,FUN=function(x){x!=''})])

## now for SNPs
SEu.f=unique(paste(SEu.f[,1],SEu.f[,2],sep='_'))
SEu.f.in3rp=unique(paste(SEu.f.in3rp[,1],SEu.f.in3rp[,2],sep='_'))
SAf.f=unique(paste(SAf.f[,1],SAf.f[,2],sep='_'))
SAf.f.in3rp=unique(paste(SAf.f.in3rp[,1],SAf.f.in3rp[,2],sep='_'))
SAu.f=unique(paste(SAu.f[,1],SAu.f[,2],sep='_'))
SAu.f.in3rp=unique(paste(SAu.f.in3rp[,1],SAu.f.in3rp[,2],sep='_'))
SNA.f=unique(paste(SNA.f[,1],SNA.f[,2],sep='_'))
SNA.f.in3rp=unique(paste(SNA.f.in3rp[,1],SNA.f.in3rp[,2],sep='_'))

# get overlapping SNPs & genes
FullnoAu=Reduce(intersect,list(SEu.f,SAf.f,SNA.f))
FullnoAu.in3p=Reduce(intersect,list(SEu.f.in3rp,SAf.f.in3rp,SNA.f.in3rp))

Full=Reduce(intersect,list(SEu.f,SAf.f,SNA.f,SAu.f))
Full.in3p=Reduce(intersect,list(SEu.f.in3rp,SAf.f.in3rp,SNA.f.in3rp,SAu.f.in3rp))

FullnoAu.genes=Reduce(intersect,list(SAf.f.genes,SEu.f.genes,SNA.f.genes))
FullnoAu.in3p.genes=Reduce(intersect,list(SEu.f.in3rp.genes,SAf.f.in3rp.genes,SNA.f.in3rp.genes))

Full.genes=Reduce(intersect,list(SEu.f.genes,SAf.f.genes,SNA.f.genes,SAu.f.genes))
Full.in3p.genes=Reduce(intersect,list(SEu.f.in3rp.genes,SAf.f.in3rp.genes,SNA.f.in3rp.genes,SAu.f.in3rp.genes))

dir.create('/data/PopGen/candidates')

# candidates
SEu.c=read.table(paste('/data/PopGen/PI-PS_cand_ann.weir.fst',sep=''),sep='\t', quote = '')
SEu.c.in3rp=subset(SEu.c,SEu.c[,1]=='3R' & SEu.c[,2]> 16432209 & SEu.c[,2] < 24744010)
SEu.c.genes=unique(SEu.c$V5[sapply(SEu.c$V5,FUN=function(x){x!=''})])
SEu.c.in3rp.genes=unique(SEu.c.in3rp$V5[sapply(SEu.c.in3rp$V5,FUN=function(x){x!=''})])
SEu.c=unique(paste(SEu.c[,1],SEu.c[,2],sep='_'))
SEu.c.in3rp=unique(paste(SEu.c.in3rp[,1],SEu.c.in3rp[,2],sep='_'))
SEu.c.ovNoAu.in3rp=SEu.c[SEu.c %in% FullnoAu.in3p]
SEu.c.ov.in3rp=SEu.c[SEu.c %in% Full.in3p]
SEu.c.ovNoAu=SEu.c[SEu.c %in% FullnoAu]
SEu.c.ov=SEu.c[SEu.c %in% Full]

SAf.c=read.table(paste('/data/PopGen/ZI-ZS_cand_ann.weir.fst',sep=''),sep='\t', quote = '')
SAf.c.in3rp=subset(SAf.c,SAf.c[,1]=='3R' & SAf.c[,2]> 16432209 & SAf.c[,2] < 24744010)
SAf.c.genes=unique(SAf.c$V5[sapply(SAf.c$V5,FUN=function(x){x!=''})])
SAf.c.in3rp.genes=unique(SAf.c.in3rp$V5[sapply(SAf.c.in3rp$V5,FUN=function(x){x!=''})])
SAf.c=unique(paste(SAf.c[,1],SAf.c[,2],sep='_'))
SAf.c.in3rp=paste(SAf.c.in3rp[,1],SAf.c.in3rp[,2],sep='_')
SAf.c.ovNoAu.in3rp=SAf.c[SAf.c %in% FullnoAu.in3p]
SAf.c.ov.in3rp=SAf.c[SAf.c %in% Full.in3p]
SAf.c.ovNoAu=SAf.c[SAf.c %in% FullnoAu]
SAf.c.ov=SAf.c[SAf.c %in% Full]

SAu.c=read.table(paste('/data/PopGen/II-IS_cand_ann.weir.fst',sep=''),sep='\t', quote = '')
SAu.c.in3rp=subset(SAu.c,SAu.c[,1]=='3R' & SAu.c[,2]> 16432209 & SAu.c[,2] < 24744010)
SAu.c.genes=unique(SAu.c$V5[sapply(SAu.c$V5,FUN=function(x){x!=''})])
SAu.c.in3rp.genes=unique(SAu.c.in3rp$V5[sapply(SAu.c.in3rp$V5,FUN=function(x){x!=''})])
SAu.c=paste(SAu.c[,1],SAu.c[,2],sep='_')
SAu.c.in3rp=paste(SAu.c.in3rp[,1],SAu.c.in3rp[,2],sep='_')
SAu.c.ov.in3rp=SAu.c[SAu.c %in% Full.in3p]
SAu.c.ov=SAu.c[SAu.c %in% Full]


SNA.c=read.table(paste('/data/PopGen/FI-FS_cand_ann.weir.fst',sep=''),sep='\t', quote = '')
SNA.c.in3rp=subset(SNA.c,SNA.c[,1]=='3R' & SNA.c[,2]> 16432209 & SNA.c[,2] < 24744010)
SNA.c.genes=unique(SNA.c$V5[sapply(SNA.c$V5,FUN=function(x){x!=''})])
SNA.c.in3rp.genes=unique(SNA.c.in3rp$V5[sapply(SNA.c.in3rp$V5,FUN=function(x){x!=''})])
SNA.c=paste(SNA.c[,1],SNA.c[,2],sep='_')
SNA.c.in3rp=paste(SNA.c.in3rp[,1],SNA.c.in3rp[,2],sep='_')
SNA.c.ovNoAu.in3rp=SNA.c[SNA.c %in% FullnoAu.in3p]
SNA.c.ov.in3rp=SNA.c[SNA.c %in% Full.in3p]
SNA.c.ovNoAu=SNA.c[SNA.c %in% FullnoAu]
SNA.c.ov=SNA.c[SNA.c %in% Full]


## overlap in candidate SNPs from all continents except Australia and in all genomic regions
sink(paste('/data/PopGen/candidates/All-noAu.txt',sep=''))
pdf(paste('/data/PopGen/candidates/All-noAu.pdf',sep=''),width=8,height=8)
sets=list('Europe'=unique(SEu.c.ovNoAu),'NorthAmerica'=unique(SNA.c.ovNoAu),'Africa'=unique(SAf.c.ovNoAu))
res=supertest(sets,n=length(FullnoAu))
plot(res,Layout='landscape',degree=2:4,sort.by = 'size',x.pos=c(0.1,0.9),minMinusLog10PValue=2,maxMinusLog10PValue=300)
dev.off()
print(summary(res))
sink()
write.table(data.frame('NAME'=summary(res)[1],'GENES'=summary(res)[8]),file=paste('/data/PopGen/candidates/All-noAU.list',sep=''),row.names=F,quote = F,sep='\t')

## overlap in candidate SNPs from all continents except Australia and in In(3R)P only
sink(paste('/data/PopGen/candidates/All-noAu-in3rpOnly.txt',sep=''))
pdf(paste('/data/PopGen/candidates/All-noAu-in3rpOnly.pdf',sep=''),width=8,height=8)
sets.in3rp=list('Europe'=SEu.c.ovNoAu.in3rp,'NorthAmerica'=SNA.c.ovNoAu.in3rp,'Africa'=SAf.c.ovNoAu.in3rp)
res.in3rp=supertest(sets.in3rp,n=length(FullnoAu.in3p))
plot(res.in3rp,Layout='landscape',degree=2:4,sort.by = 'size',x.pos=c(0.1,0.9),minMinusLog10PValue=2,maxMinusLog10PValue=300)
dev.off()
print(summary(res.in3rp))
sink()
write.table(data.frame('NAME'=summary(res.in3rp)[1],'GENES'=summary(res.in3rp)[8]),file=paste('/data/PopGen/candidates/All-noAu-in3rpOnly.list',sep=''),row.names=F,quote = F,sep='\t')

## overlap in candidate SNPs from all continents in In(3R)P only
sink(paste('/data/PopGen/candidates/All-in3rpOnly.txt',sep=''))
pdf(paste('/data/PopGen/candidates/All-in3rpOnly.pdf',sep=''),width=8,height=8)
sets.in3rp.AU=list('Europe'=SEu.c.ov.in3rp,'NorthAmerica'=SNA.c.ov.in3rp,'Africa'=SAf.c.ov.in3rp,'Australia'=SAu.c.ov.in3rp)
res.in3rp.AU=supertest(sets.in3rp.AU,n=length(Full.in3p))
plot(res.in3rp.AU,Layout='landscape',degree=2:8,sort.by = 'size',x.pos=c(0.1,0.9),minMinusLog10PValue=2,maxMinusLog10PValue=50)
dev.off()
print(summary(res.in3rp.AU))
sink()
write.table(data.frame('NAME'=summary(res.in3rp.AU)[1],'GENES'=summary(res.in3rp.AU)[8]),file=paste('/data/PopGen/candidates/All-in3rpOnly.list',sep=''),row.names=F,quote = F,sep='\t')

## overlap in candidate SNPs from all continents and in all genomic regions
sink(paste('/data/PopGen/candidates/All.txt',sep=''))
pdf(paste('/data/PopGen/candidates/All.pdf',sep=''),width=8,height=8)
sets.AU=list('Europe'=SEu.c.ov,'NorthAmerica'=SNA.c.ov,'Africa'=SAf.c.ov,'Australia'=SAu.c.ov)
res.AU=supertest(sets.AU,n=length(Full.in3p))
plot(res.AU,Layout='landscape',degree=2:8,sort.by = 'size',x.pos=c(0.1,0.9),minMinusLog10PValue=0,maxMinusLog10PValue=50)
dev.off()
print(summary(res.AU))
sink()
write.table(data.frame('NAME'=summary(res.AU)[1],'GENES'=summary(res.AU)[8]),file=paste('/data/PopGen/candidates/All.list',sep=''),row.names=F,quote = F,sep='\t')

# now gene-based

## get overlapping genes
SAf.c.ovNoAu.in3rp.genes=SAf.c.genes[SAf.c.genes %in% FullnoAu.in3p.genes]
SAf.c.ov.in3rp.genes=SAf.c.genes[SAf.c.genes %in% Full.in3p.genes]
SAf.c.ovNoAu.genes=SAf.c.genes[SAf.c.genes %in% FullnoAu.genes]
SAf.c.ov.genes=SAf.c.genes[SAf.c.genes %in% Full.genes]

SEu.c.ovNoAu.in3rp.genes=SEu.c.genes[SEu.c.genes %in% FullnoAu.in3p.genes]
SEu.c.ov.in3rp.genes=SEu.c.genes[SEu.c.genes %in% Full.in3p.genes]
SEu.c.ovNoAu.genes=SEu.c.genes[SEu.c.genes %in% FullnoAu.genes]
SEu.c.ov.genes=SEu.c.genes[SEu.c.genes %in% Full.genes]

SNA.c.ovNoAu.in3rp.genes=SNA.c.genes[SNA.c.genes %in% FullnoAu.in3p.genes]
SNA.c.ov.in3rp.genes=SNA.c.genes[SNA.c.genes %in% Full.in3p.genes]
SNA.c.ovNoAu.genes=SNA.c.genes[SNA.c.genes %in% FullnoAu.genes]
SNA.c.ov.genes=SNA.c.genes[SNA.c.genes %in% Full.genes]

SAu.c.ov.in3rp.genes=SAu.c.genes[SAu.c.genes %in% Full.in3p.genes]
SAu.c.ov.genes=SAu.c.genes[SAu.c.genes %in% Full.genes]

## overlap in candidate genes from all continents except Australia and in all genomic regions
sink(paste('/data/PopGen/candidates/All-noAu-genes.txt',sep=''))
pdf(paste('/data/PopGen/candidates/All-noAu-genes.pdf',sep=''),width=8,height=8)
sets=list('Europe'=unique(SEu.c.ovNoAu.genes),'NorthAmerica'=unique(SNA.c.ovNoAu.genes),'Africa'=unique(SAf.c.ovNoAu.genes))
res=supertest(sets,n=length(FullnoAu))
plot(res,Layout='landscape',degree=2:4,sort.by = 'size',x.pos=c(0.1,0.9),minMinusLog10PValue=2,maxMinusLog10PValue=300)
dev.off()
print(summary(res))
sink()
write.table(data.frame('NAME'=summary(res)[1],'GENES'=summary(res)[8]),file=paste('/data/PopGen/candidates/All-nAU-genes.list',sep=''),row.names=F,quote = F,sep='\t')

## overlap in candidate genes from all continents except Australia in In(3R)P only
sink(paste('/data/PopGen/candidates/All-noAu-in3rpOnly-genes.txt',sep=''))
pdf(paste('/data/PopGen/candidates/All-noAu-in3rpOnly-genes.pdf',sep=''),width=8,height=8)
sets.in3rp=list('Europe'=SEu.c.ovNoAu.in3rp.genes,'NorthAmerica'=SNA.c.ovNoAu.in3rp.genes,'Africa'=SAf.c.ovNoAu.in3rp.genes)
res.in3rp=supertest(sets.in3rp,n=length(FullnoAu.in3p))
plot(res.in3rp,Layout='landscape',degree=2:4,sort.by = 'size',x.pos=c(0.1,0.9),minMinusLog10PValue=2,maxMinusLog10PValue=300)
dev.off()
print(summary(res.in3rp))
sink()
write.table(data.frame('NAME'=summary(res.in3rp)[1],'GENES'=summary(res.in3rp)[8]),file=paste('/data/PopGen/candidates/All-noAU-in3rpOnly-genes.list',sep=''),row.names=F,quote = F,sep='\t')

## overlap in candidate genes from all continents in In(3R)P only
sink(paste('/data/PopGen/candidates/All-in3rpOnly-genes.txt',sep=''))
pdf(paste('/data/PopGen/candidates/All-in3rpOnly-genes.pdf',sep=''),width=8,height=8)
sets.in3rp.AU=list('Europe'=SEu.c.ov.in3rp.genes,'NorthAmerica'=SNA.c.ov.in3rp.genes,'Africa'=SAf.c.ov.in3rp.genes,'Australia'=SAu.c.ov.in3rp.genes)
res.in3rp.AU=supertest(sets.in3rp.AU,n=length(Full.in3p))
plot(res.in3rp.AU,Layout='landscape',degree=2:8,sort.by = 'size',x.pos=c(0.1,0.9),minMinusLog10PValue=2,maxMinusLog10PValue=200)
dev.off()
print(summary(res.in3rp.AU))
sink()
write.table(data.frame('NAME'=summary(res.in3rp.AU)[1],'GENES'=summary(res.in3rp.AU)[8]),file=paste('/data/PopGen/candidates/All-in3rpOnly-genes.list',sep=''),row.names=F,quote = F,sep='\t')

## overlap in candidate genes from all continents and all genomic regions
sink(paste('/data/PopGen/candidates/All-genes.txt',sep=''))
pdf(paste('/data/PopGen/candidates/All-genes.pdf',sep=''),width=8,height=8)
sets.AU=list('Europe'=SEu.c.ov.genes,'NorthAmerica'=SNA.c.ov.genes,'Africa'=SAf.c.ov.genes,'Australia'=SAu.c.ov.genes)
res.AU=supertest(sets.AU,n=length(Full.in3p))
plot(res.AU,Layout='landscape',degree=2:8,sort.by = 'size',x.pos=c(0.1,0.9),minMinusLog10PValue=0,maxMinusLog10PValue=200)
dev.off()
print(summary(res.AU))
sink()
write.table(data.frame('NAME'=summary(res.AU)[1],'GENES'=summary(res.AU)[8]),file=paste('/data/PopGen/candidates/All-genes.list',sep=''),row.names=F,quote = F,sep='\t')
```

### 3.3) African origin

Here, we tested if inversion-specific candidate alleles in Florida are also found in the standard arrangement in flies for the ancestral African origin.

```bash 

mkdir /data/InvAlleles

continent=( "Africa" "America" )
region=( "Zambia" "Florida"  )
Pop1=( "Africa_Zambia_Inv" "America_Florida_Inv" )
Pop2=( "Africa_Zambia_Std" "America_Florida_Std" )
FST=( "ZI-ZS" "FI-FS" )

for index in ${!FST[@]}
do
    ## Identify inversion-specific allele in candidates from specific region
    python3 /scripts/InvAf_from_VCF.py \
        --input /data/consensus/${continent[index]}.vcf.gz \
        --inv /data/${Pop1[index]}.samples \
        --std /data/${Pop2[index]}.samples \
        | python3 /scripts/v3/OverlapSNPs.py \
        --source /data/PopGen/${FST[index]}_cand_ann.weir.fst \
        --target - \
        > /data/InvAlleles/${continent[index]}.alleles
    
    Sample=( "Africa_Zambia_Inv" "Africa_Zambia_Std" "Europe_Portugal_Inv" "Europe_Portugal_Std" "Europe_Sweden_Std" "America_Florida_Inv" "America_Florida_Std" "America_Maine_Std" "Australia_Queensland_Inv" "Australia_Queensland_Std" "Australia_Victoria_Std" )

    for index2 in ${!Sample[@]}
    do
        ## get allele frequency of inversion-specific allele in subset of data
        python3 /scripts/AF_by_alleleNind.py \
            --input /data/consensus/${continent[index]}.vcf.gz \
            --AFl /data/InvAlleles/${continent[index]}.alleles \
            --pops /data/${Sample[index2]}.samples \
            > /data/InvAlleles/${continent[index]}_${Sample[index2]}.af

    done
done
```
Now plot with _R_ for inversion-specific alleles in North America

```R

FlS=read.table("/data/InvAlleles/America_America_Florida_S.af",header=F)
FlI=read.table("/data/InvAlleles/America_America_Florida_I.af",header=F)
MaS=read.table("/data/InvAlleles/America_America_Maine_S.af",header=F)
ZaI=read.table("/data/InvAlleles/America_Africa_Zambia_I.af",header=F)
ZaS=read.table("/data/InvAlleles/America_Africa_Zambia_S.af",header=F)
PoI=read.table("/data/InvAlleles/America_Europe_Portugal_I.af",header=F)
PoS=read.table("/data/InvAlleles/America_Europe_Portugal_S.af",header=F)
SwS=read.table("/data/InvAlleles/America_Europe_Sweden_S.af",header=F)

INTER<-Reduce(intersect,list(FlS$V2,FlI$V2,MaS$V2,ZaI$V2,ZaS$V2,PoI$V2,PoS$V2,SwS$V2))

NEW=list("Florida-I"=FlI[,4][FlI$V2 %in% INTER],"Florida-S"=FlS[,4][FlS$V2 %in% INTER],"Maine-S"=MaS[,4][MaS$V2 %in% INTER],"Zambia-I"=ZaI[,4][ZaI$V2 %in% INTER],"Zambia-S"=ZaS[,4][ZaS$V2 %in% INTER],"Portugal-I"=PoI[,4][PoI$V2 %in% INTER],"Portugal-S"=PoS[,4][PoS$V2 %in% INTER],"Sweden-S"=SwS[,4][SwS$V2 %in% INTER])

## plot AFs in all populations
pdf("/data/InvAlleles/Inversion-specific-AF.pdf",width=12,height=6)
boxplot(NEW,ylab="Allele Frequencies")
dev.off()

## add position information and write to file
NEW=data.frame(list("POS"=FlI[,2][FlI$V2 %in% INTER],"Florida-I"=FlI[,4][FlI$V2 %in% INTER],"Florida-S"=FlS[,4][FlS$V2 %in% INTER],"Maine-S"=MaS[,4][MaS$V2 %in% INTER],"Zambia-I"=ZaI[,4][ZaI$V2 %in% INTER],"Zambia-S"=ZaS[,4][ZaS$V2 %in% INTER],"Portugal-I"=PoI[,4][PoI$V2 %in% INTER],"Portugal-S"=PoS[,4][PoS$V2 %in% INTER],"Sweden-S"=SwS[,4][SwS$V2 %in% INTER]))

write.table(NEW,"/data/InvAlleles/Inversion-specific-AF.txt",quote=F,row.names=F)

## insolate positions with ABS(STD-INV)>0.5 in Africa
NEW2=NEW[abs(NEW[,5]-NEW[,6])>0.5,]
write.table(NEW2,"/data/InvAlleles/Inversion-specific-AF0.5.txt",quote=F,row.names=F)

H=hist(NEW2[,1],breaks=100,xlim=c(0,32000000))

# plot genomic positions of these SNPS
pdf("/data/InvAlleles/Inversion-specific-AF_dist.pdf",width=12,height=6)
plot(H$mids/1000000,H$counts,type="l",xlim=c(15,27),xlab="Genomic Position (Mbp)",ylab="No. of SNPs")
rect(16.432209,0,24.744010,max(H$counts),border="grey",lty=2)
dev.off()

## test if significant allelic difference between STD and INV in Africa
sink("/data/InvAlleles/Inversion-specific-AF.stat")
wilcox.test(NEW[,5],NEW[,6])
sink()
```


### 3.4) Linkage Disequilibrium

Here, we assessed how linkage disequilibrium decays within the genomic region spanned by _In(3R)Payne_ in inverted and non-inverted haplotypes. 

#### 3.4.1) SNP-wise LD with the inversion

First, we calculated and plotted LD of each SNP with _In(3R)Payne_ along the 3R chromosomal arm with the Inversion.
```bash

mkdir /data/LD

## Africa
python3 /scripts/v3/LD_bar.py \
--input /data/Africa/Africa.cons.gz \
--output /data/LD/Zambia_Bar \
--ind ZI81,ZI508,ZI505,ZI486,ZI448,ZI446,ZI405,ZI397N,ZI384,ZI362,ZI341,ZI319,ZI293,ZI292,ZI291,ZI264,ZI237,ZI233,ZI228,ZI226,ZI221,ZI220,ZI99,ZI91,ZI86,ZI85,ZI76,ZI68,ZI61,ZI530,ZI527,ZI524,ZI517,ZI514N,ZI50N,ZI504,ZI491,ZI490,ZI488,ZI477,ZI472,ZI471,ZI468,ZI467,ZI466,ZI460,ZI458,ZI457,ZI456,ZI455N,ZI453,ZI447,ZI445,ZI444,ZI433,ZI431,ZI421,ZI418N,ZI413,ZI402,ZI398,ZI395,ZI394N,ZI392,ZI386,ZI378,ZI377,ZI374,ZI373,ZI370,ZI365,ZI364,ZI359,ZI358,ZI357N,ZI348,ZI344,ZI342,ZI339,ZI336,ZI335,ZI332,ZI33,ZI329,ZI324,ZI321,ZI320,ZI31N,ZI316,ZI314,ZI313,ZI311N,ZI303,ZI295,ZI286,ZI281,ZI28,ZI279,ZI276,ZI271,ZI27,ZI269,ZI268,ZI265,ZI263,ZI261,ZI26,ZI255,ZI254N,ZI253,ZI252,ZI251N,ZI250,ZI240,ZI239,ZI235,ZI232,ZI231,ZI225,ZI219,ZI216N,ZI213,ZI212,ZI211,ZI210,ZI207,ZI206,ZI200,ZI198,ZI196,ZI194,ZI193,ZI191,ZI185,ZI184,ZI182,ZI181,ZI179,ZI178,ZI176,ZI173,ZI172,ZI170,ZI167,ZI165,ZI164,ZI161,ZI157,ZI152,ZI138,ZI136,ZI134N,ZI129,ZI126,ZI118N,ZI117,ZI114N,ZI104,ZI10,ZI56,ZI420,ZI388 \
--chromosome 3R \
--N-cutoff 0.1 \
--color blue,green,yellow,red,purple \
--min-allele 0.1 \
--Inv /data/Africa_Zambia_Inv.samples \
--Std /data/Africa_Zambia_Std.samples \

## Europe
python3 /scripts/v3/LD_bar.py \
--input /data/Europe/Europe.cons.gz \
--output /data/LD/Portugal_Bar \
--ind F0-8-1,F0-11-1,F0-16-1,F0-21-1,F0-24-1,F0-25-1,hot-168-1,cold-21-0,cold-89-0,cold-91-0,F0-2-0,F0-3-0,F0-4-0,F0-5-0,F0-6-0,F0-7-0,F0-9-0,F0-10-0,F0-12-0,F0-13-0,F0-14-0,F0-15-0,F0-17-0,F0-18-0,F0-27-0,F0-29-0,SU02n,SU05n,SU07n,SU08,SU21n,SU25n,SU26n,SU29,SU37n,SU58n,SU75n,SU81n,SU93n,SU94 \
--chromosome 3R \
--N-cutoff 0.1 \
--color blue,green,yellow,red,purple \
--min-allele 0.1 \
--Inv /data/Europe_Portugal_Inv.samples \
--Std /data/Europe_Portugal_Std.samples \

## America
python3 /scripts/v3/LD_bar.py \
--input /data/USA/USA.cons.gz \
--output /data/LD/USA_Bar \
--ind Florida_S_1142,Florida_S_1145,Florida_S_1153,Florida_S_1155,Florida_S_1156,Florida_S_1157,Florida_S_1163,Florida_S_1164,Florida_S_1167,Florida_S_1218,Florida_S_1189,Florida_S_1170,Florida_S_1178,Florida_S_1203,Florida_S_1204,Florida_S_1158,Florida_S_1149,Florida_S_1174,Florida_S_1160,Florida_I_1153,Florida_I_1165,Florida_I_1169,Florida_I_1203,Florida_I_1218,Florida_I_1142,Florida_I_1146,Florida_I_1147,Florida_I_1149,Florida_I_1150,Florida_I_1178,Florida_I_1143,Florida_I_1156,Florida_I_1160,Florida_I_1161,Florida_I_1162,Florida_I_1164,Florida_I_1174,Florida_I_1152,Florida_I_1158,Maine_S_10-96,Maine_S_10-95,Maine_S_10-82,Maine_S_10-53,Maine_S_10-73,Maine_S_10-24,Maine_S_10-72,Maine_S_10-12,Maine_S_10-77,Maine_S_10-89,Maine_S_10-76,Maine_S_10-69,Maine_S_10-93,Maine_S_10-57,Maine_S_10-58,Maine_S_10-60,Maine_S_10-67,Maine_S_10-84,Maine_S_10-79,Maine_S_10-81 \
--chromosome 3R \
--N-cutoff 0.1 \
--color blue,green,yellow,red,purple \
--min-allele 0.1 \
--Inv /data/America_Florida_Inv.samples \
--Std /data/America_Florida_Std.samples \

## Australia
python3 /scripts/v3/LD_bar.py \
--input /data/Australia/Australia.cons.gz \
--output /data/LD/Australia_Bar \
--ind Ii1,Ii2,Ii3,Ii4,Ii5,Ii6,Ii7,Ii8,Ii9,Ii10,Ii11,Ii12,Ii13,Ii14,Ii15,Ii16,Ii17,Ii18,Ii19,Is1,Is2,Is3,Is4,Is5,Is6,Is7,Is8,Is9,Is10,Is11,Is12,Is13,Is14,Is15,Is16,Is17,Is18,YS1,YS2,YS3,YS4,YS5,YS6,YS7,YS8,YS9,YS10,YS11,YS12,YS13,YS14,YS15,YS16,YS17,YS18 \
--chromosome 3R \
--N-cutoff 0.1 \
--color blue,green,yellow,red,purple \
--min-allele 0.1 \
--Inv /data/Australia_Queensland_Inv.samples \
--Std /data/Australia_Queensland_Std.samples \

```

#### 3.4.2) LD heatmaps

Next, we calculated LD heatmaps within an among karyotypes focusing on a subset of 5000 SNPs from populations with mixed karyotypes in different continents. The numbers of the parameter `ind` correspond to the (zero-based) position of the chosen individuals (from \*.samples files) in the consensus block. See [Kapun _et al._ 2014](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.12594) for more details. Below, we are showing the code for the LD analysis within karyotypes, the `ind` parameter needs to be adjusted (based on the \*.sample files in data/) when repeating the analyses for LD within karyotypes.

```bash

## Zambia - full
python3 /scripts/LD_heatmap.py \
    --input /data/Africa/Africa.cons.gz \
    --output /data/LD/Zambia \
    --subsample 5000 \
    --ind 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,31,32,33,51,54,72,74,76,90,94,105,114,116,119,120,126,135,137,146,151,158 \
    --chromosome 3R \
    --N-cutoff 0.1 \
    --color purple,red,yellow,green,blue \
    --min-allele 0.3

## Portugal - full 
python3 /scripts/LD_heatmap.py \
    --input /data/Europe/Europe.cons.gz \
    --output /data/LD/Portugal \
    --subsample 5000 \
    --ind 0,1,2,3,4,5,6,10,11,13,17,22,23,25 \
    --chromosome 3R \
    --N-cutoff 0.1 \
    --color purple,red,yellow,green,blue \
    --min-allele 0.3

## Florida - full 
python3 /scripts/LD_heatmap.py \
    --input /data/USA/USA.cons.gz \
    --output /data/LD/Florida \
    --subsample 5000 \
    --ind 0,1,2,3,4,5,6,7,8,9,10,11,13,14,15,17,18,19,21,22,23,24,25,26,27,29,30,31,32,33,34,35,36,37,38 \
    --chromosome 3R \
    --N-cutoff 0.1 \
    --color purple,red,yellow,green,blue \
    --min-allele 0.3

## Australia - full 
python3 /scripts/LD_heatmap.py \
    --input /data/Australia/Australia.cons.gz \
    --output /data/LD/Queensland \
    --subsample 5000 \
    --ind 0,1,3,6,9,10,11,15,16,17,25,32,35,36,12,13,14,19,21,24,26,27,29,31,33 \
    --chromosome 3R \
    --N-cutoff 0.1 \
    --color purple,red,yellow,green,blue \
    --min-allele 0.3
```

#### 3.4.3) Decay of LD

Next, we estimated differences in the decay of LD with respect to different karyotypes based on the `*.dist.gz` output files from the analyses above and tested for significant differences with non-linear mixed models.

```bash

## loop through samples
continent=( "Africa" "Europe" "America" "Australia" )
region=( "Zambia" "Portugal" "Florida" "Queensland" )
Pop1=( "Africa_Zambia_Inv" "Europe_Portugal_Inv"  "America_Florida_Inv"  "Australia_Queensland_Inv"  )
Pop2=( "Africa_Zambia_Std" "Europe_Portugal_Std" "America_Florida_Std" "Australia_Queensland_Std" )

for index in ${!Pop1[@]}
do  

    ## plot the decay
    python3 /scripts/v3/LD_decay.py \
        --vcf /data/consensus/${continent[index]}.vcf.gz \
        --samples /data/${Pop1[index]}.samples,/data/${Pop2[index]}.samples \
        --names I,S \
        --combsample 50000 \
        --SNPs 5000 \
        --maxdist 100000 \
        --threshold 0.1 \
        --N-threshold 0.5 \
        --region 3R:16432209-24744010 \
        --out /data/LD/${region[index]}_decay

    ## test for significant differences
    Rscript /scripts/NLS.R /data/LD/${region[index]}_decay

done
```

### 4) Transcriptomic analysis

#### 4.1) Trim and map RNASeq data 

```bash
## loop through all libraries and trim

mkdir -p /data/RNASeq/trimmed
for i in /data/RNASeq/*.fq.gz
do

    ## use only a subset of the variable i as the name
    name=${i##*/}
    n="$( cut -d '.' -f 1 <<< "$name" )"

    ## trim with cutadapt
    /scripts/cutadapt \
    -q 18 \
    --minimum-length 75 \
    -o /data/RNASeq/trimmed/$n.fq.gz \
    ${i} 

done
```
Get reference transcriptome

```bash
cd /data/RNASeq/
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.17_FB2017_04/fasta/dmel-all-CDS-r6.17.fasta.gz
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.17_FB2017_04/fasta/dmel-all-r6.17.gtf.gz
```

Next, we pseudomapped the trimmed reads with [Kallisto](https://github.com/pachterlab/kallisto)

```bash
## index the transcriptome
/scripts/kallisto index \
    -i /data/RNASeq/dmel-all-CDS-r6.17.idx \
    /data/RNASeq/dmel-all-CDS-r6.17.fasta.gz

mkdir /data/RNASeq/kallisto

## loop through all libraries and map with kallisto
for i in /data/RNASeq/*.fq.gz
do

    name=${i##*/}
    n="$( cut -d '.' -f 1 <<< "$name" )"

    kallisto quant \
    -i /data/RNASeq/dmel-all-CDS-r6.17.idx \
    -o /data/RNASeq/kallisto/$n \
    --single -l 101 -s 10 -b 100 --rf-stranded -t 24 \
    /data/RNASeq/trimmed/$name \

done
```
After that, we merged the read counts at gene level across all datasets

```bash

## merge datasets and sum up Pseudocounts for gene-level estimates
mkdir -p /data/RNASeq/count

python3 /scripts/v3/mergeKallisto.py \
    /data/RNASeq/kallisto \
    /data/RNASeq/count/analyses \
    /data/RNASeq/dmel-all-r6.17.gtf.gz
```

#### 4.2) Test for differential expression 

Then, we used [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) and [limma](https://bioconductor.org/packages/release/bioc/html/limma.html) to test for differential expression (DE) among different factors and test for enrichment of Gene Ontology (GO) terms.

```R
### Analysis of RNA-Seq data of inversion lines

## install packages
#install.packages("readr")
#install.packages("tximport")
#install.packages("edgeR")
#install.packages("RColorBrewer")
#install.packages("limma")
#install.packages("pheatmap")

I="/data/RNASeq/count/analyses_counts.txt"
O="/data/RNASeq/"

#Factors

Treat=factor(c('FI','FI','FI','FI','FI','FI','FI','FI','FI','FI','FI','FI','FI','FI','FI','FI','FI','FI','FS','FS','FS','FS','FS','FS','FS','FS','FS','FS','FS','FS','FS','FS','FS','FS','FS','FS','MS','MS','MS','MS','MS','MS','MS','MS','MS','MS','MS','MS','MS','MS','MS','MS','MS','MS'))
Treat=relevel(Treat,ref="FS")
Temp=factor(c('18','18','18','18','18','18','18','18','18','25','25','25','25','25','25','25','25','25','18','18','18','18','18','18','18','18','18','25','25','25','25','25','25','25','25','25','18','18','18','18','18','18','18','18','18','25','25','25','25','25','25','25','25','25'))
Rep=factor(c(10,14,17,18,22,23,26,5,9,5,9,10,14,17,18,22,23,26,13,16,2,21,27,3,4,6,7,2,3,4,6,7,13,16,21,27,1,11,12,15,19,20,24,25,8,1,8,11,12,15,19,20,24,25))

## 1a) read Data
rawData<-read.table(I,header=T)
rawCounts=rawData[,3:56]
rownames(rawCounts)<-rawData[,1]

## 2) create DGE object using edgeR
library(edgeR)
FullDGE<-DGEList(rawCounts,group=names(rawCounts))

## 3) add factors to object
FullDGE$samples$treat <- Treat
FullDGE$samples$temp <- Temp

## Transform counts (Counts per million: CPM) and normalize the counts using the TMM method in edgeR
FullDGE <- calcNormFactors(FullDGE, method = "TMM")
cpmFullDGE<-cpm(FullDGE,log=T)

out=paste(O,"/countdata/CPM.txt",sep="")
TAB=cpmFullDGE
write.table(data.frame("Gene"=rownames(TAB),TAB),file=out,row.names=F,quote = F)

## remove genes covered by less than 2 read per million in at least nine samples
keep <- rowSums(cpmFullDGE >= 2) >=9
FullDGE <- FullDGE[keep,] 
cpmTFullDGE <- cpmFullDGE[keep,] 

## Plot density histograms to quantify the effect of filtering out lowly expressed genes
dir.create(paste(O,"/quality",sep=""))
pdf(paste(O,"/quality/CPM-hist.pdf",sep=""),width=15,height=5)
nsamples <- ncol(FullDGE)
col <- topo.colors(nsamples)
par(mfrow=c(1,2))
plot(density(cpmFullDGE[,1]), col=col[1], lwd=2, ylim=c(0,0.25), las=2,
     main="", xlab="")
abline(v=0, lty=3)

for (i in 2:nsamples){
  den <- density(cpmFullDGE[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
} 
legend("topright", legend=colnames(FullDGE), text.col=col, bty="n",ncol = 3,cex=0.5)
title(main="A. Raw data", xlab="Log-cpm")

plot(density(cpmTFullDGE[,1]), col=col[1], lwd=2, ylim=c(0,0.25), las=2,
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm") 
abline(v=0, lty=3)

for (i in 2:nsamples){
  den <- density(cpmTFullDGE[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
} 
legend("topright", legend=colnames(FullDGE), text.col=col, bty="n",ncol = 3,cex=0.5)
dev.off()

## explore Data with MDS to identify clustering
library(RColorBrewer)
pdf(paste(O,"/quality/MDS.pdf",sep=""),width=14,height=12)
par(mfrow=c(2,2))

col.Treat <- Treat
levels(col.Treat) <- brewer.pal(nlevels(col.Treat), "Set1")
plotMDS(cpmTFullDGE, col=as.character(col.Treat))
title(main="Karypop: Dim 1/2")

col.Temp <- Temp
levels(col.Temp) <- brewer.pal(nlevels(col.Temp), "Set1")
plotMDS(cpmTFullDGE, col=as.character(col.Temp))
title(main="Temperature: Dim 1/2")

plotMDS(cpmTFullDGE, col=as.character(col.Treat), dim=c(3,4))
title(main="Karypop: Dim 3/4")

plotMDS(cpmTFullDGE, col=as.character(col.Temp), dim=c(3,4))
title(main="Temperature: Dim 3/4")

dev.off()

### explore data with PCA to explore clustering

pdf(paste(O,"/quality/PCA.pdf",sep=""),width=12,height=8)
par(mfrow=c(2,3))
newcpm=t(cpmTFullDGE)
PCA=prcomp(newcpm)
colKP=rainbow(length(levels(Treat)))[Treat]
colTemp=rainbow(length(levels(Temp)))[as.factor(Temp)]
colPop=rainbow(length(levels(Rep)))[as.factor(Rep)]

plot(PCA$x[,1],PCA$x[,2],col=colKP,main="Karypop",pch=16,cex=3)
#text(PCA$x[,1],PCA$x[,2],labels=interaction(Treat,Temp,Rep),col=colKP)

plot(PCA$x[,1],PCA$x[,2],col=colTemp,main="Temp",pch=16,cex=3)
#text(PCA$x[,1],PCA$x[,2],labels=interaction(Treat,Temp,Rep),col=colTemp)

plot(PCA$x[,1],PCA$x[,2],col=colPop,main="Replicate",pch=16,cex=3)
#text(PCA$x[,1],PCA$x[,2],labels=interaction(Treat,Temp,Rep),col=colPop)

plot(PCA$x[,3],PCA$x[,4],col=colKP,main="Karypop",pch=16,cex=3)
#text(PCA$x[,3],PCA$x[,4],labels=interaction(Treat,Temp,Rep),col=colKP)

plot(PCA$x[,3],PCA$x[,4],col=colTemp,main="Temp",pch=16,cex=3)
#text(PCA$x[,3],PCA$x[,4],labels=interaction(Treat,Temp,Rep),col=colTemp)

plot(PCA$x[,3],PCA$x[,4],col=colPop,main="Replicate",pch=16,cex=3)
#text(PCA$x[,3],PCA$x[,4],labels=interaction(Treat,Temp,Rep),col=colPop)

dev.off()

### Limma analysis with singe grouping factor

library(limma)

Treat=relevel(Treat,ref="FS")
dir.create(paste(O,"/limma",sep=""))
G=interaction(Treat,Temp)

design.SP=model.matrix(~0+G)
colnames(design.SP)=levels(G)

## Do VOOM
voomFullDGE.SP <- voom(FullDGE, design.SP, plot=F)

# write transformed data to file
out=paste(O,"/countdata/voomCounts.txt",sep="")
TAB=voomFullDGE.SP$E
write.table(data.frame("Gene"=rownames(TAB),TAB),file=out,row.names=F,quote = F)

corfit.SP <- duplicateCorrelation(voomFullDGE.SP,design.SP,block=Rep)
vfit.SP <- lmFit(voomFullDGE.SP, design.SP, block=Rep, correlation=corfit.SP$consensus)

############ Karypop

cont.matrix.SP.Treat <- makeContrasts(
  "FS-FI"=(FS.18+FS.25)/2-(FI.18+FI.25)/2,
  "FS-MS"=(FS.18+FS.25)/2-(MS.18+MS.25)/2,
  "MS-FI"=(MS.18+MS.25)/2-(FI.18+FI.25)/2,
  levels=design.SP)

# now calculate the contrasts:
vfit.Treat <- contrasts.fit(vfit.SP, cont.matrix.SP.Treat)
efit.Treat <- eBayes(vfit.Treat)

# now export the tables with the p-values and log2-fc for each of the contrasts
for (name in colnames(cont.matrix.SP.Treat)){
  # print the name
  print(name)
  
  # make an output file in outputfolder
  out=paste(O,"/limma/",name,".txt",sep="")
  
  # generate summary table with Benjamini Hochberg adjsuted p-values for all genes
  TAB=topTable(efit.Treat,coef=name, adjust.method="BH", n=Inf)
  
  # export table to text file
  write.table(data.frame("Gene"=rownames(TAB),TAB),file=out,row.names=F,quote = F)
}

# another cool thing about eBayes is that it calculates an F-value testing if the contrasts are different from 0. In the above case we are comparing all levels for the factor "Treat" only. Ebays combines the t-statistics for all the contrasts into an overall test of significance for that gene. The resulting F-value will thus tell us if "Treat" is significantly differnet from 0. See limma manual page 61. Note that I append another column to the dataframe with the P-value adjusted a la Benjamini-Hochberg
efit.F.Treat=data.frame("F"=efit.Treat$F,"Pval"=efit.Treat$F.p.value,"Pval.BH"=p.adjust(efit.Treat$F.p.value, method="BH"))
rownames(efit.F.Treat)=rownames(efit.Treat)

# export table as above for the contrasts

# make an output file in outputfolder
out=paste(O,"/limma/Treat-F-value.txt",sep="")

# export table to text file
write.table(data.frame("Gene"=rownames(efit.F.Treat),efit.F.Treat),file=out,row.names=F,quote = F)

############ Temp

cont.matrix.SP.Temp <- makeContrasts(
  "T18-25"=(FS.18+FI.18+MS.18)/3-(FS.25+FI.25+MS.25)/3,
  levels=design.SP)

# now calculate the contrasts:
vfit.Temp <- contrasts.fit(vfit.SP, cont.matrix.SP.Temp)
efit.Temp <- eBayes(vfit.Temp)

# now export the tables with the p-values and log2-fc for each of the contrasts
for (name in colnames(cont.matrix.SP.Temp)){
  # print the name
  print(name)
  
  # make an output file in outputfolder
  out=paste(O,"/limma/",name,".txt",sep="")
  
  # generate summary table with Benjamini Hochberg adjsuted p-values for all genes
  TAB=topTable(efit.Temp,coef=name, adjust.method="BH", n=Inf)
  
  # export table to text file
  write.table(data.frame("Gene"=rownames(TAB),TAB),file=out,row.names=F,quote = F)
}

efit.F.Temp=data.frame("F"=efit.Temp$F,"Pval"=efit.Temp$F.p.value,"Pval.BH"=p.adjust(efit.Temp$F.p.value, method="BH"))
rownames(efit.F.Temp)=rownames(efit.Temp)

# export table as above for the contrasts

# make an output file in outputfolder
out=paste(O,"/limma/Temp-F-value.txt",sep="")

# export table to text file
write.table(data.frame("Gene"=rownames(efit.F.Temp),efit.F.Temp),file=out,row.names=F,quote = F)

cont.matrix.SP.interactions <- makeContrasts(
  "FS-FI_int"=(FS.18-FS.25)-(FI.18-FI.25),
  "FI-MS_int"=(FI.18-FI.25)-(MS.18-MS.25),
  "FS-MS_int"=(FS.18-FS.25)-(MS.18-MS.25),
  levels=design.SP)

# now calculate the contrasts:
vfit.interactions <- contrasts.fit(vfit.SP, cont.matrix.SP.interactions)
efit.interactions <- eBayes(vfit.interactions)

# now export the tables with the p-values and log2-fc for each of the contrasts
for (name in colnames(cont.matrix.SP.interactions)){
  # print the name
  print(name)
  
  # make an output file in outputfolder
  out=paste(O,"/limma/",name,".txt",sep="")
  
  # generate summary table with Benjamini Hochberg adjsuted p-values for all genes
  TAB=topTable(efit.interactions,coef=name, adjust.method="BH", n=Inf)
  
  # export table to text file
  write.table(data.frame("Gene"=rownames(TAB),TAB),file=out,row.names=F,quote = F)
}

######## Interactions

efit.F.interactions=data.frame("F"=efit.interactions$F,"Pval"=efit.interactions$F.p.value,"Pval.BH"=p.adjust(efit.interactions$F.p.value, method="BH"))
rownames(efit.F.interactions)=rownames(efit.interactions)

# export table as above for the contrasts

# make an output file in outputfolder
out=paste(O,"/limma/interactions-F-value.txt",sep="")

# export table to text file
write.table(data.frame("Gene"=rownames(efit.F.interactions),efit.F.interactions),file=out,row.names=F,quote = F)

############ Other comparisons

cont.matrix.SP.other <- makeContrasts(
  "FS-FI_18"=FS.18-FI.18,
  "FS-MS_18"=FS.18-MS.18,
  "MS-FI_18"=MS.18-FI.18,
  "FS-FI_25"=FS.25-FI.25,
  "FS-MS_25"=FS.25-MS.25,
  "MS-FI_25"=MS.25-FI.25,
  levels=design.SP)

# now calculate the contrasts:
vfit.other <- contrasts.fit(vfit.SP, cont.matrix.SP.other)
efit.other <- eBayes(vfit.other)

# now export the tables with the p-values and log2-fc for each of the contrasts
for (name in colnames(cont.matrix.SP.other)){
  # print the name
  print(name)
  
  # make an output file in outputfolder
  out=paste(O,"/limma/",name,".txt",sep="")
  
  # generate summary table with Benjamini Hochberg adjsuted p-values for all genes
  TAB=topTable(efit.other,coef=name, adjust.method="BH", n=Inf)
  
  # export table to text file
  write.table(data.frame("Gene"=rownames(TAB),TAB),file=out,row.names=F,quote = F)
}

######## make heatmaps

dir.create(paste(O,"/heatmaps",sep=""))

library(gplots)

for (name in colnames(efit.Temp)) {
  col.pan <- colorpanel(100, "red", "yellow", "blue")
  To <- topTable(efit.Temp,coef=name, adjust.method="BH",n=Inf)
  logCPM <- cpmFullDGE[rownames(To[1:50,]),]
  logCPM <- t(scale(t(logCPM)))
  pdf(paste(O,"/heatmaps/Heatmap-",name,".pdf",sep=""), width=10, height=10)
  heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",trace="none", dendrogram="column", cexRow=1, cexCol=1.4, density.info="none",margin=c(10,9), lhei=c(2,10), lwid=c(2,6),main = paste(name,': Heatmap of top 50 genes',sep=""))
  dev.off()
}

for (name in colnames(efit.Treat)) {
  col.pan <- colorpanel(100, "red", "yellow", "blue")
  To <- topTable(efit.Treat,coef=name, adjust.method="BH",n=Inf)
  logCPM <- cpmFullDGE[rownames(To[1:50,]),]
  logCPM <- t(scale(t(logCPM)))
  pdf(paste(O,"/heatmaps/Heatmap-",name,".pdf",sep=""), width=10, height=10)
  heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",trace="none", dendrogram="column", cexRow=1, cexCol=1.4, density.info="none",margin=c(10,9), lhei=c(2,10), lwid=c(2,6),main = paste(name,': Heatmap of top 50 genes',sep=""))
  dev.off()
}

for (name in colnames(efit.interactions)) {
  col.pan <- colorpanel(100, "red", "yellow", "blue")
  To <- topTable(efit.interactions,coef=name, adjust.method="BH",n=Inf)
  logCPM <- cpmFullDGE[rownames(To[1:50,]),]
  logCPM <- t(scale(t(logCPM)))
  pdf(paste(O,"/heatmaps/Heatmap-",name,".pdf",sep=""), width=10, height=10)
  heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",trace="none", dendrogram="column", cexRow=1, cexCol=1.4, density.info="none",margin=c(10,9), lhei=c(2,10), lwid=c(2,6),main = paste(name,': Heatmap of top 50 genes',sep=""))
  dev.off()
}

for (name in colnames(efit.other)) {
  col.pan <- colorpanel(100, "red", "yellow", "blue")
  To <- topTable(efit.other,coef=name, adjust.method="BH",n=Inf)
  logCPM <- cpmFullDGE[rownames(To[1:50,]),]
  logCPM <- t(scale(t(logCPM)))
  pdf(paste(O,"/heatmaps/Heatmap-",name,".pdf",sep=""), width=10, height=10)
  heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none",trace="none", dendrogram="column", cexRow=1, cexCol=1.4, density.info="none",margin=c(10,9), lhei=c(2,10), lwid=c(2,6),main = paste(name,': Heatmap of top 50 genes',sep=""))
  dev.off()
}

  
## GO-term enrichment 
library(goseq)
dir.create(paste(O,"/GO",sep=""))
for (name in colnames(efit.Treat$lods)) {
  print(name)
  CS=topTable(efit.Treat,coef=name, adjust.method="BH", n=Inf)
  genes=as.integer(CS$adj.P.Val<0.05)
  names(genes)=row.names(CS)
  pwf=nullp(genes,"dm6","ensGene")
  GO.wall=goseq(pwf,"dm6","ensGene")
  GO.wall=data.frame(GO.wall,"P.adjus_BH"=p.adjust(GO.wall$over_represented_pvalue,method="BH"))
  name=gsub(":","-",name)
  out=paste(O,"/GO/GO-",name,".txt",sep="")
  write.table(GO.wall,file=out,row.names=F,quote = F)
}
for (name in colnames(efit.Treat$lods)) {
goana(efit.Treat, coef = name, FDR = 0.05, trend = FALSE,species="Dm",convert=F)
}

require("org.Dm.eg.db", character.only = TRUE)
GO2ALLEGS <- "org.Dm.egGO2ALLEGS"
EG.GO <- AnnotationDbi::toTable(get(GO2ALLEGS))
d <- duplicated(EG.GO[, c("gene_id", "go_id", "Ontology")])
EG.GO <- EG.GO[!d, ]
de.by.go <- split(EG.GO$gene_id, paste(EG.GO$go_id, EG.GO$Ontology, sep="."))

dir.create(paste(O,"/topgo",sep=""))
for (name in colnames(efit.Treat$lods)) {
  print(name)
  CS=topTable(efit.Treat,coef=name, adjust.method="BH", n=Inf)
  genes=as.integer(CS$adj.P.Val<0.05)
  E=mapIds(org.Dm.eg.db,
           keys=row.names(CS),
           column="ENTREZID",
           keytype="FLYBASE",
           multiVals="first")
  names(genes)=E
  cand=subset(genes,genes==1)
  GO=goana(names(cand),universe=names(genes),species="Dm",trend = "Length")
  
  GO=data.frame(GO,"P.adjus_BH"=p.adjust(GO$P.DE,method="BH"))
  
  DATAB <- lapply(de.by.go, FUN=function(x) { x[x %in% names(cand)] })
  for (i in 1:nrow(GO)){
    x=paste(rownames(GO)[i],GO$Ont[i],sep=".")
    if (length(DATAB[[x]])!=0){
      E=mapIds(org.Dm.eg.db,
               keys=DATAB[[x]],
               column="FLYBASE",
               keytype="ENTREZID",
               multiVals="first")
      GO$gene[i]=paste(E,collapse=",")
    }
    else {
      GO$gene[i]=""  
    }}
  
  name=gsub(":","-",name)
  out=paste(O,"/topgo/GO-",name,".txt",sep="")
  write.table(GO[ order(GO[,6]),],file=out,row.names=T,quote = F,sep="$")
}

for (name in colnames(efit.Temp$lods)) {
  print(name)
  CS=topTable(efit.Temp,coef=name, adjust.method="BH", n=Inf)
  genes=as.integer(CS$adj.P.Val<0.05)
  E=mapIds(org.Dm.eg.db,
           keys=row.names(CS),
           column="ENTREZID",
           keytype="FLYBASE",
           multiVals="first")
  names(genes)=E
  cand=subset(genes,genes==1)
  GO=goana(names(cand),universe=names(genes),species="Dm",trend = "Length")
  
  GO=data.frame(GO,"P.adjus_BH"=p.adjust(GO$P.DE,method="BH"))
  
  DATAB <- lapply(de.by.go, FUN=function(x) { x[x %in% names(cand)] })
  for (i in 1:nrow(GO)){
    x=paste(rownames(GO)[i],GO$Ont[i],sep=".")
    if (length(DATAB[[x]])!=0){
      E=mapIds(org.Dm.eg.db,
               keys=DATAB[[x]],
               column="FLYBASE",
               keytype="ENTREZID",
               multiVals="first")
      GO$gene[i]=paste(E,collapse=",")
    }
    else {
      GO$gene[i]=""  
    }}
  
  name=gsub(":","-",name)
  out=paste(O,"/topgo/GO-",name,".txt",sep="")
  write.table(GO[ order(GO[,6]),],file=out,row.names=T,quote = F,sep="$")
}


for (name in colnames(efit.interactions$lods)) {
  print(name)
  CS=topTable(efit.interactions,coef=name, adjust.method="BH", n=Inf)
  genes=as.integer(CS$adj.P.Val<0.05)
  E=mapIds(org.Dm.eg.db,
           keys=row.names(CS),
           column="ENTREZID",
           keytype="FLYBASE",
           multiVals="first")
  names(genes)=E
  cand=subset(genes,genes==1)
  GO=goana(names(cand),universe=names(genes),species="Dm",trend = "Length")
  
  GO=data.frame(GO,"P.adjus_BH"=p.adjust(GO$P.DE,method="BH"))
  
  DATAB <- lapply(de.by.go, FUN=function(x) { x[x %in% names(cand)] })
  for (i in 1:nrow(GO)){
    x=paste(rownames(GO)[i],GO$Ont[i],sep=".")
    if (length(DATAB[[x]])!=0){
      E=mapIds(org.Dm.eg.db,
               keys=DATAB[[x]],
               column="FLYBASE",
               keytype="ENTREZID",
               multiVals="first")
      GO$gene[i]=paste(E,collapse=",")
    }
    else {
      GO$gene[i]=""  
    }}
  
  name=gsub(":","-",name)
  out=paste(O,"/topgo/GO-",name,".txt",sep="")
  write.table(GO[ order(GO[,6]),],file=out,row.names=T,quote = F,sep="$")
}



for (name in colnames(efit.other$lods)) {
  print(name)
  CS=topTable(efit.other,coef=name, adjust.method="BH", n=Inf)
  genes=as.integer(CS$adj.P.Val<0.05)
  E=mapIds(org.Dm.eg.db,
           keys=row.names(CS),
           column="ENTREZID",
           keytype="FLYBASE",
           multiVals="first")
  names(genes)=E
  cand=subset(genes,genes==1)
  GO=goana(names(cand),universe=names(genes),species="Dm",trend = "Length")
  
  GO=data.frame(GO,"P.adjus_BH"=p.adjust(GO$P.DE,method="BH"))
  
  DATAB <- lapply(de.by.go, FUN=function(x) { x[x %in% names(cand)] })
  for (i in 1:nrow(GO)){
    x=paste(rownames(GO)[i],GO$Ont[i],sep=".")
    if (length(DATAB[[x]])!=0){
      E=mapIds(org.Dm.eg.db,
               keys=DATAB[[x]],
               column="FLYBASE",
               keytype="ENTREZID",
               multiVals="first")
      GO$gene[i]=paste(E,collapse=",")
    }
    else {
      GO$gene[i]=""  
    }}
  
  name=gsub(":","-",name)
  out=paste(O,"/topgo/GO-",name,".txt",sep="")
  write.table(GO[ order(GO[,6]),],file=out,row.names=T,quote = F,sep="$")
}

```

Then, we added additional information, such as gene name, length and genomic coordinates to the DE datasets and isolated candidate genes based on p-value

```bash

cd /data/RNASeq/limma
mkdir -p /data/RNASeq/limma/FullInfo/candidates

for i in *.txt
do
    ## add info about genes
    python3 /scripts/v3/AddGeneInfo.py \
        /data/RNASeq/limma/${i} \
        /data/RNASeq/dmel-all-r6.17.gtf.gz \
        > /data/RNASeq/limma/FullInfo/${i}

    ## isolate candidate genes with p<0.05
    printf 'Gene\tFBgn\tlogFC\tPval\n' > /data/RNASeq/limma/FullInfo/candidates/${i}
    awk '$12<0.05 {print $1"\t"$2"\t"$8"\t"$12}' ${i} >> /data/RNASeq/limma/FullInfo/candidates/${i}

done
```

#### 4.2) Test for positional overrepresentation

```bash
## test for overrepresentation of candidates with Fisher's Exact tests
mkdir /data/RNASeq/overrepresentation

printf "Dataset\tCand-Inv\tCand-NonInv\tNonCand-Inv\tNonCand-NonInv\tFET-pval\tCand-Inv_Breakpoint\tCand-Inv_Body\tNonCand-Inv_Breakpoint\tNonCand-Inv_Body\tFET-pval-InsideInv\n" \
> /data/RNASeq/overrepresentation/FET.txt

for i in /data/RNASeq/limma/FullInfo/*
do

    python3 /scripts/FET.py --input ${i} \
    >> /data/RNASeq/overrepresentation/FET.txt

done

## make manhattan plots
mkdir -p /data/RNASeq/overrepresentation/manhattan/logFC
mkdir -p /data/RNASeq/overrepresentation/manhattan/P

for i in /data/RNASeq/limma/FullInfo/*.txt
do

    ## Ignore datasets based on F-statistics generated above.
    if [[ ${i} != *"-F-"* ]]
    then

        ## make manhattan plots based on the -log10 P-values
        python3 /scripts/plot-manhattan.py \
            --full ${i} \
            --CF -2 \
            --log \
            --out /data/RNASeq/overrepresentation/manhattan/P

        ## make manhattan plots based on the log2-flod difference
        python3 /TransInv/scripts/plot-manhattan.py \
            --full ${i} \
            --CF 7 \
            --out /data/RNASeq/overrepresentation/manhattan/logFC

    fi

done
```

#### 4.3) Compare by temperarture

Next, we tested how the p-value distribution differed by temperature. We therefore manually compiled a table with all p-values from the different analyses (see data/)

```R

dir.create('/data/RNASeq/CompareByTemp')

## load pvalue data from limma
Data=read.table("/data/CompareByTemp_Pval.txt",header=T)

## Compare p-value distributions of 18 and 25 temp 

pdf("/data/RNASeq/CompareByTemp/Pval-genewise-byTemp.pdf",width=12,height=8)
par(mfrow=(c(3,3)))

plot(Data$FS.FI.25,Data$FS.FI.18,col=rgb(0,0,0,0.2),pch=16,main = "FS vs. FI - Temp",xlab="25C", ylab="18C")
plot(Data$FS.MS.25,Data$FS.MS.18,col=rgb(0,0,0,0.2),pch=16,main = "FS vs. MS - Temp",xlab="25C", ylab="18C")
plot(Data$MS.FI.25,Data$MS.FI.18,col=rgb(0,0,0,0.2),pch=16,main = "MS vs. FI - Temp",xlab="25C", ylab="18C")
plot(Data$FS.FI.25,Data$FS.FI,col=rgb(0,0,0,0.2),pch=16,main = "FS vs. FI - 25 vs. Average",xlab="25C", ylab="average")
plot(Data$FS.MS.25,Data$FS.MS,col=rgb(0,0,0,0.2),pch=16,main = "MS vs. FI - 25 vs. Average",xlab="25C", ylab="average")
plot(Data$MS.FI.25,Data$MS.FI,col=rgb(0,0,0,0.2),pch=16,main = "MS vs. FI - 25 vs. Average",xlab="25C", ylab="average")
plot(Data$FS.FI.18,Data$FS.FI,col=rgb(0,0,0,0.2),pch=16,main = "FS vs. FI - 18 vs. Average",xlab="18C", ylab="average")
plot(Data$FS.MS.18,Data$FS.MS,col=rgb(0,0,0,0.2),pch=16,main = "MS vs. FI - 18 vs. Average",xlab="18C", ylab="average")
plot(Data$MS.FI.18,Data$MS.FI,col=rgb(0,0,0,0.2),pch=16,main = "MS vs. FI - 18 vs. Average",xlab="18C", ylab="average")

dev.off()

pdf("/data/RNASeq/CompareByTemp/Pval-genewise-byTemp-QQ.pdf",width=12,height=8)
par(mfrow=(c(3,3)))

plot(sort(Data$FS.FI.25),sort(Data$FS.FI.18),col=rgb(0,0,0,0.2),pch=16,main = "FS vs. FI - Temp",xlab="25C", ylab="18C")
abline(0,1,col="red",lty=2,lwd=2)
plot(sort(Data$FS.MS.25),sort(Data$FS.MS.18),col=rgb(0,0,0,0.2),pch=16,main = "FS vs. MS - Temp",xlab="25C", ylab="18C")
abline(0,1,col="red",lty=2,lwd=2)
plot(sort(Data$MS.FI.25),sort(Data$MS.FI.18),col=rgb(0,0,0,0.2),pch=16,main = "MS vs. FI - Temp",xlab="25C", ylab="18C")
abline(0,1,col="red",lty=2,lwd=2)
plot(sort(Data$FS.FI.25),sort(Data$FS.FI),col=rgb(0,0,0,0.2),pch=16,main = "FS vs. FI - 25 vs. Average",xlab="25C", ylab="average")
abline(0,1,col="red",lty=2,lwd=2)
plot(sort(Data$FS.MS.25),sort(Data$FS.MS),col=rgb(0,0,0,0.2),pch=16,main = "FS vs. MS - 25 vs. Average",xlab="25C", ylab="average")
abline(0,1,col="red",lty=2,lwd=2)
plot(sort(Data$MS.FI.25),sort(Data$MS.FI),col=rgb(0,0,0,0.2),pch=16,main = "MS vs. FI - 25 vs. Average",xlab="25C", ylab="average")
abline(0,1,col="red",lty=2,lwd=2)
plot(sort(Data$FS.FI.18),sort(Data$FS.FI),col=rgb(0,0,0,0.2),pch=16,main = "FS vs. FI - 18 vs. Average",xlab="18C", ylab="average")
abline(0,1,col="red",lty=2,lwd=2)
plot(sort(Data$FS.MS.18),sort(Data$FS.MS),col=rgb(0,0,0,0.2),pch=16,main = "FS vs. MS - 18 vs. Average",xlab="18C", ylab="average")
abline(0,1,col="red",lty=2,lwd=2)
plot(sort(Data$MS.FI.18),sort(Data$MS.FI),col=rgb(0,0,0,0.2),pch=16,main = "MS vs. FI - 18 vs. Average",xlab="18C", ylab="average")
abline(0,1,col="red",lty=2,lwd=2)

dev.off()

```

#### 4.4) Overlap between Genomic and RNASeq candidates

Finally, we tested for overlap between genomic and RNASeq candidates based on rank-rank hypergeometric tests as implemented in the [RRHO2](https://github.com/RRHO2/RRHO2) package in _R_. First perpare the two datasets

```bash

mkdir -p /data/candidates/RRHO

## rank genes based on average FST in the dataset from Florida
python3 /scripts/v3/RankGenesByStat.py \
    --input /data/PopGen/FI-FS_ann.weir.fst \
    --perc 0.25 \
    --gtf /data/RNASeq/dmel-all-r6.17.gtf.gz \
    > /data/candidates/RRHO/Genes.fst


## loop through all RNASeq comparisons
for i in FS-FI FS-FI_18 FS-FI_25 FS-FI_int
do

    mkdir /data/candidates/RRHO/RRHO-${i}

    ## make input matrix for RRHO
    python3 /scripts/v3/Input4RRHO.py \
        --inputs /data/candidates/RRHO/Genes.fst,/data/RNASeq/limma/FullInfo/${i}.txt \
        --genes 1,1 \
        --stats -2,-2 \
        --type FST,P \
        > /data/candidates/RRHO/Genes-${i}.rrho

    ## execute R code for each dataset 
    printf """
    library(RRHO2)
    Rdata=read.table('/data/candidates/RRHO/Genes-${i}.rrho',header=T)

    DNA=data.frame('gene'=Rdata\$genes1,'val'=Rdata\$rank1)
    RNA=data.frame('gene'=Rdata\$genes2,'val'=Rdata\$rank2)
    Test=RRHO2(DNA,
        RNA,
        plots = T,
        labels=c('DNA','RNA'),
        outputdir = '/data/candidates/RRHO/RRHO-${i}',
        alternative = 'enrichment')

    """ > /data/candidates/RRHO/RRHO-${i}/test.r

    Rscript /data/candidates/RRHO/RRHO-${i}/test.r

done

```

## References

Kapopoulou, A., Kapun, M., Pieper, B., Pavlidis, P., Wilches, R., Duchen, P., et al. 2020. Demographic analyses of a new sample of haploid genomes from a Swedish population of _Drosophila melanogaster_. Sci Rep 10: 22415. Nature Publishing Group.

Kapun, M., Schalkwyk, H. van, McAllister, B., Flatt, T. & Schlötterer, C. 2014. Inference of chromosomal inversion dynamics from Pool-Seq data in natural and laboratory populations of _Drosophila melanogaster_. Molecular Ecology 23: 1813–1827.

Kapun, M., Barrón, M.G., Staubach, F., Obbard, D.J., Wiberg, R.A.W., Vieira, J., et al. 2020. Genomic Analysis of European _Drosophila melanogaster_ Populations Reveals Longitudinal Structure, Continent-Wide Selection, and Previously Unknown DNA Viruses. Mol Biol Evol 37: 2661–2678. Oxford Academic.

Kapun, M., Mitchell, E.D., Kawecki, T.J., Schmidt, P. & Flatt, T. 2023. An Ancestral Balanced Inversion Polymorphism Confers Global Adaptation. bioRxiv.

Kofler, R., Pandey, R.V. & Schlotterer, C. 2011. PoPoolation2: identifying differentiation between populations using sequencing of pooled DNA samples (Pool-Seq). Bioinformatics 27: 3435–3436.

Lack, J.B., Lange, J.D., Tang, A.D., B, C.-D., Russell & Pool, J.E. 2016. A Thousand Fly Genomes: An Expanded _Drosophila_ Genome Nexus. Mol Biol Evol 33: msw195-3313.

Rane, R.V., Rako, L., Kapun, M., Lee, S.F. & Hoffmann, A.A. 2015. Genomic evidence for role of inversion 3RP of _Drosophila melanogaster_ in facilitating climate change adaptation. Molecular Ecology 24: 2423–2432.
