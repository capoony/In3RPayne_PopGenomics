# Bioinformatic pipeline for population genomic inference of _In(3R)Payne_

see Material and Methods in [Kapun _et al._ (2023)]() for more details

### 1) Map sequencing data from USA

For each of the newly sequenced library in the USA, test the quality with FASTQC, trim the raw reads with cutadapt, map the reads with bbmap, sort and deduplicate the BAM file with picard and realign around InDels with GATK

```bash
## loop through raw reads in folder /data and store BAM in /data/mapping
for i in <samplenames>
do

sh shell/mapping.sh \
/data/USA/${i}_R1.fq.gz	\
/data/USA/${i}_R2.fq.gz	\	
${i} \
/data/USA/mapping/${i} \
bbmap

done 
```

### 2) Generate phased data for samples from the USA

Following the approach in [Kapun _et al._ (2014)]() bioninformatically obtain haploid genomes from hemiclones

```bash
## compile list of all BAM files. IMPORTANTLY, make sure that NG9, the reference strain is in the first row. 
for i in /data/USA/mapping/*.bam
do
echo $i >> /data/USA/haplotypes/USA_bamlist
done 

## synchronize as MPILEUP based on hologenome from Kapun et al. (2020) only including 3L and 3R
samtools mpileup \
-B \
-f reference/Dmel_6.04_hologenome_v2.fasta \
-b /data/USA/haplotypes/USA_bamlist \
-l data /regions.bed.txt > /data/USA/haplotypes/USA.mpileup

## use PoPoolation2 (Kofler et al. 2011) to covert MPILEUP to SYNC format
java -jar /scripts/popoolation2_1201/mpileup2sync.jar \
--input /data/USA/haplotypes/USA.mpileup \
--threads 20 \
--output /data/USA/haplotypes/USA.sync

## calculate max coverage threshold (see Kapun et al. 2014)
python2.7 /scripts/max_coverage.py \
--input /data/USA/haplotypes/USA.sync \
--max-coverage 0.05 \
--output /data/USA/haplotypes/USA.cov

## use GNU parallel to parallelize the phasing (see Kapun et al. 2014 for details)
parallel -a /data/USA/haplotypes/USA.sync \
 -k \
 --pipepart \
 --cat python2.7 /scripts/extract_consensus.py \
 --input {} \
 --min-coverage 10 \
 --min-count 10 \
 --max-coverage /data/USA/haplotypes/USA.cov \
 --output /data/consensus/usa_min10_max005_mc10 \
 | gzip > ta/consensus/usa_min10_max005_mc10.consensus.gz
```

### 3) Obtain and map data from Zambia

```bash
## download raw data from SRA and map with same pipeline as above
while IFS=',' read -r name SRA
do 
    sh /shell/obtain-n-map-africa.sh \
    /data/Zambia/
    ${name} \
    ${SRA}
done < /data/DGN_SRA.txt
```







## References

Kapun, M., Schalkwyk, H. van, McAllister, B., Flatt, T. & Schlötterer, C. 2014. Inference of chromosomal inversion dynamics from Pool-Seq data in natural and laboratory populations of Drosophila melanogaster. Molecular Ecology 23: 1813–1827.

Kapun, M., Barrón, M.G., Staubach, F., Obbard, D.J., Wiberg, R.A.W., Vieira, J., et al. 2020. Genomic Analysis of European Drosophila melanogaster Populations Reveals Longitudinal Structure, Continent-Wide Selection, and Previously Unknown DNA Viruses. Mol Biol Evol 37: 2661–2678. Oxford Academic.

Kofler, R., Pandey, R.V. & Schlotterer, C. 2011. PoPoolation2: identifying differentiation between populations using sequencing of pooled DNA samples (Pool-Seq). Bioinformatics 27: 3435–3436.
