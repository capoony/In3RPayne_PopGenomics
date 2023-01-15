#!/bin/sh

#  Master.sh
#
#
#  Created by Martin Kapun on 08/02/18.
#

###### PROCESS DATA

## at first map and call SNPs in pooled  North American and European Data jointly

sh /Volumes/MartinResearch1/Inv_sequencing/pooled_karyos/shell/map-n-call.sh

### get & map single individual seqencing data from USA

sh /Volumes/MartinResearch1/single_karyos/shells/download_florida_maine.sh

sh /Volumes/MartinResearch1/single_karyos/shells/map_USA.sh

## isolate haplotypes and combine datasets

sh /Volumes/MartinResearch1/single_karyos/shells/haplotypes.sh

## process African samples

sh /Volumes/MartinResearch1/Inv_sequencing/single_karyos/shells/mapAfrica.sh

sh /Volumes/MartinResearch1/Inv_sequencing/pooled_karyos/shell/Africa.sh

## process Australian samples

sh /Volumes/MartinResearch1/Inv_sequencing/single_karyos/shells/mapAustralia.sh

sh /Volumes/MartinResearch1/Inv_sequencing/pooled_karyos/shell/Australia.sh

## merge pseudo-pooled datasets

sh /Volumes/MartinResearch1/Inv_sequencing/shell/mergeData.sh

###### GENEALOGY

sh /Volumes/MartinResearch1/Inv_sequencing/single_karyos/shells/genealogy.sh

###### POPULATION GENETICS

# PI

sh /Volumes/MartinResearch1/Inv_sequencing/shell/pi.sh

# FST

sh /Volumes/MartinResearch1/Inv_sequencing/shell/FST.sh

# GLM

#sh /Volumes/MartinResearch1/Inv_sequencing/pooled_karyos/shell/GLM.sh

# candiates based on FST cutoffs, compare across continents and test for overlaps, compare Pool and Singles

sh /Volumes/MartinResearch1/Inv_sequencing/shell/FST-candidates.sh

# Linkage

sh /Volumes/MartinResearch1/Inv_sequencing/shell/LD.sh

##### TRANSCRIPTOMICS

sh /Volumes/MartinResearch1/Inv_sequencing/TransInv/shell/Master.sh

#### identify inv alleles of NA in other datasets

sh /Volumes/MartinResearch1/Inv_sequencing/shell/NA-alleles.sh

#### to test if the AF pattern in STD African strains is an artefact of the analyses, repeat the same with African candidates

sh /Volumes/MartinResearch1/Inv_sequencing/shell/Z-alleles.sh
## calculate dxy for all single hap data

sh /Volumes/MartinResearch1/Inv_sequencing/shells/dxy.sh

## repeat total PopGen analysis using PopGenome in R

sh /Volumes/MartinResearch1/Inv_sequencing/shells/PopGenome.sh

## make Table with final candidate genes from differnt sources

python3 /Volumes/MartinResearch1/Inv_sequencing/scripts/MakeGeneList.py \
--input /Volumes/MartinResearch1/Inv_sequencing/analyses/genelists2.txt \
--gtf /Volumes/MartinResearch1/Inv_sequencing/Transinv/RNASeq_new/reference/dmel-all-r6.17.gtf.gz \
--conv /Volumes/MartinResearch1/Inv_sequencing/analyses/candidates/Trans-Genomics/DROID/database/fbgn_annotation_ID_fb_2018_06.tsv.gz \
> /Volumes/MartinResearch1/Inv_sequencing/analyses/genelists_table2.txt
