#!/bin/sh

#  obtain-n-map-africa.sh
#
#  Created by Martin Kapun on 03/02/16.
#

SRA=$3
name=$2
out=$1

mkdir $out/data

## Use Aspera Connect to speed up things
/Aspera\ Connect.app/Contents/Resources/ascp \
    -i /Aspera\ Connect.app/Contents/Resources/asperaweb_id_dsa.openssh \
    -T anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${SRA:0:6}/$SRA/$SRA.sra \
    $out/data

## Use SRAtoolkit 2.5.2 to convert SRA tp FASTQ
/sratoolkit.2.5.2-mac64/bin/fastq-dump.2.5.2 \
    --gzip \
    --split-3 \
    --outdir $out/data \
    -A $name \
    $out/data/$SRA.sra

rm $out/data/$SRA.sra

## map data
mkdir -p $out/mapping/$name

sh /shell/mapping.sh \
    $out/data/$name\_1.fastq.gz \
    $out/data/$name\_2.fastq.gz \
    $name \
    $out/mapping/$name bbmap
