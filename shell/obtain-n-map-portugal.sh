#!/bin/sh

#  obtain-n-map-portugal.sh
#
#
#  Created by Martin Kapun on 31/01/16.
#
out=$1
name=$2
link1=$3
link2=$4

#download data

mkdir $out/data

curl -o $out/data/$name-1.fq.gz $link1

#
curl -o $out/data/$name-2.fq.gz $link2
#
#wait

## map data

mkdir -p $out/mapping/$name

sh /shell/mapping.sh \
    $out/data/$name-1.fq.gz \
    $out/data/$name-2.fq.gz \
    $name \
    $out/mapping/$name \
    bbmap
