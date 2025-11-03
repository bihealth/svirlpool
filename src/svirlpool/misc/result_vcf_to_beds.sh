#!/bin/bash

input=$1
# create temporary directory
wdir=$(mktemp -d)

for metric in fn fp tp-call; do

    bcftools query -i 'INFO/SVTYPE=="DEL"' -f '%CHROM\t%POS\t%INFO/SVLEN\n' $input |
        awk -v OFS='\t' ' { print $1,$2,$2-$3,"del:"$3} ' |
        sed 's/-//g' >$wdir/$metric.proto.bed

    bcftools query -i 'INFO/SVTYPE=="INS"' -f '%CHROM\t%POS\t%INFO/SVLEN\n' $input |
        awk -v OFS='\t' ' { print $1,$2,$2+1,"ins:"$3} ' |
        sed 's/-//g' >>$wdir/$metric.proto.bed

    bcftools query -i 'INFO/SVTYPE=="DUP"' -f '%CHROM\t%POS\t%INFO/SVLEN\n' $input |
        awk -v OFS='\t' ' { print $1,$2,$2+$3,"dup:"$3} ' |
        sed 's/-//g' >>$wdir/$metric.proto.bed

    bcftools query -i 'INFO/SVTYPE=="INV"' -f '%CHROM\t%POS\t%INFO/SVLEN\n' $input |
        awk -v OFS='\t' ' { print $1,$2,$2+$3,"inv:"$3} ' |
        sed 's/-//g' >>$wdir/$metric.proto.bed

    bcftools query -i 'INFO/SVTYPE=="CON"' -f '%CHROM\t%POS\t%INFO/SVLEN\n' $input |
        awk -v OFS='\t' ' { print $1,$2,$2+$3,"con:"$3} ' |
        sed 's/-//g' >>$wdir/$metric.proto.bed

    # print to std out via sort
    bedtools sort -i $wdir/$metric.proto.bed
    rm $wdir/$metric.proto.bed
done
# remove temporary directory
rm -rf $wdir
