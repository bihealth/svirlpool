#!/bin/bash

ref=$1
bm=$2
vcf=$3
outd=$4
bed=$5

export TMPDIR=$HOME/scratch/tmp
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

truvari bench \
    --passonly \
    --giabreport \
    -f $ref \
    -b $bm \
    -c $vcf \
    -o $TMPDIR/i_exist \
    --includebed $bed

cp $TMPDIR/i_exist/* $outd/.
