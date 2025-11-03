for ubam in $(ls *.bam); do
    bn=$(basename $ubam .css.bam)
    bnn=${bn%.ccs.bam}
    echo $bnn
    samtools fastq -n -@ 30 $ubam >$bnn.fq
done

cat *.fq >all_reads.fastq
