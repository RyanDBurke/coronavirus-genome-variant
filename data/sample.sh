#!/bin/bash

## build the index for Bowtie2
bowtie2-build 2019-nCoV.fa 2019-nCoV_bt2_idx

## sample 1000 reads at (pseudo) random using seqtk (very specific path on my system below)
~/miniconda3/bin/seqtk sample -s 423 reads.fa 1000 > reads_s1000.fa

## align with bowtie2 (note, it uses an *affine* gap penalty, which you are not required to implement)
## so the alignments you get may not all be identical to those returned by Bowtie2, but your alignments should
## be optimal under your scoring function
bowtie2 --rdg 0,2 --rfg 0,2 --mp 2 -x 2019-nCoV_bt2_idx -f -U reads_s1000.fa -S mapped_s1000.sam

## convert to bam (if we want)
samtools view -b -o mapped_s1000.bam mapped_s1000.sam

## sort (if we want)
samtools sort mapped_s1000.bam > mapped_s1000_sorted.bam

## index (if we want)
samtools index mapped_s1000_sorted.bam
