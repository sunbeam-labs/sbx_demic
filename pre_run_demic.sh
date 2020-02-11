#!/usr/bin/env bash

#$ -cwd
#$ -r n
#$ -l h_vmem=4G
#$ -j y
#$ -pe smp 12

source ~/.bashrc
conda activate base
source config.sh

set -x 
set -e

find $SBM_OUT/qc/decontam -iname "*.fastq.gz" > ${PRJ_DIR}/decontam_list

#mkdir -p $SORTED_DIR
mkdir -p $BINNED_DIR
#
#for file in $(cat ${DATA_DIR}/mapped_list); do
#
#    SAMPLE=$(basename $file _R1.mapped.fasta)
#
#    if [[ ! -e ${OUTPUT_DIR}/"$SAMPLE"_sorted.bam ]]; then
#        samtools view -b -@ 12 ${OUTPUT_DIR}/$SAMPLE.sam | samtools sort -@ 12 - -o ${OUTPUT_DIR}/"$SAMPLE"_sorted.bam
#    fi
#
#    if [[ ! -e ${SORTED_DIR}/"$SAMPLE"_sorted.sam ]]; then
#        samtools view -h ${OUTPUT_DIR}/"$SAMPLE"_sorted.bam > ${SORTED_DIR}/"$SAMPLE"_sorted.sam 
#    fi
#
#done
#
run_MaxBin.pl -thread 10 -contig $CONTIGS -out $BINNED_DIR -reads_list ${PRJ_DIR}/decontam_list




