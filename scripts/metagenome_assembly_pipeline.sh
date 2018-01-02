#!/bin/bash
set -e

READ1=$1
READ2=$2
SAMPLE=$3
WORKING_DIR=$4
THREADS=$5

usage() {
    echo; echo "Usage: $0 /Path/to/read_pair_1.fastq /Path/to/read_pair_2.fastq Sample_Name /Path/to/working_directory/ #THREADS"
    echo;
    exit 1;
}

## Make sure you passed 5 arguments
if [ $# -gt 5 ] || [ $# -lt 5 ]
then
    usage
fi

## Create the output directories
mkdir -p $WORKING_DIR
mkdir -p $WORKING_DIR/00-trimmomatic
mkdir -p $WORKING_DIR/00-trimmomatic/$SAMPLE
mkdir -p $WORKING_DIR/01-flash
mkdir -p $WORKING_DIR/01-flash/$SAMPLE
mkdir -p $WORKING_DIR/02-spades
mkdir -p $WORKING_DIR/02-spades/$SAMPLE
mkdir -p $WORKING_DIR/03-bowtie2
mkdir -p $WORKING_DIR/03-bowtie2/$SAMPLE
mkdir -p $WORKING_DIR/04-metagene
mkdir -p $WORKING_DIR/04-metagene/$SAMPLE

## Trimmomatic
echo -n " Running trimmomatic ............... "
java -jar ~/bin/trimmomatic-0.36.jar PE -phred33 $READ1 $READ2 \
    $WORKING_DIR/00-trimmomatic/$SAMPLE/output_forward_paired.fq \
    $WORKING_DIR/00-trimmomatic/$SAMPLE/output_forward_unpaired.fq \
    $WORKING_DIR/00-trimmomatic/$SAMPLE/output_reverse_paired.fq \
    $WORKING_DIR/00-trimmomatic/$SAMPLE/output_reverse_unpaired.fq \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:60 \
    > $WORKING_DIR/00-trimmomatic/$SAMPLE/trim.log 2>&1
echo -n " Done "; date;

## FLASh
echo -n " Running FLASh ..................... "
flash -t ${THREADS} \
    --max-overlap 100 \
    --output-directory=$WORKING_DIR/01-flash/$SAMPLE \
    $WORKING_DIR/00-trimmomatic/$SAMPLE/output_forward_paired.fq \
    $WORKING_DIR/00-trimmomatic/$SAMPLE/output_reverse_paired.fq \
    > $WORKING_DIR/01-flash/$SAMPLE/flash.log 2>&1

cat $WORKING_DIR/00-trimmomatic/$SAMPLE/output_forward_unpaired.fq \
    $WORKING_DIR/00-trimmomatic/$SAMPLE/output_reverse_unpaired.fq \
    >> $WORKING_DIR/01-flash/$SAMPLE/out.extendedFrags.fastq
rm $WORKING_DIR/00-trimmomatic/$SAMPLE/*.fq ## clean up the trimmomatic results
echo -n " Done "; date;

## SPAdes
echo -n " Running MetaSPAdes ............... "
metaspades.py -o $WORKING_DIR/02-spades/$SAMPLE/ \
    --only-assembler \
    -t $THREADS \
    -1 $WORKING_DIR/01-flash/$SAMPLE/out.notCombined_1.fastq \
    -2 $WORKING_DIR/01-flash/$SAMPLE/out.notCombined_2.fastq \
    -s $WORKING_DIR/01-flash/$SAMPLE/out.extendedFrags.fastq \
    > $WORKING_DIR/01-flash/$SAMPLE/spades_run.log 2>&1
echo -n " Done "; date;

## Bowtie2
echo -n " Running recruitment ............... "
ln -s ../../02-spades/$SAMPLE/contigs.fasta $WORKING_DIR/03-bowtie2/$SAMPLE/contigs.fasta
bowtie2-build --threads $THREADS \
    -f $WORKING_DIR/03-bowtie2/$SAMPLE/contigs.fasta \
    $WORKING_DIR/03-bowtie2/$SAMPLE/METAGENOME \
    > $WORKING_DIR/03-bowtie2/$SAMPLE/bowtie2_build.log 2>&1

bowtie2 --very-sensitive-local -p $THREADS \
    -x $WORKING_DIR/03-bowtie2/$SAMPLE/METAGENOME \
    -1 $WORKING_DIR/01-flash/$SAMPLE/out.notCombined_1.fastq \
    -2 $WORKING_DIR/01-flash/$SAMPLE/out.notCombined_2.fastq \
    -U $WORKING_DIR/01-flash/$SAMPLE/out.extendedFrags.fastq 2> $WORKING_DIR/03-bowtie2/$SAMPLE/bowtie2.log | samtools view -Sb -F 4 - 2> $WORKING_DIR/03-bowtie2/$SAMPLE/samtools_view.log | samtools sort - -o $WORKING_DIR/03-bowtie2/$SAMPLE/out_sorted.bam 2> $WORKING_DIR/03-bowtie2/$SAMPLE/samtools_sort.log
echo -n " Done "; date;

## MetaGene
echo -n " Running MetaGene .................. "
mga_linux_ia64 -m $WORKING_DIR/02-spades/$SAMPLE/contigs.fasta > $WORKING_DIR/04-metagene/$SAMPLE/contigs.mga

mga2fasta.pl -f $WORKING_DIR/02-spades/$SAMPLE/contigs.fasta \
    -mg $WORKING_DIR/04-metagene/$SAMPLE/contigs.mga \
    -p $SAMPLE \
    -o $WORKING_DIR/04-metagene/$SAMPLE/
echo -n " Done "; date;

exit 0
