#!/bin/bash
set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

mkdir BabyViromeDB
cd BabyViromeDB

## Download SEED subsystem information
wget --ftp-user=anonymous --ftp-password=dnasko@umiacs.umd.edu "ftp://ftp.theseed.org/subsystems/subsystems2peg.gz"
wget --ftp-user=anonymous --ftp-password=dnasko@umiacs.umd.edu "http://ftp.theseed.org/subsystems/subsys.txt"
wget --ftp-user=anonymous --ftp-password=dnasko@umiacs.umd.edu "ftp://ftp.theseed.org/genomes/SEED/all.faa.gz"
gunzip subsystems2peg.gz
gunzip all.faa.gz
mv all.faa.gz SEED_original
${SCRIPT_DIR}/rebuild_seed_headers.pl -s2 subsystems2peg -s subsys.txt -f SEED_original
mv SEED_original_updated.fasta SEED
rm SEED_original

## Make the blast database
makeblastdb -in SEED -dbtype prot

## Download the taxonomy information
mkdir tax
cd tax
wget "http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
tar -xvf taxdump.tar.gz
rm taxdump.tar.gz

"${SCRIPT_DIR}/extract_viral_taxonomy.pl" -no nodes.dmp -na names.dmp -o ../taxonomy_lookup.txt
cd ../
rm -rf tax
cd ../
