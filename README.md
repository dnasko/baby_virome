# Baby Virome

##### A small virome analysis pipeline

The scripts in this repository serve to run the beta version of the Baby Virome Pipeline. The goal is to get taxonomic and functioanl information from shotgun viral metagenome (virome). This is done by performing a BLASTp of virome peptide open reading frames (ORFs) against the SEED database.

### Building the Baby Virome database

Building the Baby Virome database is easy, run the following command:

```bash
./scripts/create_database.sh
```

And a "BabyViromeDB" directory will be created with all of the files needed for the analysis. Specifically, it will:

1. Download the SEED database amino acid FASTA file and two files that link PEG's to SubSystems
2. Download the Phage SEED database amino acid FASTA file
3. Link PEGs to SubSystems and re-write the headers of the FASTA file
4. Make a BLASTp database from the FASTA file
5. Download the NCBI Taxonomy database
6. Parse these files to create the taxonomy strings for tax results

### Running the Baby Virome analysis

Once the database is downloaded you will want to predict ORFs on your contigs:

```bash
./software/mga_linux_ia64 -m input.fasta > output.mga
./scripts/mga2fasta.pl -f input.fasta -mg output.mga -p prefix -o ./
```

Next, BLAST the peptide ORFs against the SEED database. Feel free to use the para_blastp.pl script, which is usually faster than blastp:

```bash
./scripts/para_blastp.pl -q prefix.pep -d ./BabyViromeDB/SEED -o output.btab --outfmt='"6 std salltitles"' -e 1e-3 -t 30
```

Next we'll parse through the tabular output for results.

*Rev DJN 17Nov2017*
