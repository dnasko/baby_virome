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
2. Link PEGs to SubSystems and re-write the headers of the FASTA file
3. Make a BLASTp database from the FASTA file
4. Download the NCBI Taxonomy database
5. Parse these files to create the taxonomy strings for tax results

### Running the Baby Virome analysis

Blah.

*Rev DJN 07Nov2017*
