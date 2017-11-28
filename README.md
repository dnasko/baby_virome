# Baby Virome

##### A small virome analysis pipeline

The scripts in this repository serve to run the beta version of the Baby Virome Pipeline. The goal is to get taxonomic and functional information from shotgun viral metagenome (virome). This is done by performing a BLASTp of virome peptide open reading frames (ORFs) against the SEED and Phage SEED databases.

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

## Running the Baby Virome analysis

The input to the Baby Virome pipeline is assembeld contigs (presumably from a virome). Next, you need to predict open reading frames (ORFs) from your contigs:

```bash
./software/mga_linux_ia64 -m input.fasta > output.mga
./scripts/mga2fasta.pl -f input.fasta -mg output.mga -p prefix -o ./
```

### Abundance information

The Baby Virome pipeline relies on abundnace information for each contig. I mean, without that information, what are you counting? Some *k*-mer based assemblers produce a coverage metric in the header (e.g. SPAdes). These abundance estimates aren't very good. It's better to recruit your QC'd reads back to the contigs you just assembled. You can do this using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) or [BBMap](https://sourceforge.net/p/bbmap/wiki/Home/). Here's an example using Bowtie2 and SAMtools:

```bash
bowtie2 -p 30 --very-sensitive-local -x MD10_VIRAL_OCT \
	-1 ./01-flash/out.notCombined_1.fastq.gz \
	-2 ./01-flash/out.notCombined_2.fastq.gz \
	-U ./01-flash/out.extendedFrags.fastq.gz 2> bowtie2.log | samtools view -Sb -F 4 - | samtools sort -o out_sorted.bam -
```

Using either tool you'll need to produce a sorted BAM file to get the abundance information we want. Using that sorted BAM file you can get abundance information for contigs and predicted ORFs, like so:

```bash
./scripts/bam2abundance_spades.pl --bam=./out_sorted.bam --out=./contig_abundance.txt
./scripts/bam2orf_abundance.pl --bam=./out_sorted.bam --orfs=./orfs.pep --out=./orf_abundance.txt
```

### BLAST against SEED and Phage SEED

Next, BLAST the peptide ORFs against the SEED and Phage SEED databases. Feel free to use the para_blastp.pl script, which is usually faster than blastp:

```bash
./scripts/para_blastp.pl -q prefix.pep -d ./BabyViromeDB/SEED -o seed.btab --outfmt='"6 std salltitles"' -e 1e-3 -t 30
./scripts/para_blastp.pl -q prefix.pep -d ./BabyViromeDB/Phage_SEED -o phage_seed.btab --outfmt='"6 std salltitles"' -e 1e-3 -t 30
```

Next we'll parse through the tabular output for results.

### Krona plots

Krona can be a great tool for exploring hierarchical data. Fortunatly your taxonomy and functional results are hierarchical. We can generate a functional Krona plot using the following:

```bash
./scripts/function_assignment_2_krona.pl --btab=seed.btab \
	--subsys2peg=./BabyViromeDB/subsystems2peg \
	--subsys=./BabyViromeDB/subsys.txt \
	--abundance=orf_abundance.txt \
	--out=dec_virus_functions.txt
```

Below is an example Krona plot using SEED subsystems:

![alt text](https://github.com/dnasko/baby_virome/tree/master/images/example_1.png "Example 1")

Be sure to cite Ondov, B. D., Bergman, N. H., & Phillippy, A. M. (2011). Interactive metagenomic visualization in a Web browser. BMC bioinformatics, 12(1), 385.

*Rev DJN 27Nov2017*
