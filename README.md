# Baby Virome

##### A small virome analysis pipeline

The scripts in this repository serve to run the beta version of the Baby Virome Pipeline. The goal is to get taxonomic and functional information from shotgun viral metagenome (virome). This is done by performing a BLASTp of virome peptide open reading frames (ORFs) against the [SEED](http://www.theseed.org/wiki/Main_Page) and [Phage SEED](http://www.phantome.org) databases.

### Building the Baby Virome database

Building the Baby Virome database is easy, run the following command:

```bash
./scripts/create_database.sh
```

And a "BabyViromeDB" directory will be created with all of the files needed for the analysis. Specifically, it will:

1. Download the SEED database amino acid FASTA file and two files that link PEG's to SubSystems
2. Download the Phage SEED database amino acid FASTA file
3. Link PEGs to SubSystems and re-write the headers of the SEED FASTA file
4. Make a BLASTp database from the FASTA files
5. Download the NCBI Taxonomy database
6. Parse these files to create the taxonomy strings for tax results

## Running the Baby Virome analysis

The input to the Baby Virome pipeline is assembeld contigs (presumably from a virome). Next, you need to predict open reading frames (ORFs) from your contigs:

```bash
./software/mga_linux_ia64 -m input.fasta > output.mga
./scripts/mga2fasta.pl -f input.fasta -mg output.mga -p prefix -o ./
```

### Abundance information

The Baby Virome pipeline relies on abundance information for each contig and ORF. I mean, without that information, what are you counting? Some *k*-mer based assemblers produce a coverage metric in the header (e.g. SPAdes). However, these abundance estimates aren't very good. It's better to recruit your QC'd reads back to the contigs you just assembled. You can do this using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (Langmead and Salzberg 2012) or [BBMap](https://sourceforge.net/p/bbmap/wiki/Home/). Here's an example using Bowtie2 and SAMtools:

```bash
bowtie2 -p 30 --very-sensitive-local -x BOWTIE_DATABASE \
	-1 ./01-flash/out.notCombined_1.fastq.gz \
	-2 ./01-flash/out.notCombined_2.fastq.gz \
	-U ./01-flash/out.extendedFrags.fastq.gz 2> bowtie2.log | samtools view -Sb -F 4 - | samtools sort -o out_sorted.bam -
```

Using either tool you'll need to produce a sorted BAM file (which the above example does) to get the abundance information we want. Using that sorted BAM file you can get abundance information for contigs and predicted ORFs. The output files will contain abundance (coverage) for each contig or orf and normalized abundance (coverage divided by giga bases recruited to whole assembly). Run like so:

```bash
./scripts/bam2abundance_spades.pl --bam=./out_sorted.bam --out=./contig_abundance.txt
./scripts/bam2orf_abundance.pl --bam=./out_sorted.bam --orfs=./orfs.pep --out=./orf_abundance.txt
```

Note that the contig abundance can only be calculated when the input contigs are from a SPAdes assembly. Happy to change this if someone would like me to...

### BLAST against SEED and Phage SEED

Next, BLASTp the peptide ORFs against the SEED and Phage SEED databases. Choose any expecancy value (E value) cutoff you'd like. It's required that you produce a special tabular file for the output though (`-outfmt "6 std salltitles"`). Feel free to use the para_blastp.pl script, which is usually faster than blastp:

```bash
./scripts/para_blastp.pl -q prefix.pep -d ./BabyViromeDB/SEED -o seed.btab --outfmt='"6 std salltitles"' -e 1e-3 -t 30
./scripts/para_blastp.pl -q prefix.pep -d ./BabyViromeDB/Phage_SEED -o phage_seed.btab --outfmt='"6 std salltitles"' -e 1e-3 -t 30
```

Next we'll parse through the tabular output for results.

### Taxonomic assignment

Each contig that has an ORF that produced a BLASTp hit will be assigned a taxonomy ID. Run the following command to get two tab-delimmited output files:

```bash
taxonomic_assignment.pl \
	--btab seed.btab,phage_seed.btab \
	--abundance ./bowtie/contig_abundance.txt \
	--tax ./BabyViromeDB/taxonomy_lookup.txt \
	--out taxonomy_results
```

Note that you need to pass both SEED and Phage SEED results to this script (separated by one comma and **no space**). This will produce the "whole virome" and "per query" output files that summarize the abundance of taxonomy throughout the metagenome and the taxonomic assignment for each query. I'd show you an example of what the output files should look like, but they contain many columns, typically one for abundance, one for the NCBI taxonomy numeric identifier, and one field for each level of taxonomy (domain, phylum, ..., genus, species).

### Functional assignment

Each peptide ORF with a hit will be assiend a SEED subsystem. Run the following command to get two tab-delimmited output files:

```bash
function_assignment.pl \
	--btab seed.btab \
	--abundance orf_abundance.txt \
	--out viral_functional
```

Two output files will be written, a "whole virome" and a "per query" file. The whole virome file contains two columns showing the overall abundance of various SEED subsystems throughout the entire metagenome, e.g.:

```bash
43526.42        DNA Metabolism
23079.33        Phages, Prophages, Transposable elements
18775.30        Carbohydrates
17644.79        Protein Metabolism
16199.73        Amino Acids and Derivatives
13130.15        Motility and Chemotaxis
12808.94        Cofactors, Vitamins, Prosthetic Groups, Pigments
11545.92        Cell Wall and Capsule
[...]           [...]
```

The per query file will contain three fields (sequnece ID, sequence abundance, functional assignment) with one row for each ORF assigned a function, e.g.:

```bash
NODE_115186_length_353_cov_2.60403_70_353_1     1.31        DNA Metabolism
NODE_115189_length_353_cov_2.5604_353_1_1       2.55        Protein Metabolism
NODE_115196_length_353_cov_2.41275_353_77_1     2.42        Protein Metabolism
NODE_115200_length_353_cov_2.36577_1_163_1      1.50        Sulfur Metabolism
NODE_115200_length_353_cov_2.36577_163_353_2    2.28        Carbohydrates
```

By default functions will be assigned using both bacterial and viral reference sequences from the SEED database. If you would like only viral functions this can be done using the --viruses_only flag and passing the taxonomy_lookup file:

```bash
	function_assignment.pl \
	--btab seed.btab \
    --abundance orf_abundance.txt \
    --out viral_functional
	--viruses_only ./BabyViromeDB/taxonomy_lookup.txt
```

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

![alt text](https://github.com/dnasko/baby_virome/blob/master/images/example_1.jpg)

And here's another where we're only getting functions from viral hits:

![alt text](https://github.com/dnasko/baby_virome/blob/master/images/example_2.jpg)

### References

BLASTp: Altschul et al. 1997

Bowtie2: Langmead and Salzberg 2012

SAMtools: Li et al. 2009

MetaGene: Noguchi et al. 2006

Kraken: Ondov, B. D., Bergman, N. H., & Phillippy, A. M. (2011). Interactive metagenomic visualization in a Web browser. BMC bioinformatics, 12(1), 385.

SEED: Overbeek et al. 2005

*Rev DJN 29Nov2017*
