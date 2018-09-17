#!/usr/bin/perl

# MANUAL FOR bam2contig_abundance.pl

=pod

=head1 NAME

bam2contig_abundance.pl -- create an abundance table from a sorted BAM file and FASTA file of contigs

=head1 SYNOPSIS

 bam2contig_abundance.pl --bam=/Path/to/input_sorted.bam --fasta=/Path/to/input_contigs.fasta --out=/Path/to/output.txt
                     [--help] [--manual]

=head1 DESCRIPTION

 Go through a sorted bam file and calculate the abundance for each sequence.
 Needs the contig FASTA file so we can calculate per-contig lengths.

 Five fields are printed out: contig_id, bases_recruited_to_contig, contig_len, coverage, normalized_coverage

=head1 OPTIONS

=over 3

=item B<-b, --bam>=FILENAME

Sorted BAM file (Required).

=item B<-f, --fasta>=FILENAME

Contig FASTA file used for the recruitment (Required).

=item B<-o, --out>=FILENAME

Output file in txt format. (Required)

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-m, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries.



=head1 AUTHOR

Written by Daniel Nasko, 
Center for Bioinformatics and Computational Biology, University of Maryland.

=head1 REPORTING BUGS

Report bugs to dan.nasko@gmail.com

=head1 COPYRIGHT

Copyright 2017 Daniel Nasko.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's 
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut


use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Pod::Usage;

#ARGUMENTS WITH NO DEFAULT
my($bam,$fasta,$outfile,$help,$manual);

GetOptions (	
                                "b|bam=s"	=>	\$bam,
                                "f|fasta=s"     =>      \$fasta,
				"o|out=s"	=>	\$outfile,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --bam not found.\n\n", -exitval => 2, -verbose => 1)  if (! $bam );
pod2usage( -msg  => "\n\n ERROR!  Required argument --fasta not found.\n\n", -exitval => 2, -verbose => 1)  if (! $fasta );
pod2usage( -msg  => "\n\n ERROR!  Required argument --out not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile);

## Make sure SAMtools is installed
my $samtools = `which samtools`;
unless ($samtools =~ m/samtools/ && $samtools !~ m/which/) { die "\n\n ERROR: External dependency samtools not installed in system PATH\n\n";}

my %Hash;
my %Size;
my $giga_bases=0;
my $sum_cov = 0; ## This holds the sum coverage, it's what's used in the TPM normalization
my ($h,$s) = ("","");

## Parse through the FASTA file to get the size for each sequence
open(IN,"<$fasta") || die "\n Cannot open the FASTA file: $fasta\n";
while(<IN>) {
    chomp;
    if ($_ =~ m/^>/) {
	unless ($h eq "") {
	    $Size{$h} = length($s);
	}
	$h = get_header($_);
	$s = "";
    }
    else { $s = $s . $_; }
}
close(IN);
$Size{$h} = length($s); ## One more time for the last sequence

## Time to go through the BAM file...
if (-e $bam) {
    open(my $cmd, '-|', 'samtools', 'depth', "$bam") or die $!;
    while (my $line = <$cmd>) {
	chomp($line);
	my @a = split(/\t/, $line);
	$Hash{$a[0]} += $a[2];
	$giga_bases += $a[2];
    }
    close $cmd;
}
else {
    die "\n Error: The BAM file is not there: $bam \n";
}

## Go through it once to get the sum_cov filled
foreach my $i (keys %Hash) {
    my $len;
    if (exists $Size{$i}) { $len = $Size{$i}; }
    else { die "\n Error: A sequence in the BAM file isnt in the contig FASTA file you provided: $i\n\n"; }
    my $cov = $Hash{$i}/$len;
    $sum_cov += $cov;
}

$sum_cov /= 1000000; ## Transcripts per *million* mapped. So divide by 1 million
$giga_bases /= 1000000000; ## This is for the RPKM-like normalization. Coverage per giga-base.

open(OUT,">$outfile") || die "\n Cannot open the file: $outfile\n";
print OUT join("\t", "#contig_id","bases_recruited","contig_length","coverage","cov_norm_rpkm","cov_norm_tpm") . "\n";
foreach my $i (keys %Hash) {
    my $len;
    if (exists $Size{$i}) { $len = $Size{$i}; }
    else { die "\n Error: A sequence in the BAM file isnt in the contig FASTA file you provided: $i\n\n"; }
    my $cov = $Hash{$i}/$len;
    my $norm_rpkm = $cov / $giga_bases;
    my $norm_tpm  = $cov / $sum_cov;
    print OUT join("\t", $i,$Hash{$i},$len,$cov,$norm_rpkm,$norm_tpm) . "\n";
}
close(OUT);

sub get_header
{
    my $s = $_[0];
    $s =~ s/^>//;
    $s =~ s/ .*//;
    return $s;
}

exit 0;
