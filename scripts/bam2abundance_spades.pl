#!/usr/bin/perl

# MANUAL FOR bam2abundance_spades.pl

=pod

=head1 NAME

bam2abundance_spades.pl -- create an abundance table from a sorted BAM file with SPAdes contigs

=head1 SYNOPSIS

 bam2abundance_spades.pl --bam=/Path/to/input_sorted.bam --out=/Path/to/output.txt
                     [--help] [--manual]

=head1 DESCRIPTION

 Go through a sorted bam file and calculate the abundance for each sequence.
 Can only be reliable with SPAdes contigs because we need the contig's length.

 Five fields are printed out: contig_id, bases_recruited_to_contig, contig_len, coverage, normalized_coverage

=head1 OPTIONS

=over 3

=item B<-b, --bam>=FILENAME

Sorted BAM file (Required).

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

Report bugs to dnasko@umiacs.umd.edu

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
my($bam,$outfile,$help,$manual);

GetOptions (	
                                "b|bam=s"	=>	\$bam,
				"o|out=s"	=>	\$outfile,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --bam not found.\n\n", -exitval => 2, -verbose => 1)  if (! $bam );
pod2usage( -msg  => "\n\n ERROR!  Required argument --out not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile);

my %Hash;
my $giga_bases=0;

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

$giga_bases /= 1000000000;

open(OUT,">$outfile") || die "\n Cannot open the file: $outfile\n";
print OUT "#contig_id\tbases_recruited\tcontig_length\tabundance\tnormalized_abundance\n";
foreach my $i (keys %Hash) {
    my $len = get_len($i);
    my $cov = $Hash{$i}/$len;
    my $norm = $cov / $giga_bases;
    print OUT $i . "\t" . $Hash{$i} . "\t" . $len . "\t" . $cov . "\t" . $norm . "\n";
}
close(OUT);

sub get_len
{
    my $s = $_[0];
    $s =~ s/.*length_//;
    $s =~ s/_.*//;
    return $s;
}

exit 0;
