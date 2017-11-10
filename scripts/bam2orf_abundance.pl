#!/usr/bin/perl

# MANUAL FOR bam2orf_abundance.pl

=pod

=head1 NAME

bam2orf_abundance.pl -- create an abundance table from a sorted BAM file with ORF header info

=head1 SYNOPSIS

 bam2orf_abundance.pl --bam=/Path/to/input_sorted.bam --orfs=/Path/to/orfs.fasta --out=/Path/to/output.txt
                     [--help] [--manual]

=head1 DESCRIPTION

 Go through a sorted bam file and calculate the abundance for each sequence.
 Need the ORFs FASTA file so the abundance of ORFs is calculated.

 Five fields are printed out: contig_id, bases_recruited_to_contig, contig_len, coverage, normalized_coverage 

=head1 OPTIONS

=over 3

=item B<-b, --bam>=FILENAME

Sorted BAM file (Required).

=item B<-or, --orfs>=FILENAME

Input ORF FASTA file (Required).

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
my($bam,$orfs,$outfile,$help,$manual);

GetOptions (	
                                "b|bam=s"	=>	\$bam,
                                "or|orfs=s"     =>      \$orfs,
				"o|out=s"	=>	\$outfile,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --bam not found.\n\n", -exitval => 2, -verbose => 1)  if (! $bam );
pod2usage( -msg  => "\n\n ERROR!  Required argument --orfs not found.\n\n", -exitval => 2, -verbose => 1)  if (! $orfs );
pod2usage( -msg  => "\n\n ERROR!  Required argument --out not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile);

my %Coord;
my %Abun;
my $giga_bases=0;

open(IN,"<$orfs") || die "\n Cannot open the file: $orfs\n";
while(<IN>) {
    chomp;
    if ($_ =~ m/^>/) {
	my $orfid = $_;
	$orfid =~ s/^>//;
	$orfid =~ s/ .*//;
	my ($root,$start,$stop) = parse_header($orfid);
	for (my $i=$start; $i<=$stop; $i++) {
	    $Coord{$root}{$i}=$orfid;
	}
    }
}
close(IN);


if (-e $bam) {
    open(my $cmd, '-|', 'samtools', 'depth', "$bam") or die $!;
    while (my $line = <$cmd>) {
	chomp($line);
	my @a = split(/\t/, $line);
	if (exists $Coord{$a[0]}{$a[1]}) {
	    $Abun{$Coord{$a[0]}{$a[1]}} += $a[2];
	    $giga_bases += $a[2];
	}
    }
    close $cmd;
}
else {
    die "\n Error: The BAM file is not there: $bam \n";
}
$giga_bases /= 1000000000;

open(OUT,">$outfile") || die "\n Cannot open the file: $outfile\n";
foreach my $i (keys %Abun) {
    my $orf_len = orf_len($i);
    my $cov = $Abun{$i}/$orf_len;
    my $norm = $cov/$giga_bases;

    print OUT $i . "\t" . $Abun{$i} . "\t" . $orf_len . "\t" . $cov . "\t" . $norm ."\n";
}
close(OUT);

sub orf_len
{
    my $s = $_[0];
    my @a = split(/_/, $s);
    pop(@a);
    my $stop = pop(@a);
    my $start = pop(@a);
    if ($start > $stop) { my $t = $start; $start = $stop; $stop = $t; }
    my $len = $stop - $start + 1;
    return $len;
}
sub parse_header
{
    my $s = $_[0];
    my @a = split(/_/, $s);
    pop(@a);
    my $stop = pop(@a);
    my $start = pop(@a);
    if ($start > $stop) { my $t = $start; $start = $stop; $stop = $t; }
    return (join("_", @a), $start, $stop);
}

sub get_len
{
    my $s = $_[0];
    $s =~ s/.*length_//;
    $s =~ s/_.*//;
    return $s;
}

exit 0;
