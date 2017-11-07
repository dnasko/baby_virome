#!/usr/bin/perl

# MANUAL FOR function_assignment.pl

=pod

=head1 NAME

function_assignment.pl -- Calculate per-ORF and overall functional counts from BLAST output

=head1 SYNOPSIS

 function_assignment.pl --btab=/Path/to/input.btab --out=/Path/to/output.txt [--abundance=/Path/to/abun.txt]
                     [--help] [--manual]

=head1 DESCRIPTION

 Calculates per-ORF functional assignment using the btab output of a BLASTp
 against the Phage SEED database. Function is assigned based on max(sum(bit score))
 for each function that a query has.

 If an abundnace file is passed (--abundnace) then the abundnaces in that 2-column file
 will be used to claculate abundance of a given function. Otherwise raw counts will be used.
 
=head1 OPTIONS

=over 3

=item B<-b, --btab>=FILENAME

BLAST tabular output from a search against Phage SEED. (Required)

=item B<-a, --abundance>=FILENAME

2-column file with <header> [TAB] <abundance>. Will be used to calculate abundance (Optional).

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
my($btab,$abundance,$outfile,$help,$manual);

GetOptions (	
                                "b|btab=s"	=>	\$btab,
                                "a|abundance=s" =>      \$abundance,
				"o|out=s"	=>	\$outfile,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --btab not found.\n\n", -exitval => 2, -verbose => 1)  if (! $btab );
pod2usage( -msg  => "\n\n ERROR!  Required argument -outfile not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile);

my %Abundance;
if ($abundance) {
    open(IN,"<$abundance") || die "\n Cannot open the abundance file: $abundance\n";
    while(<IN>) {
	chomp;
	my @a = split(/\t/, $_);
	$Abundance{$a[0]} = $a[1];
    }
    close(IN);
}

open(IN,"<$btab") || die "\n Error: Cannot open the file: $btab\n";
while(<IN>) {
    chomp;
    my @a = split(/\t/, $_);
    
}
close(IN);

exit 0;
