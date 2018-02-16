#!/usr/bin/perl

# MANUAL FOR extract_sequences.pl

=pod

=head1 NAME

extract_sequences.pl -- Extract a set of sequences from a FASTA file

=head1 SYNOPSIS

 extract_sequences.pl --lookup=/Path/to/input.btab --fasta=/Path/to/input.fasta --out=/Path/to/output.txt
                     [--help] [--manual]

=head1 DESCRIPTION

 Extract the sequences specific in the "lookup" file from a FASTA file and write
 them to the "out" file.
 
=head1 OPTIONS

=over 3

=item B<-l, --lookup>=FILENAME

File of sequence IDs. One per line (Required).

=item B<-f, --fasta>=FILENAME

Input FASTA file (Required).

=item B<-o, --out>=FILENAME

Output file in FASTA format. (Required)

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

=cut


use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Pod::Usage;

#ARGUMENTS WITH NO DEFAULT
my($lookup,$fasta,$outfile,$help,$manual);

GetOptions (	
                                "l|lookup=s"	=>	\$lookup,
                                "f|fasta=s"     =>      \$fasta,
				"o|out=s"	=>	\$outfile,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --lookup not found.\n\n", -exitval => 2, -verbose => 1)  if (! $lookup );
pod2usage( -msg  => "\n\n ERROR!  Required argument --fasta not found.\n\n", -exitval => 2, -verbose => 1)  if (! $fasta );
pod2usage( -msg  => "\n\n ERROR!  Required argument --out not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile );

my %Lookup;
my ($seqs2pull,$seqs_pulled) = (0,0);
my $print_flag = 0;

open(IN,"<$lookup") || die "\n Cannot open the file: $lookup\n";
while(<IN>) {
    chomp;
    $Lookup{$_} = $_;
    $seqs2pull++;
}
close(IN);

open(OUT,">$outfile") || die "\n Cannot write to $outfile\n";
open(IN,"<$fasta") || die "\n Cannot open the file: $fasta\n";
while(<IN>) {
    chomp;
    if ($_ =~ m/^>/) {
	$print_flag = 0;
	my $h = $_;
	$h =~ s/^>//;
	$h =~ s/ .*//; $h =~ s/\t.*//;
	if (exists $Lookup{$h}) {
	    $seqs_pulled++;
	    print OUT $_ . "\n";
	    $print_flag = 1;
	}
    }
    elsif ($print_flag == 1) {
	print OUT $_ . "\n";
    }
}
close(IN);
close(OUT);

print "\n $seqs_pulled / $seqs2pull extracted\n\n";

exit 0;
