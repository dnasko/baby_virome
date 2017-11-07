#!/usr/bin/perl

# MANUAL FOR rebuild_seed_headers.pl

=pod

=head1 NAME

rebuild_seed_headers.pl -- reformat the SEED DB headers

=head1 SYNOPSIS

 rebuild_seed_headers.pl --subsys2peg=/Path/to/subsystems2peg --subsys=/Path/to/subsys --fasta=/Path/to/phage_seed.fasta
                     [--help] [--manual]

=head1 DESCRIPTION

 Update the header in the phage seed FASTA so blast results have everything we need.
 
=head1 OPTIONS

=over 3

=item B<-s2, --subsys2peg>=FILENAME

http://ftp.theseed.org/subsystems/subsystems2peg.gz. (Required) 

=item B<-s, --subsys>=FILENAME

http://ftp.theseed.org/subsystems/subsys.txt. (Required)

=item B<-f, --fasta>=FILENAME

SEED FAA FASTA file. (Required) 

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
my($subsys2peg,$subsys,$fasta,$help,$manual);

GetOptions (	
                                "s2|subsys2peg=s" =>	\$subsys2peg,
                                "s|subsys=s"    =>      \$subsys,
				"f|fasta=s"	=>	\$fasta,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --subsys2peg not found.\n\n", -exitval => 2, -verbose => 1) if (! $subsys2peg );
pod2usage( -msg  => "\n\n ERROR!  Required argument --subsys not found.\n\n", -exitval => 2, -verbose => 1)     if (! $subsys );
pod2usage( -msg  => "\n\n ERROR!  Required argument --fasta not found.\n\n", -exitval => 2, -verbose => 1)      if (! $fasta );

my %Peg;
my %Subsys;

open(IN,"<$subsys2peg") || die "\n Cannot open the file: $subsys2peg\n";
while(<IN>) {
    chomp;
    my @a = split(/\t/, $_);
    $Peg{$a[2]} = $a[0];
}
close(IN);

open(IN,"<$subsys") || die "\n Cannot open the file: $subsys\n";
while(<IN>) {
    chomp;
    my @a = split(/\t/, $_);
    $Subsys{$a[2]} = $a[0];
}
close(IN);

my $out = $fasta . "_updated.fasta";
my $print=0;
open(OUT,">$out") || die "\n Cannot write to $out\n";
open(IN,"<$fasta") || die "\n Cannot open the file: $fasta\n";
while(<IN>) {
    chomp;
    if ($_ =~ m/^>/) {
	$print = 0;
	my $fig = parse_fig($_);
	if (exists $Peg{$fig}) {
	    if (exists $Subsys{$Peg{$fig}}) {
	    	print OUT ">" . $fig . " " . $Subsys{$Peg{$fig}} . "\n";
	    	$print = 1;
	    }
	    # else { print OUT ">" .$fig . " [$Peg{$fig}] " . $all_else . "\n"; }
	    # print OUT ">" . $fig . " " . $Peg{$fig} . "\n";
	    # $print = 1;
	}
	# else { print OUT ">" .$fig . " [miss2] " . $all_else . "\n"; }
    }
    elsif ($print == 1) {
	print OUT $_ . "\n";
    }
}
close(IN);
close(OUT);

sub parse_fig
{
    my $s = $_[0];
    $s =~ s/^>//;
    $s =~ s/ .*//;
    return $s;
}

exit 0;
