#!/usr/bin/perl

# MANUAL FOR function_assignment_2_korona.pl

=pod

=head1 NAME

function_assignment_2_korona.pl -- Calculate overall functional counts from BLAST output and produce Krona input

=head1 SYNOPSIS

 function_assignment_2_korona.pl --subsys2peg=/Path/to/subsys2peg.txt --subsys=/Path/to/subsys.txt --btab=/Path/to/input.btab --out=/Path/to/output.txt [--abundance=/Path/to/abun.txt] [--viruses_only=/Path/to/taxonomy_lookup.txt]
                     [--help] [--manual]

=head1 DESCRIPTION

 Calculates functional assignment using the btab output of a BLASTp
 against the Phage SEED database. Function is assigned based on max(sum(bit score))
 for each function that a query has.

 If an ORF abundance file is passed (--abundnace) then the abundances for each ORF will be calculated
 and reported for each function. Otherwise raw counts will be used.

 If you just want results for viruses, you pass the taxonomy_lookup file with the --viruses_only flag
 
=head1 OPTIONS

=over 3

=item B<-b, --btab>=FILENAME

BLAST tabular output from a search against Phage SEED. (Required)

=item B<-sp, --subsys2peg>=FILENAME

The Subsystem to peg lookup. (Required).

=item B<-s, --subsys>=FILENAME

The subsys.txt file. (Required).

=item B<-a, --abundance>=FILENAME

ORF abundance file (Optional).

=item B<-v, --viruses_only>=FILENAME

Results for just viral hits. Point it to the taxonomy_lookup.txt file. (Optional)

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
my($btab,$subsys2peg,$subsys,$viruses_only,$abundance,$outfile,$help,$manual);

GetOptions (	
                                "b|btab=s"	=>	\$btab,
                                "sp|subsys2peg=s" =>      \$subsys2peg,
                                "s|subsys=s"      =>      \$subsys,
                                "v|viruses_only=s" =>     \$viruses_only,
                                "a|abundance=s" =>      \$abundance,
				"o|out=s"	=>	\$outfile,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --btab not found.\n\n", -exitval => 2, -verbose => 1)  if (! $btab );
pod2usage( -msg  => "\n\n ERROR!  Required argument --subsys2peg not found.\n\n", -exitval => 2, -verbose => 1)  if (! $subsys2peg );
pod2usage( -msg  => "\n\n ERROR!  Required argument --subsys not found.\n\n", -exitval => 2, -verbose => 1)  if (! $subsys );
pod2usage( -msg  => "\n\n ERROR!  Required argument -outfile not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile);
if ($abundance) {
    if ($abundance !~ m/orf/i) {
	print STDERR "\nWARNING: Make sure your abundance file is an ORF abundance file: $abundance\n\n";
    }
}

my %Peg2Subsys;
my %SubsysString; ## holds the 4 seed levels
my %Abundance; # Holding abundance information
my @Order; # Order of the query sequences
my %Results; 
my %ViromeResults;
my %Viruses;

if ($viruses_only) {
    open(IN,"<$viruses_only") || die "\n Cannot open the file: $viruses_only\n";
    while(<IN>) {
        chomp;
        my @a = split(/\t/, $_);
        if ($a[1] eq "Viruses") {
            $Viruses{$a[0]} = 1;
        }
    }
    close(IN);
}

open(IN,"<$subsys2peg") || die "\n Cannot open the file: $subsys2peg\n";
while(<IN>) {
    chomp;
    my @a = split(/\t/, $_);
    $Peg2Subsys{$a[2]} = $a[0] . "\t" . $a[1];
}
close(IN);

open(IN,"<$subsys") || die "\n Cannot open the file: $subsys\n";
while(<IN>) {
    chomp;
    my @a = split(/\t/, $_);
    $SubsysString{$a[2] . "\t" . $a[3]} = $_;
    pop(@a);
    $SubsysString{$a[2] . "\t"} = join("\t", @a);
}
close(IN);

open(IN,"<$btab") || die "\n Error: Cannot open the file: $btab\n";
while(<IN>) {
    chomp;
    my @a = split(/\t/, $_);
    my $fxn;
    if ($viruses_only) {
	my $taxid = get_taxid($a[1]);
	if (exists $Viruses{$taxid}) {
	    unless (exists $Results{$a[0]}) {
		push(@Order, $a[0]);
	    }
	    if (exists $Peg2Subsys{$a[1]}) {
		$fxn = $Peg2Subsys{$a[1]};
		$Results{$a[0]}{$fxn} += $a[11];
	    }
	    else { die "\n Cannot find subsys for the peg: $a[1]\n\n"; }
	}
    }
    else {
	unless (exists $Results{$a[0]}) {
	    push(@Order, $a[0]);
	}
	if (exists $Peg2Subsys{$a[1]}) {
	    $fxn = $Peg2Subsys{$a[1]};
	    $Results{$a[0]}{$fxn} += $a[11];
	}
	else { die "\n Cannot find subsys for the peg: $a[1]\n\n"; }
    }
}
close(IN);

my $line_count=0;
if ($abundance) {
    open(IN,"<$abundance") || die "\n Cannot open the abundance file: $abundance\n";
    while(<IN>) {
        chomp;
	if ($line_count > 0) {
	    my @a = split(/\t/, $_);
	    $Abundance{$a[0]} = $a[4];
	}
	$line_count++;
    }
    close(IN);
}

foreach my $i (@Order) {
    my $max=0;
    my $fxn;
    foreach my $j (keys %{$Results{$i}}) {
	if ($Results{$i}{$j} > $max) {
	    $Results{$i}{$j} = $max;
	    $fxn = $j;
	}
    }
    if ($abundance) {
	if (exists $Abundance{$i}) {
	    $ViromeResults{$fxn} += $Abundance{$i};
	}
	else {  $ViromeResults{$fxn} += 0; } # If it isn't there, then the abundance must be zero
    }
    else {
	$ViromeResults{$fxn}++;
    }
}

open(OUT,">$outfile") || die "\n Cannot open the file: $outfile\n";
foreach my $i (sort { $ViromeResults{$b} <=> $ViromeResults{$a} } keys %ViromeResults) {
    if (exists $SubsysString{$i}) {
	print OUT $ViromeResults{$i} . "\t" . $SubsysString{$i} . "\n";
    }
    else { die "\n Cannot find the string for: --> $i <--\n\n"; }
}
close(OUT);

sub get_taxid
{
    my $s = $_[0];
    $s =~ s/fig\|//;
    $s =~ s/\..*//;
    return $s;
}

sub get_fxn
{
    my $s = $_[0];
    $s =~ s/.*? //;
    return $s;
}

exit 0;
