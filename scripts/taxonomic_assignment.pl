#!/usr/bin/perl

# MANUAL FOR taxonomic_assignment.pl

=pod

=head1 NAME

taxonomic_assignment.pl -- Calculate per-ORF and overall taxonomic counts from BLAST output

=head1 SYNOPSIS

 taxonomic_assignment.pl --btab=/Path/to/input1.btab,/to/input2.btab --tax=/Path/to/taxa.lookup --out=/Path/to/output.txt [--abundance=/Path/to/abun.txt]
                     [--help] [--manual]

=head1 DESCRIPTION

 Calculates per-ORF taxonomic assignment using the btab output of a BLASTp
 against the Phage SEED database. Taxonomy is assigned based on max(sum(bit score))
 for each taxon that a contig's collection of ORF's have.

 If an abundnace file is passed (--abundnace) then the abundnaces in that 2-column file
 will be used to claculate abundance of a given taxon. Otherwise raw counts will be used.
 
=head1 OPTIONS

=over 3

=item B<-b, --btab>=FILENAME

BLAST tabular output(s) from a search against SEED/Phage SEED. Multiple BTAB's can be passed using commas. (Required)

=item B<-a, --abundance>=FILENAME

2-column file with <header> [TAB] <abundance>. Will be used to calculate abundance (Optional).

=item B<-t, --tax>=FILENAME

Taxonomy lookup file made by the build_taxonomy_table.pl script (Required).

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
my($btab,$taxonomy,$abundance,$outfile,$help,$manual);

GetOptions (	
                                "b|btab=s"	=>	\$btab,
                                "t|tax=s"       =>      \$taxonomy,
                                "a|abundance=s" =>      \$abundance,
				"o|out=s"	=>	\$outfile,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --btab not found.\n\n", -exitval => 2, -verbose => 1)  if (! $btab );
pod2usage( -msg  => "\n\n ERROR!  Required argument --taxonomy not found.\n\n", -exitval => 2, -verbose => 1)  if (! $taxonomy );
pod2usage( -msg  => "\n\n ERROR!  Required argument -outfile not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile);

my %Abundance;
my %Taxa;
my @Order;
my %Results;
my %ViromeResults;

my $out_per_query = $outfile . "_per_query.txt";
my $out_whole_set = $outfile . "_whole_virome.txt";
my $line_count=0;
if ($abundance) {
    open(IN,"<$abundance") || die "\n Cannot open the abundance file: $abundance\n";
    while(<IN>) {
	chomp;
	my @a = split(/\t/, $_);
	if ($line_count > 0) { ## if not on the first line
	    $Abundance{$a[0]} = $a[4];
	}
	$line_count++;
    }
    close(IN);
}

open(IN,"<$taxonomy") || die "\n Cannot open the taxonomy file: $taxonomy\n";
while(<IN>) {
    chomp;
    my @a = split(/\t/, $_);
    my $taxid = shift(@a);
    $Taxa{$taxid} = join("\t", @a);
}
close(IN);

my @Btabs = split(/,/, $btab);
foreach my $btab_file (@Btabs) {
    open(IN,"<$btab_file") || die "\n Error: Cannot open the btab file: $btab_file\n";
    while(<IN>) {
	chomp;
	my @a = split(/\t/, $_);
	my $taxid = get_tax($a[1]);
	my $root = get_root($a[0]);
	unless (exists $Results{$root}) {
	    push(@Order, $root);
	}
	$Results{$root}{$taxid} += $a[11];
    }
    close(IN);
}

open(OUT,">$out_per_query") || die "\n Cannot open the file: $out_per_query\n";
foreach my $i (@Order) {
    my $max=0;
    my $tax;
    foreach my $j (keys %{$Results{$i}}) {
	if ($Results{$i}{$j} > $max) {
	    $Results{$i}{$j} = $max;
	    $tax = $j;
	}
    }
    print OUT $i . "\t";
    if ($abundance) {
	if (exists $Abundance{$i}) {
	    print OUT $Abundance{$i} . "\t";
	    $ViromeResults{$tax} += $Abundance{$i};
	}
	else { print OUT "0\t"; }
    }
    else {
	print OUT "1\t";
	$ViromeResults{$tax}++;
    }
    print OUT $tax . "\t";
    if (exists $Taxa{$tax}) { print OUT $Taxa{$tax}; }
    print OUT "\n";
}
close(OUT);

open(OUT,">$out_whole_set") || die "\n Cannot open the file: $out_whole_set\n";
foreach my $i (sort { $ViromeResults{$b} <=> $ViromeResults{$a} } keys %ViromeResults) {
    print OUT $ViromeResults{$i} . "\t" . $i . "\t";
    if (exists $Taxa{$i}) {
	print OUT $Taxa{$i};
    }
    else { print OUT "Unknown" }
    print OUT "\n";
}
close(OUT);

sub get_root
{
    my $s = $_[0];
    my @a = split(/_/, $s);
    pop(@a); pop(@a); pop(@a);
    return(join("_", @a));
}

sub get_tax
{
    my $s = $_[0];
    $s =~ s/^fig\|//;
    $s =~ s/\..*//;
    return $s;
}

exit 0;
