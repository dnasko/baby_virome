#!/usr/bin/perl

# MANUAL FOR taxonomic_assignment.pl

=pod

=head1 NAME

taxonomic_assignment_update.pl -- Calculate per-contig and overall taxonomic counts from BLAST output

=head1 SYNOPSIS

 taxonomic_assignment.pl --btab=/Path/to/input1.btab,/to/input2.btab --names /Path/to/names.dmp --nodes /Path/to/nodes.dmp --out /Path/to/output.txt [--abundance=/Path/to/abun.txt]
                     [--help] [--manual]

=head1 DESCRIPTION

 Calculates per-contig taxonomic assignment using the btab output of a BLASTp
 against the UniRef, SEED, or Phage SEED databases. Taxonomy is assigned 
 based on max(sum(bit score)) for each taxon that a contig's collection of
 ORFs have.

 If an abundnace file is passed (--abundance) then the abundnaces in the tpm column
 will be used to claculate abundance of a given taxon. Otherwise raw counts will be used.
 
=head1 OPTIONS

=over 3

=item B<-b, --btab>=FILENAME

BLAST tabular output(s) from a search against UniRef or SEED/Phage SEED. Multiple BTAB's can be passed using commas. (Required)

=item B<-a, --abundance>=FILENAME

The contig abundance file created from the bam2contig_abundance.pl script. Will use the `tpm` field (Optional).

=item B<-na, --names>=FILENAME

Path to the NCBI names.dmp file from the taxdump directory. (Required)

=item B<-no, --nodes>=FILENAME

Path to the NCBI nodes.dmp file from the taxdump directory. (Required) 

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

Copyright 2018 Daniel Nasko.  
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
my($btab,$names,$nodes,$abundance,$outfile,$help,$manual);

GetOptions (	
                                "b|btab=s"	=>	\$btab,
                                "na|names=s"    =>      \$names,
                                "no|nodes=s"    =>      \$nodes,
                                "a|abundance=s" =>      \$abundance,
				"o|out=s"	=>	\$outfile,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --btab not found.\n\n", -exitval => 2, -verbose => 1)  if (! $btab );
pod2usage( -msg  => "\n\n ERROR!  Required argument --names not found.\n\n", -exitval => 2, -verbose => 1)  if (! $names );
pod2usage( -msg  => "\n\n ERROR!  Required argument --nodes not found.\n\n", -exitval => 2, -verbose => 1)  if (! $nodes );
pod2usage( -msg  => "\n\n ERROR!  Required argument --out not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile);

## Build the taxonomt lookup
my %Names;
my %Tree;
my %Taxa;

open(IN,"<$names") || die "\n Cannot open the file: $names\n";
while(<IN>) {
    chomp;
    my @a = split(/\t/, $_);
    if ($a[6] eq "scientific name") {
	$Names{$a[0]} = $a[2];
    }
}
close(IN);

open(IN,"<$nodes") || die "\n Cannot open the file: $nodes\n";
while(<IN>) {
    chomp;
    my @a = split(/\t/, $_);
    my $parent = $a[2];
    my $child  = $a[0];
    $Tree{$child} = $parent;
}
close(IN);

foreach my $child (keys %Tree) {
    my $parent = $Tree{$child};
    my @Taxa = get_lineage($child);
    my $lineage = join("\t", @Taxa);
    $Taxa{$child} = $lineage;
}

## Load the abundance information
my %Abundance;
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
	    $Abundance{$a[0]} = $a[5]; ## we want the TPM column
	}
	$line_count++;
    }
    close(IN);
}

## Parse through the btab files
my @Btabs = split(/,/, $btab);
foreach my $btab_file (@Btabs) {
    if ($btab_file =~ m/\.gz$/) {
	open(IN,"gunzip -c $btab_file|") || die "\n Error: Cannot open the btab file: $btab_file\n";
    }
    else {
	open(IN,"<$btab_file") || die "\n Error: Cannot open the btab file: $btab_file\n";
    }
    while(<IN>) {
	chomp;
	my @a = split(/\t/, $_);
	my $taxid = "Unknown";
	if ($a[1] =~ m/^fig/) { ## If its a SEED subject sequence
	    $taxid = get_tax_seed($a[1]);
	}
	elsif ($a[1] =~ m/^UniRef100_/) { ## If its a UniRef subject sequence
	    $taxid = get_tax_uniref($a[-1]); ## Need the last field. It has the TaxID
	}
	else { ## Fail, because this is a sequence we didn't expect
	    die "\n Error: this line contains neither a SEED nor a UniRef sequence ID: $_\n\n";
	}
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

exit 0;

sub get_root
{
    my $s = $_[0];
    my @a = split(/_/, $s);
    pop(@a); pop(@a); pop(@a);
    return(join("_", @a));
}

sub get_tax_seed
{
    my $s = $_[0];
    $s =~ s/^fig\|//;
    $s =~ s/\..*//;
    return $s;
}

sub get_tax_uniref
{
    my $s = $_[0];
    $s =~ s/.*TaxID=//;
    $s =~ s/ .*//;
    return $s;
}

sub get_lineage
{
    my $s = $_[0];
    my $parent = $Tree{$s};
    my @a = ($s, $parent);
    while($parent != 1) {
	$parent = $Tree{$parent};
	push(@a, $parent);
    }
    @a = reverse(@a);
    my @b;
    foreach my $i (@a) {push(@b, $Names{$i}); }
    shift(@b);
    return(@b);
}
