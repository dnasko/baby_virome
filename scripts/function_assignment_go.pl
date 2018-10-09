#!/usr/bin/perl

# MANUAL FOR function_assignment_go.pl

=pod

=head1 NAME

function_assignment.pl -- Calculate per-ORF and overall functional counts from UniRef BLAST output

=head1 SYNOPSIS

 function_assignment.pl --btab /Path/to/input.btab --out /Path/to/output.txt [--cutoff 0.05] [--abundance /Path/to/abun.txt] [--viruses_only /Path/to/taxonomy_lookup.txt]
                     [--help] [--manual]

=head1 DESCRIPTION

 Calculates per-ORF functional assignment using the btab output of a BLASTp
 against the UniRef database. Function is reported as GO terms that are assigned
 to subjects that have hits to the query within --cut-off of the top-hit bit score.

 If an ORF abundance file is passed (--abundnace) then the abundances for each ORF will be calculated
 and reported for each function. Otherwise raw counts will be used.
 
 If you just want results for viruses, you pass the taxonomy_lookup file with the --viruses_only flag
=head1 OPTIONS

=over 3

=item B<-b, --btab>=FILENAME

BLAST tabular output from a search against Phage SEED. (Required)

=item B<-a, --abundance>=FILENAME

ORF abundance file (Optional).

=item B<-c, --cutoff>=FLOAT

Include GO terms within the cut-off of the top-hit bit score. (Default=0.05)

=item B<-o, --out>=FILENAME

Output file in txt format. (Required)

=item B<-v, --viruses_only>=FILENAME

Results for just viral hits. Point it to the taxonomy_lookup.txt file. (Optional)

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
my($btab,$abundance,$viruses_only,$outfile,$help,$manual);

## Arguments with defaults
my $cutoff = 0.05;

GetOptions (	
                                "b|btab=s"	=>	\$btab,
                                "a|abundance=s" =>      \$abundance,
                                "c|cutoff=s"    =>      \$cutoff,
                                "v|viruses_only=s" =>     \$viruses_only,
				"o|out=s"	=>	\$outfile,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --btab not found.\n\n", -exitval => 2, -verbose => 1)  if (! $btab );
pod2usage( -msg  => "\n\n ERROR!  Required argument --out not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile);
pod2usage( -msg  => "\n\n ERROR!  --cutoff must be between 0 and 0.5\n\n", -exitval => 2, -verbose => 1)  if ($cutoff < 0 || $cutoff > 0.05);
if ($abundance) {
    if ($abundance !~ m/orf/i) {
        print STDERR "\nWARNING: Make sure your abundance file is an ORF abundance file: $abundance\n\n";
    }
}

my %Abundance;
my @Order;
my %Results;
my %ViromeResults;
my %Viruses;
my %Tophit;


my $out_per_query = $outfile . "_per_query.txt";
my $out_whole_set = $outfile . "_whole_virome.txt";

if ($viruses_only) {
    open(my $vfh,'<', $viruses_only) || die "\n Cannot open the file: $viruses_only\n";
    while(<$vfh>) {
	chomp;
	my @a = split(/\t/, $_);
	if ($a[1] eq "Viruses") {
	    $Viruses{$a[0]} = 1;
	}
    }
    close($vfh);
}

my $bfh;
if ($btab =~ m/\.gz$/) {
    open($bfh,"gunzip -c $btab|") || die "\n Error: Cannot open the file: $btab\n";
}
else {
    open($bfh,'<', $btab) || die "\n Error: Cannot open the file: $btab\n";
}
while(<$bfh>) {
    chomp;
    my @a = split(/\t/, $_);
    my $subject = $a[1];
    my $bitscore = $a[11];
    my $cutoff_met = 0;
    if (exists $Tophit{$a[0]}) {
	if (($Tophit{$a[0]} * $cutoff) <= $bitscore) {
	    $cutoff_met = 1;
	}
    }
    else {
	$Tophit{$a[0]} = $bitscore;
	$cutoff_met = 1;
    }
    if ($cutoff_met == 1) {
	my @Fxn; ## needs to be an array, because sometimes we will have multiple functions
	if ($subject =~ m/^UniRef/) {
	    if ($a[-1] =~ m/GO=/) {
		@Fxn = get_fxn_uniref($a[-1]);
	    }
	}
	else { die "\n Error: This does not look like a UniRef file\n"; }
	if ($viruses_only) {
	    my $taxid;
	    if ($subject =~ m/^UniRef/) { $taxid = get_taxid_uniref($a[-1]); }
	    else { $taxid = get_taxid_seed($a[1]); }
	    if (exists $Viruses{$taxid}) {
		if ($a[-1] =~ m/GO=/) {
		    unless (exists $Results{$a[0]}) {
			push(@Order, $a[0]);
		    }
		    foreach my $f (@Fxn) {
			$Results{$a[0]}{$f} = 1;
		    }
		}
	    }
	}
	else {
	    if ($a[-1] =~ m/GO=/) {
		unless (exists $Results{$a[0]}) {
		    push(@Order, $a[0]);
		}
		foreach my $f (@Fxn) {
		    $Results{$a[0]}{$f} = 1;
		}
	    }
	}
    }
}
close($bfh);

my $line_count=0;
if ($abundance) {
    open(my $afh,'<', $abundance) || die "\n Cannot open the abundance file: $abundance\n";
    while(<$afh>) {
        chomp;
	if ($line_count > 0) {
	    my @a = split(/\t/, $_);
	    $Abundance{$a[0]} = $a[5]; ## We want the TPM abundance
	}
	$line_count++;
    }
    close($afh);
}

open(my $ofh,'>', $out_per_query) || die "\n Cannot open the file: $out_per_query\n";
foreach my $i (@Order) {
    my @Fxn;
    foreach my $j (sort keys %{$Results{$i}}) {
	push(@Fxn, $j);
	if (exists $Abundance{$i}) { $ViromeResults{$j} += $Abundance{$i}; }
	else { $ViromeResults{$j}++; }
    }
    print $ofh $i . "\t";
    if ($abundance) {
	if (exists $Abundance{$i}) {
	    print $ofh $Abundance{$i} . "\t";
	}
	else { print $ofh "0\t"; }
    }
    else {
	print $ofh "1\t";
    }
    print $ofh join(';',@Fxn) . "\n";
}
close($ofh);

open(my $owfh,'>', $out_whole_set) || die "\n Cannot open the file: $out_whole_set\n";
foreach my $i (sort { $ViromeResults{$b} <=> $ViromeResults{$a} } keys %ViromeResults) {
    print $owfh $ViromeResults{$i} . "\t" . $i . "\n";
}
close($owfh);

exit 0;

sub get_fxn_uniref
{
    my $s =$_[0];
    $s =~ s/.*GO=//;
    $s =~ s/ .*//;
    my @Ret = split(/;/, $s);
    return @Ret;
}

sub get_taxid_uniref
{
    my $s = $_[0];
    $s =~ s/.* TaxID=//;
    $s =~ s/ .*//;
    return $s;
}

sub get_taxid_seed
{
    my $s = $_[0];
    $s =~ s/fig\|//;
    $s =~ s/\..*//;
    return $s;
}

sub get_fxn_seed
{
    my $s = $_[0];
    my @Ret;
    $s =~ s/.*? //;
    push(@Ret, $s);
    return @Ret;
}
