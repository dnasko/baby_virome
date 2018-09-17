#!/usr/bin/perl

# MANUAL FOR taxonomic_summary.pl

=pod

=head1 NAME

 taxonomic_summary.pl -- Summarize the whole-library taxonomic abundance at a taxonomic level

=head1 SYNOPSIS

 taxonomic_summary.pl --taxass=/Path/to/whole_library_taxonomy.txt --names=/Path/to/names.dmp --nodes=/Path/to/nodes.dmp --level=phylum --out=/Path/to/output.txt
                     [--help] [--manual]

=head1 DESCRIPTION

 Calculate the abundance of taxons at a taxonomic level.
 
=head1 OPTIONS

=over 3

=item B<-ta, --taxass>=FILENAME

The whole library taxonomic assignment file from the `taxonomic_assignment.pl` script (Required).

=item B<-na, --names>=FILENAME

The names.dmp file from NCBI taxonomy database (Required).

=item B<-no, --nodes>=FILENAME

The nodes.dmp file from NCBI taxonomy database (Required).

=item B<-l, --level>=TYPE

The taxonomy leve you want the summary to be at, must be: phylum, family, genus, or species (Required)

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
my($taxassfile,$namesfile,$nodesfile,$level,$outfile,$help,$manual);

GetOptions (	
                                "ta|taxass=s"	=>	\$taxassfile,
                                "na|names=s"    =>      \$namesfile,
                                "no|nodes=s"    =>      \$nodesfile,
                                "l|level=s"     =>      \$level,
				"o|out=s"	=>	\$outfile,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --taxass not found.\n\n", -exitval => 2, -verbose => 1)  if (! $taxassfile );
pod2usage( -msg  => "\n\n ERROR!  Required argument --names not found.\n\n", -exitval => 2, -verbose => 1)  if (! $namesfile );
pod2usage( -msg  => "\n\n ERROR!  Required argument --nodes not found.\n\n", -exitval => 2, -verbose => 1)  if (! $nodesfile );
pod2usage( -msg  => "\n\n ERROR!  Required argument --level not found.\n\n", -exitval => 2, -verbose => 1)  if (! $level );
pod2usage( -msg  => "\n\n ERROR!  Required argument --out not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile);
my %ValidLevels = ('phylum'=>0,'class'=>0,'family'=>0,'genus'=>0,'species'=>0);
unless (exists $ValidLevels{$level}) { die "\n Error: The level you provided is not valid!\n\n"; }
my %Names;  ## The scientific name of a taxid
my %Level;  ## Tells you the level a taxid is at
my %Tree;   ## The taxonomic tree
my %Domain; ## The domain of a given taxid
my %Lookup; ## The --level of a given taxid
my %Results;

open(IN,"<$namesfile") || die "\n Error: Cannot open the file: $namesfile\n";
while(<IN>) {
    chomp;
    my @a = split(/\t/, $_);
    my $id = $a[0];
    my $name = $a[2];
    my $type = $a[6];
    if ($type eq "scientific name") {
	$Names{$id} = $name;
    }
}
close(IN);

open(IN,"<$nodesfile") || die "\n Error: Cannot open the file: $nodesfile\n";
while(<IN>) {
    chomp;
    my @a = split(/\t/, $_);
    my $child = $a[0];
    my $parent = $a[2];
    my $level = $a[4];
    $Tree{$child} = $parent;
    $Level{$child} = $level;
}

foreach my $c (keys %Tree) {
    my ($level_you_want,$domain) = find_level($c, $level);
    $Lookup{$c} = $level_you_want;
    $Domain{$c} = $domain;
    # print join("\t", $c, $level_you_want, $domain) . "\n";
}

open(IN,"<$taxassfile") || die "\n Error: Cannot open the file: $taxassfile\n\n";
while(<IN>) {
    chomp;
    my @a = split(/\t/, $_);
    my $abun  = $a[0];
    my $taxid = $a[1];
    my ($lev,$domain) = ("none","none");
    if (exists $Domain{$taxid}) { $domain = $Domain{$taxid}; }
    if (exists $Lookup{$taxid})  { $lev = $Lookup{$taxid}; }
    $Results{$lev}{$domain} += $abun;
}
close(IN);

## Need to read the results into this temporary hash to sort the results
my %Format;
foreach my $i (keys %Results) {
    foreach my $j (keys %{$Results{$i}}) {
	my $form = join("\t", $i, $j);
        $Format{$form} = $Results{$i}{$j}
    }
}

## Sort the temporary hash on the values numerically descending
open(OUT,">$outfile") || die "\n Error: Cannot write to $outfile\n";
print OUT join("\t", "#abundance", "$level", "domain" ) . "\n";
foreach my $i (sort { $Format{$b} <=> $Format{$a} } keys %Format) {
    print OUT join("\t", $Format{$i}, $i) . "\n";
}
close(OUT);

exit 0;

sub find_level
{
    my $taxid = $_[0];
    my $lev = $_[1];
    my $current_level = $Level{$taxid};

    my $level_return = "none";
    my $dom = "none";
    while($taxid != 1) {
	# print join("\t", $taxid, $current_level, $dom) . "\n";
	if ($current_level eq $lev) {
	    $level_return = $Names{$taxid};
	}
	elsif ($current_level eq "superkingdom") {
	    $dom = $Names{$taxid};
	}
	$taxid = $Tree{$taxid};
	$current_level = $Level{$taxid};
    }
    return ($level_return, $dom);
}
