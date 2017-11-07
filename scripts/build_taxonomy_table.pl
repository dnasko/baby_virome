#!/usr/bin/perl

# MANUAL FOR build_taxonomy_table.pl

=pod

=head1 NAME

build_taxonomy_table.pl -- Created a txt file of taxonomy

=head1 SYNOPSIS

 build_taxonomy_table.pl --nodes=/Path/to/nodes.dmp --names=/Path/to/names.dmp --out=/Path/to/output.txt
                     [--help] [--manual]

=head1 DESCRIPTION

 Dumps all taxonomy.
 
=head1 OPTIONS

=over 3

=item B<-no, --nodes>=FILENAME

NCBI nodes.dmp file. (Required) 

=item B<-na, --names>=FILENAME

NCBI names.dmp file. (Required)

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

=cut


use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Pod::Usage;

#ARGUMENTS WITH NO DEFAULT
my($nodes,$names,$outfile,$help,$manual);

GetOptions (	
                                "no|nodes=s"	=>	\$nodes,
                                "na|names=s"    =>      \$names,
				"o|out=s"	=>	\$outfile,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --nodes not found.\n\n", -exitval => 2, -verbose => 1)  if (! $nodes );
pod2usage( -msg  => "\n\n ERROR!  Required argument --names not found.\n\n", -exitval => 2, -verbose => 1)  if (! $names );
pod2usage( -msg  => "\n\n ERROR!  Required argument -outfile not found.\n\n", -exitval => 2, -verbose => 1)  if (! $outfile);

my %Names;
my %Tree;

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

open(OUT,">$outfile") || die "\n Cannot open the file: $outfile\n";
foreach my $child (keys %Tree) {
    my $parent = $Tree{$child};
    my @Taxa = get_lineage($child);
    # if ($Taxa[0] eq "Viruses") {
	print OUT $child . "\t" . join("\t", @Taxa) . "\n";
    # }
}
close(OUT);

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

exit 0;
