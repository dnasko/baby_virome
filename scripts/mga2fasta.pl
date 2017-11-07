#!/usr/bin/perl -w

# MANUAL FOR mga2fasta.pl

=pod

=head1 NAME

mga2fasta.pl -- converts a metagene mga file to a pep and nuc fasta file

=head1 SYNOPSIS

 mga2fasta.pl --fasta=/Path/to/infile.fasta --mga=/Path/to/infile.mga --prefix=PREFIX --outdir=/Path/to/output_dir
                     [--help] [--manual]

=head1 DESCRIPTION

 This is an updated script to the old mga2seq_pep.pl. Should be faster
 and easier to install / run. Converts MetaGene files to peptide and
 nucleotide orf files.
 
=head1 OPTIONS

=over 3

=item B<-f, --fasta>=FILENAME

Input file in FASTA format. (Required) 

=item B<-mg, --mga>=FILENAME

Input file in MGA format. (Required)

=item B<-p, --prefix>=PREFIX

Prefix for output files. (Required)

=item B<-o, --outdir>=DIR

Output directory. (Required) 

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-m, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries.



=head1 AUTHOR

Written by Daniel Nasko, 
Center for Bioinformatics and Computational Biology, University of Delaware.

=head1 REPORTING BUGS

Report bugs to dnasko@udel.edu

=head1 COPYRIGHT

Copyright 2017 Daniel Nasko.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's 
usage <http://bioinformatics.udel.edu/Core/Acknowledge>.

=cut


use strict;
use Getopt::Long;
use File::Basename;
use Pod::Usage;

#ARGUMENTS WITH NO DEFAULT
my($fasta,$mga,$prefix,$outdir,$help,$manual);

GetOptions (	
    "f|fasta=s"	=>	\$fasta,
    "mg|mga=s" => \$mga,
    "p|prefix=s" => \$prefix,
				"o|outdir=s"	=>	\$outdir,
				"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required argument --fasta not found.\n\n", -exitval => 2, -verbose => 1)  if (! $fasta );
pod2usage( -msg  => "\n\n ERROR!  Required argument --mga not found.\n\n", -exitval => 2, -verbose => 1)    if (! $mga );
pod2usage( -msg  => "\n\n ERROR!  Required argument --prefix not found.\n\n", -exitval => 2, -verbose => 1) if (! $prefix );
pod2usage( -msg  => "\n\n ERROR!  Required argument --outdir not found.\n\n", -exitval => 2, -verbose => 1) if (! $outdir);

my %Fasta;
my ($h,$s,$l) = ("","",0);
my %Model = (
    'a' => 'archaea',
    'p' => 'phage',
    'b' => 'bacteria',
    's' => 'self');
my %Codon;
my $nseqs=0;
LoadCodonTable();
open(IN,"<$fasta") || die "\n\n Cannot open the input file: $fasta\n\n";
while(<IN>) {
    chomp;
    if ($_ =~ m/^>/) {
	unless ($l == 0) {
	    $s = uc($s);
	    $Fasta{$h} = $s;
	    $s = "";
	}
	$h = $_;
	$h =~ s/^>//;
	$h =~ s/ .*//;
    }
    else {
	$s = $s . $_;
    }
    $l++;
}
close(IN);
$s = uc($s);
$Fasta{$h} = $s;

# NODE_1_length_56803_cov_57.912_g0_i0_420_4775_1 size=1452 gc=0.302607 start=420 stop=4775 strand=+ frame=0 model=self score=844.516 type=complete caller=MetaGENE

open(PEP,">$outdir/$prefix.pep") || die "\n Cannot write to: $outdir/$prefix.pep\n";
open(NUC,">$outdir/$prefix.nuc") || die "\n Cannot write to: $outdir/$prefix.nuc\n";
open(IN,"<$mga") || die "\n Cannot open the mga file: $mga\n";
while(<IN>) {
    chomp;
    if ($_ =~ m/^gene/) {
	my @a = split(/\t/, $_);
	my $gid = get_gid($a[0]);
	if (exists $Fasta{$h}) {
	    my $nt_orf_seq = get_nt_orf($Fasta{$h}, $a[1], $a[2], $a[3], $a[4]);
	    if ($a[3] eq "-") { my $tmp=$a[1]; $a[1]=$a[2];$a[2]=$tmp; }
	    print NUC ">" . $h . "_" . $a[1] . "_" . $a[2] . "_" . $gid . " size=" . length($nt_orf_seq) . " gc=" . gc($nt_orf_seq) . " start=$a[1]" . " stop=$a[2]" . " strand=$a[3]" . " frame=$a[4]" . " model=" . $Model{$a[7]} . " score=$a[6]" . " type=" . get_type($a[5]) . " caller=MetaGENE" . "\n";
	    print NUC $nt_orf_seq . "\n";

	    print PEP ">" . $h . "_" . $a[1] . "_" . $a[2] . "_" . $gid . " size=" . orf_len($nt_orf_seq) . " gc=" . gc($nt_orf_seq) . " start=$a[1]" . " stop=$a[2]" . " strand=$a[3]" . " frame=$a[4]" . " model=" . $Model{$a[7]} . " score=$a[6]" . " type=" . get_type($a[5]) . " caller=MetaGENE" . "\n";
	    my $translated_seq = translate($nt_orf_seq);
	    if ($a[5] =~ m/^1/) {$translated_seq =~ s/^./M/;} ## This seems weird, but isnt. MetaGene allows for alternative start codons when predicting ORFs, funny thing about alternative start codons is, even if they are supposed to encode some other peptide the tRNA will always put a Met there.
	    print PEP $translated_seq . "\n";
	    $nseqs++;
	}
	else {die "\n Error! Cannot find the sequence: $h\n";}
    }
    elsif ($_ =~ m/^#/) {
	unless ($_ =~ m/^# gc = / || $_ =~ m/^# self: /) {
	    $h = $_;
	    $h =~ s/^# //;
	    $h =~ s/ .*//;
	}
    }
}
close(IN);
close(PEP);
close(NUC);

if ($nseqs == 0) { die "\n Error: There were no ORFs predicted from your input files\n FASTA file: $fasta\n MGA file: $mga\n\n"; }

sub translate
{
    my $str = $_[0];
    my $peptide = "";
    for (my $i=0; $i<length($str)-2; $i+=3) {
	my $mer = substr $str, $i, 3;
	if (exists $Codon{$mer}) { $peptide = $peptide . $Codon{$mer}; }
	else { $peptide = $peptide . 'X'; }
    }
    return $peptide;
}
sub LoadCodonTable
{   # Load the AA codon table from __END__ of program.
    my @data = <DATA>;
    foreach my $line (@data)
    {       chomp($line);
	    my @codons = split(/ /,$line);
	    my $AA = shift(@codons);
	    foreach my $nnn (@codons)
	    {       $nnn =~ s/U/T/g;
		    $Codon{$nnn} = $AA;
	    }
    }
}
sub revcomp
{
    my $str = $_[0];
    $str = scalar reverse $str;
    $str =~ tr/ATGC/TACG/;
    return $str;
}
sub get_type
{
    my $str = $_[0];
    if ($str eq "11") { return "complete";}
    elsif ($str eq "10") { return "lack_stop";}
    elsif ($str eq "01") { return "lack_start";}
    elsif ($str eq "00") { return "incomplete";}
}
sub gc
{
    my $str = $_[0];
    my $gcs = $str =~ tr/GCgc/GCGC/;
    my $gc_content = $gcs / length($str);
    return $gc_content;
}
sub orf_len
{
    my $str = $_[0];
    my $len = length($str)/3;
    return $len;
}
sub get_nt_orf
{
    my $seq = $_[0];
    my $start = $_[1];
    my $stop = $_[2];
    my $sense = $_[3];
    my $frame = $_[4];
    my $length = $stop - $start + 1;
    $start--;
    my $orf = substr $seq, $start, $length;
    if ($sense eq "-") { $orf = revcomp($orf); }
    $orf = substr $orf, $frame;
    return $orf;
}
sub get_gid
{
    my $str = $_[0];
    $str =~ s/gene_//;
    return $str;
}

__END__
A GCT GCC GCA GCG
R CGT CGC CGA CGG AGA AGG
N AAT AAC
D GAT GAC
C TGT TGC
Q CAA CAG
E GAA GAG
G GGT GGC GGA GGG
H CAT CAC
I ATT ATC ATA
L TTA TTG CTT CTC CTA CTG
K AAA AAG
M ATG
F TTT TTC
P CCT CCC CCA CCG
S TCT TCC TCA TCG AGT AGC
T ACT ACC ACA ACG
W TGG
Y TAT TAC
V GTT GTC GTA GTG
* TAA TAG TGA
