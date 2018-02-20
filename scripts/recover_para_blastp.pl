#!/usr/bin/perl

# MANUAL FOR recover_para_blastp.pl

=pod

=head1 NAME

recover_para_blastp.pl -- Recover results from a para_blastp run if it was interupted for some reason. MUST BE IN TABULAR FORMAT.

=head1 SYNOPSIS

 recover_para_blastp.pl --blast_tmp=/Path/to/blast_tmp_dir --remaining_fasta=/Path/to/output_remaining.fasta --current_blast_out=/Path/to/output.btab
                     [--help] [--manual]

=head1 DESCRIPTION

Pass this script the temporary para_blastp directory and it will create a FASTA for the remaining sequences you need to BLAST and gather the current resutls you have.

=head1 OPTIONS

=over 3

=item B<-bt, --blast_tmp>=DIRNAME

The para_blast temporary directory. blast_tmp_random_letters. (Required)

=item B<-rf, --remaining_fasta>=FILENAME

The remaining sequences that need to be BLAST'd in FASTA format. (Required)

=item B<-cb, --current_blast_out>=FILENAME

Where to save the current blast results you have. (Required)

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
use threads;
use FindBin;
use Cwd 'abs_path';
my $script_working_dir = $FindBin::Bin;

#ARGUMENTS WITH NO DEFAULT
my($blast_tmp,$remaining_fasta,$current_blast_out,$help,$manual);

GetOptions (	
				"bt|blast_tmp=s"	 => \$blast_tmp,
                                "rf|remaining_fasta=s"   => \$remaining_fasta,
                                "cb|current_blast_out=s" => \$current_blast_out,
             			"h|help"	=>	\$help,
				"m|manual"	=>	\$manual);

# VALIDATE ARGS
pod2usage(-verbose => 2)  if ($manual);
pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} )  if ($help);
pod2usage( -msg  => "\n\n ERROR!  Required arguments --blast_tmp not found.\n\n", -exitval => 2, -verbose => 1)  if (! $blast_tmp );
pod2usage( -msg  => "\n\n ERROR!  Required arguments --remaining_fasta not found.\n\n", -exitval => 2, -verbose => 1)     if (! $remaining_fasta );
pod2usage( -msg  => "\n\n ERROR!  Required arguments --current_blast_out not found.\n\n", -exitval => 2, -verbose => 1)    if (! $current_blast_out );

my @temp_fasta_files = `/bin/ls $blast_tmp/*.fsa`;
my $print_flag = 0;

if (scalar(@temp_fasta_files) < 1) { die "\n Error: The blast_tmp directory you passed in does not contain any split-N.fsa files. You expect these with a para_blast run. Make sur eyoure passing in the right blast_tmp directory\n"; }

for (my $i=1; $i<=scalar(@temp_fasta_files); $i++) {
    my $last_sequence = find_last("$blast_tmp/result_splits/split.$i.txt");
    if ($i==1) { open(OUT,">$remaining_fasta") || die "\n Cannot write to: $remaining_fasta\n"; }
    else { open(OUT,">>$remaining_fasta") || die "\n Cannot write to: $remaining_fasta\n"; }
    open(IN,"<$blast_tmp/split-$i.fsa") || die "\n Cannot open the file: $blast_tmp/split-$i.fsa\n";
    $print_flag = 0; ## Reset this for each file
    while(<IN>) {
	chomp;
	if ($_ =~ m/^>/) {
	    my $header = parse_header($_);
	    if ($header eq $last_sequence) {
		print OUT $_. "\n";
		$print_flag = 1;
	    }
	    elsif ($print_flag == 1) {
		print OUT $_ . "\n";
	    }
	}
	elsif ($print_flag == 1) {
	    print OUT $_ . "\n";
	}
    }
    close(IN);
    close(OUT);
}

`cat $blast_tmp/result_splits/split.1.txt.tmp > $current_blast_out; rm $blast_tmp/result_splits/split.1.txt.tmp`; ## write over results
for (my $i=2; $i<=scalar(@temp_fasta_files); $i++) {
    `cat $blast_tmp/result_splits/split.$i.txt.tmp >> $current_blast_out`; ## append to results
    `rm $blast_tmp/result_splits/split.$i.txt.tmp`;
}


exit 0;

sub parse_header
{
    my $s = $_[0];
    $s =~ s/^>//;
    $s =~ s/ .*//;
    return $s;
}

sub find_last
{
    my $f = $_[0];
    my $l = "";
    open(IN,"<$f") || die "\n Cannot open the file: $f\n";
    while(<IN>) {
	chomp;
	my @a = split(/\t/, $_);
	if (scalar(@a) < 6) { die "\n Make sure youre working with tabular blast output files...\n" }
	$l = $a[0];
    }
    close(IN);

    open(OUT,">$f.tmp") || die "\n Cannot write to: $f.tmp\n\n";
    open(IN,"<$f") || die "\n Cannot open the file: $f\n";
    while(<IN>) {
	chomp;
	my @a =split(/\t/, $_);
	unless ($a[0] eq $l) { print OUT $_ . "\n"; }
    }
    close(IN);
    close(OUT);
    return($l);
}
