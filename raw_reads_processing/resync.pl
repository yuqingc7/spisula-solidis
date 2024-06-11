#!/usr/bin/perl -w

#######################################################################
# resync.pl
# John Garbe
# March 2012
#
#######################################################################

=head1 NAME

resync.pl - Resynchronize a pair of paired-end fastq files.

=head1 SYNOPSIS

resync.pl sample1_R1.fastq sample1_R2.fastq [sample1_R1_synced.fastq sample1_R2_synced.fastq]

=head1 DESCRIPTION

Programs that process paired-end fastq files usually require that the Nth read in the R1 fastq file and the Nth read in the R2 fastq file are mates. Using trimming or filtering programs that aren't paired-end aware often results in reads being removed from one paired-end fastq file but not the other, resulting in "unsyncronized" files. This program reads in two unsynchronized fastq files and writes out two synchronized fastq files. The synchronized files have properly paired reads, with singleton reads removed. Casava 1.7 and 1.8 read ID formats are supported. This program shouldn't use much memory (<1GB), but maximum memory use could be equivalent to the size of one input file in a worst-case scenario.

Options:
    -h : Display usage information

=cut

use Getopt::Std;
use Pod::Usage;

our ($opt_h);

$usage = "USAGE: resync.pl sample1_R1.fastq sample1_R2.fastq [sample1_R1_synced.fastq] [sample1_R2_synced.fastq]\n";

die $usage unless getopts('h');
pod2usage(q(-verbose) => 3) if ($opt_h);
die $usage unless ($#ARGV >= 1 && $#ARGV <= 3);

# open up input files
open F1, "<$ARGV[0]" or die "cannot open $ARGV[0]\n";
open F2, "<$ARGV[1]" or die "cannot open $ARGV[1]\n";

# open up output files
if ($#ARGV > 1) {
    open O1, ">$ARGV[2]" or die "cannot open $ARGV[2]\n";
} else {
    open O1, ">$ARGV[0].out" or die "cannot open $ARGV[0].out\n";
}
if ($#ARGV == 3) {
    open O2, ">$ARGV[3]" or die "cannot open $ARGV[3]\n";
} else {
    open O2, ">$ARGV[1].out" or die "cannot open $ARGV[1].out\n";
}

$readid = `head -n 1 $ARGV[0]`;
chomp $readid;
if ($readid =~ /^@\S+\/[12]$/) { # @ at the start of the line followed by non-whitespace, a /, a 1 or 2, the end of the line
    $id_type = 1;
    print STDERR "Casava 1.7 read id style\n";
} elsif ($readid =~ /^@\S+\W[12]\S+$/) { # @ at the start of the line followed by non-whitspace, a space, a 1 or 2, non-whitespace
    $id_type = 2;
    print STDERR "Casava 1.8 read id style\n";
} else {
    print STDERR "Cannot determine read id style\n";
    print STDOUT "Unknown id style: $readid\n";
    print STDERR "Supported read id formats: Cosava versions 1.7 and 1.8\n";
    exit 1;
}

# read in a line from the first file, then a line from the second, see if the second has a match for the first, print it out, and remove from memory
$f1readcount = 0;
$f2readcount = 0;

$f1missed = 0;
$f2missed = 0;

while ((! eof(F1)) or (! eof(F2))) {

    if (! eof(F1)) {
	$f1line1 = <F1>;
	$f1readcount++;
	if ($f1readcount % 1000000 == 0) { # print progress update
	    print STDERR "$f1readcount reads processed in first file\n";
	}
	$f1line2 = <F1>;
	$f1line3 = <F1>;
	$f1line4 = <F1>;

	if ($id_type == 1) {
	    ($id, $junk) = split /\//, $f1line1;
	} else {
	    ($id, $junk) = split ' ', $f1line1;
	}

	if (defined($r2{$id})) {
	    $f2missed--;
	    print O2 "$r2{$id}";
	    print O1 "$f1line1";
	    print O1 "$f1line2";
	    print O1 "$f1line3";
	    print O1 "$f1line4";
	    delete $r2{$id};
	} else {
	    $r1{$id} = $f1line1 . $f1line2 . $f1line3 . $f1line4;
	    $f1missed++;
	}
    }

    if (! eof(F2)) {
	$f2line1 = <F2>;
	$f2readcount++;
	if ($f2readcount % 1000000 == 0) { # print progress update
	    print STDERR "$f2readcount reads processed in second file\n";
	}
	$f2line2 = <F2>;
	$f2line3 = <F2>;
	$f2line4 = <F2>;

	if ($id_type == 1) {
	    ($id, $junk) = split /\//, $f2line1;
	} else {
	    ($id, $junk) = split ' ', $f2line1;
	}

	if (defined($r1{$id})) {
	    $f1missed--;
	    print O1 "$r1{$id}";
	    print O2 "$f2line1";
	    print O2 "$f2line2";
	    print O2 "$f2line3";
	    print O2 "$f2line4";
	    delete $r1{$id};
	} else {
	    $r2{$id} = $f2line1 . $f2line2 . $f2line3 . $f2line4;
	    $f2missed++;
	}
    }
}

print STDERR "$f1readcount total reads processed in first file\n";
print STDERR "$f2readcount total reads processed in second file\n";

print "$f1readcount reads in $ARGV[0], $f1missed without mate\n";
print "$f2readcount reads in $ARGV[1], $f2missed without mate\n";


