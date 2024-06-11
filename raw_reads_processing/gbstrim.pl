#!/usr/bin/perl -w

#####################################
# gbstrim.pl 
# John Garbe
# June 2016
# 
#
#####################################

=head1 NAME

gbstrim.pl - 
 - Trim padding and/or cut site residue from beginning of reads
 - Trim sequencing adapter from ends of reads
 - Trim all reads to length

NOTE: this doesn't handle paired-end data

=head1 SYNOPSIS

gbstrim.pl --enzyme1 bamhi --enzyme2 psti --read R1 [--threads 1] [--minlength 20] --fastqfile file.fastq --outputfile out.fastq 

=head1 DESCRIPTION

Options:
    --fastqfile sample.fastq : fastq file of reads to process, can be gz compressed
    --outputfile sample.gbstrim.fastq : fastq file of processed reads
    --enzyme1 bamhi : First restriciton enzyme used in library creation
    --enzyme2 psti : Second restriciton enzyme used in library creation
    --read R1 : Specify if the provided fastq file is the R1 or R2 read
    --minlength 20 : discard reads shorter than minlength (default: 20)
    --maxlength 95 : trim all reads longer than maxlength to maxlength (default: readlength - length of longest padding sequence)
    --removecutsite : trim off cut sites as well
    --threads 1 : Number of cpu cores cutadapt should use (default: 1)
    --verbose : Print additional details while running
    --help : Display usage information

Advanced options:
    --r1padding C,TG,AAG,GCTC : comma-separated list of 5' padding sequences. To include a pad of zero start with a comma: ,C,TG,AAG,GCTC
    --r2padding C,TG,AAG,GCTC : comma-separated list of 3' padding sequences. To include a pad of zero start with a comma: ,C,TG,AAG,GCTC
    --r1overhang AGCT : Overhang sequence left by first restriction enzyme
    --r2overhang AGCT : Overhang sequence left by second restriction enzyme
    --adapter CTGTCTCTTATACACATCTCCGAG : sequencing primer adapter to trim from 3' end of reads
    --debug : Print output useful for debugging
=cut

##################### Initialize ###############################

use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::RealBin";

# set defaults
$r1adapter = "CTGTCTCTTATACACATCTCCGAG";
$r2adapter = "CTGTCTCTTATACACATCTGACGC";
$minlength = 20;
$uniformlength = 0;
GetOptions("help" => \$help,
	   "verbose" => \$verbose,
	   "threads=i" => \$threads,
	   "fastqfile=s" => \$fastqfile,
	   "outputfile=s" => \$outputfile,
	   "enzyme1=s" => \$r1enzyme,
	   "enzyme2=s" => \$r2enzyme,
	   "read=s" => \$read,
	   "uniformlength" => \$uniformlength,
	   "removecutsite" => \$removecutsite,

	   # advanced options
	   "r1overhang=s" => \$r1overhang,
	   "r2overhang=s" => \$r2overhang,
	   "r1padding=s" => \$r1padding,
	   "r2padding=s" => \$r2padding,
	   "debug" => \$debug,
	   "minlength=i" => \$minlength,
	   "maxlength=i" => \$maxlength,
	   "r1adapter=s" => \$r1adapter,
	   "r2adapter=s" => \$r2adapter,
    ) or pod2usage;
pod2usage(q(-verbose) => 3) if ($help);
if ($#ARGV >= 0) {
    print "Unknown commandline parameters: @ARGV\n";
    pod2usage;
}
#changed $enzyme1 > $r1enzyme to avoid typo?
if (! ($fastqfile and $outputfile and $read)) {
    print "--fastqfile, --outputfile, --read, --enzyme1 and enzyme2 are required\n";
    pod2usage;
}
die "Cannot find fastq file $fastqfile\n" if (! -e $fastqfile);

### Parameter handling ###
if (lc($read) eq "r1") {
    $read = "r1";
} elsif (lc($read) eq "r2") {
    $read = "r2";
} else {
    die "Unknown --read value $read\n";
}
# use simple commands to set advanced commands
if ($r1enzyme) {
    die "--enzyme2 is required, for single-digest datasets just list the enzyme again" unless ($r2enzyme);
#    $r2enzyme = $r1enzyme unless ($r2enzyme);
    die "Unknown enzyme1 $r1enzyme\n" unless ($enzymes{lc($r1enzyme)});
    die "Unknown enzyme2 $r2enzyme\n" unless ($enzymes{lc($r2enzyme)});
    $r1enzyme = lc($r1enzyme);
    $r2enzyme = lc($r2enzyme);
    if ($read eq "r2") {
	$tmp = $r1enzyme;
	$r1enzyme = $r2enzyme;
	$r2enzyme = $tmp;
    }
    $r1overhang = $enzymes{$r1enzyme}{overhang};
    $r2overhang = $enzymes{$r2enzyme}{overhang};
    if ($read eq "r1") {
	$r1padding = $padding{"$r1enzyme-r1"} unless ($r1padding);
	$r2padding = $padding{"$r2enzyme-r2"} unless ($r2padding);
    } else {
	$r1padding = $padding{"$r1enzyme-r2"} unless ($r1padding);
	$r2padding = $padding{"$r2enzyme-r1"} unless ($r2padding);
    }

} elsif ($r1overhang and $r2overhang and $r1padding and $r2padding) {
    # things are good
} else {
    print "--enzyme1 and --enzyme2 required (or the advanced overhang and padding parameters)\n";
    pod2usage;
}

$r1padding = " " . $r1padding if ($r1padding =~ /^,/);
$r2padding = " " . $r2padding if ($r2padding =~ /^,/);

$ifile = $fastqfile;

### Get read length ###
if ($ifile =~ /\.gz$/) {
    $readlength = `gunzip -c $ifile | sed '2q;d' | wc -m`;
} else {
    $readlength = `sed '2q;d' $ifile | wc -m`;
}
chomp $readlength;
$readlength--;
print "Detected read length: $readlength\n";
@padding = split /,/, $r1padding;
$longestpad = 0;
foreach $pad (@padding) {
    $length = length($pad);
    $longestpad = $length if ($length > $longestpad);
}
$maxlength = $readlength - $longestpad unless (defined($maxlength));
print "The longest padding sequence is $longestpad bp long\n" if ($verbose);


### Identify phasing padding and cut site overhang ###
$r1pattern = &adapterpattern($r1overhang, $r1padding);
$r2pattern = &adapterpattern($r2overhang, $r2padding);

### Trim any trailing adapter sequence ###
if ($read eq "r1") {
    $adapter = $r1adapter;
} else {
    $adapter = $r2adapter;
}
if ($adapter) {
    
    print "Trimming $adapter from 3' end\n" if ($adapter);

    $result = `which cutadapt`;
    if ($verbose) {
	print "cutadapt binary: $result";
	$version = `cutadapt --version`;
	print "cutadapt version: $version";
    }
    die "cutadapt is required for 3' adapter trimming\n$result\n" if ($?);

    $ofile = $outputfile . ".3trim.fastq";

    $adapttrim = ($adapter) ? "--adapter=$adapter" : "";
    $cutthreads = ($threads) ? "--cores=$threads" : "";

    $command = "cutadapt $cutthreads $adapttrim -o $ofile $ifile";
    print "$command\n" if ($verbose);
    $result = `$command`;
    die "cutadapt failed: $result" if ($?);
    print $result if ($verbose);
    @result = split /\n/, $result;

    # parse stats from output
    foreach $result (@result) {
	if ($result =~ /^Total reads processed:\s+([\d,]+)/) {
	    $total = $1;
	    $total =~ s/,//g;
#           print "Total reads\t$total\n";
	    next;
	}
	if ($result =~ /^Reads with adapters:\s+([\d,]+)/) {
	    $good = $1;
	    $good =~ s/,//g;
	    $stats{readswithtrailingadapter} = $good;
#	    print "readswithtrailingadapter\t$good\n";
	    $percent = &round100($good / $total * 100);
	    $stats{"pct-readswithtrailingadapter"} = $percent;
#	    print "pct-readswithtrailingadapter\t$percent\n";
	    next;
	}
    }
    $ifile = $ofile;
}

# trim padding/overhang sequence from 3' end of reads that are shorter than the readlength, and from 5' end of all reads
if ($r1overhang and $r1padding) {

    print "Trimming from 5' end using search pattern $r1pattern\n";
    print "Trimming from 3' end using search pattern $r2pattern\n";
    print "Trimming all reads down to $maxlength bases\n" if ($maxlength);
    print "Discarding trimmed reads shorter than $minlength bases\n" if ($minlength);

    $stats{tooshort} = 0;
    $stats{no3cutsite} = 0;
    $stats{no5cutsite} = 0;

    $ofile = $outputfile . ".5trim.fastq";
    die "Cannot find input file $ifile\n" unless (-e $ifile);

    if ($ifile =~ /\.gz$/) {
	open IFILE, "gunzip -c $ifile |" or die "cannot open gzip pipe of $ifile: $!\n";	
    } else {
	open IFILE, $ifile or die "cannot open $ifile: $!\n";
    }
    open OFILE, ">$ofile" or die "cannot open $ofile: $!\n";
    open OFILE2, ">$outputfile.nocutsite" or die "cannot open $ofile.untrimmed: $!\n";
    open OFILE3, ">$outputfile.no3cutsite" or die "cannot open $ofile.untrimmed: $!\n";

    while ($id = <IFILE>) {
	chomp $id;
	$seq = <IFILE>;
	chomp $seq;
	my $plus = <IFILE>;
	chomp $plus;
	$qual = <IFILE>;
	chomp $qual;

	$stats{inputreads}++;

	# trim 3' end if cutadapt saw and removed an 3' adapter
	if (length($seq) < $readlength) {
	    $qes = &reverse_complement($seq);
	    if (@match = $qes =~ /$r2pattern/) {
		$thispadding = $1 // "";
		$thisoverhang = $2;
		$thisadapter = $thispadding;
		if ($removecutsite) {
		    $trimlength = length($thisadapter . $thispadding);
		} else {
		    $trimlength = length($thisadapter);
		}
		if ($trimlength > 0) {
		    $seq = substr($seq,0,-$trimlength);
		    $qual = substr($qual,0,-$trimlength);
		}
		$r2cutsitecount{$thisadapter}++;
	    } else {
		print OFILE3 "$id\n$seq\n$plus\n$qual\n";	    
		$stats{no3cutsite}++;
		next;
	    }
	}
	# trim 5' end
	if (@match = $seq =~ /$r1pattern/) {
	    $thispadding = $1 // "";
	    $thisoverhang = $2;
	    $thisadapter = $thispadding;
	    if ($removecutsite) {
		$trimlength = length($thisadapter . $thispadding);
	    } else {
		$trimlength = length($thisadapter);
	    }
	    if ($maxlength) {
		$seq = substr($seq,$trimlength,$maxlength);
		$qual = substr($qual, $trimlength,$maxlength);
	    } else {
		$seq = substr($seq,$trimlength);
		$qual = substr($qual, $trimlength);
	    }
	    $cutsitecount{$thisadapter}++;
	    if (length($seq) < $minlength) {
		$stats{tooshort}++;
		next;
	    }
	    print OFILE "$id\n$seq\n$plus\n$qual\n";
	} else {
	    print OFILE2 "$id\n$seq\n$plus\n$qual\n";	    
	    $stats{no5cutsite}++;
	}
    }
    close OFILE;
    close OFILE2;

    # identify untrimmed adapters 
    if ($verbose && (-e "$outputfile.nocutsite")) {
	`cat $outputfile.nocutsite | sed -n '2~4p' | cut -c 1-11 | sort | uniq -c | sort -n | tail -n 14 > $outputfile.nocutsite.adapters`;
	print "See $outputfile.nocutsite.adapters for untrimmed adapters\n";;
    }


    if ($debug) {
	print "3' sites\n";
	foreach $sequence (sort keys %r2cutsitecount) {
	    $r2cutsitecount{$sequence} = 0 unless ($r2cutsitecount{$sequence});
	    $r2percents{$sequence} = &round100($r2cutsitecount{$sequence} / $stats{inputreads} * 100);
	    print "%$sequence\t$r2percents{$sequence}\n";
	}	
    }

    # compute stats
    $stats{"pct-no3cutsite"} = &round100($stats{no3cutsite} / $stats{inputreads} * 100);
    $stats{"pct-no5cutsite"} = &round100($stats{no5cutsite} / $stats{inputreads} * 100);
    $stats{"pct-tooshort"} = &round100($stats{tooshort} / $stats{inputreads} * 100);
    foreach $sequence (sort keys %cutsitecount) {
	$cutsitecount{$sequence} = 0 unless ($cutsitecount{$sequence});
	$percents{$sequence} = &round100($cutsitecount{$sequence} / $stats{inputreads} * 100);
    }
    $ifile = $ofile;

    foreach $sequence (sort keys %percents) {
	$pad = $sequence;
	$pad = "no_padding" if ($sequence eq "");
	$stats{"pct-" . "$pad"} = $percents{$sequence};
    }

    $removed = $stats{no3cutsite} + $stats{no5cutsite} + $stats{tooshort};
    $pctremoved = 0;
    $pctremoved = &round100($removed / $stats{inputreads} * 100) if ($stats{inputreads} > 0);
    $stats{"pct-removed"} = $pctremoved;
    $stats{"removed"} = $removed;
}

# print out stats
print "STATS Start\n";
foreach $key (sort (keys %stats)) {
    print "$key\t$stats{$key}\n";
}
print "STATS End\n";

### name output file properly, clean up tmp files ###
`mv $ofile $outputfile`;
#print "mv $ofile $outputfile\n" if ($verbose);
foreach $file ("$outputfile.3trim.fastq","$outputfile.5trim.fastq") {
    `rm $file` if (-e $file);
}

# basic error check
if ($pctremoved > 20) {
    print "WARNING: GBStrim discarded ${pctremoved}% the reads. Higher than 20% is unusual. Reads are removed if no 5' cut site is found, if sequencing adapter is removed from the 3' end but no 3' cut site is found, or if the trimmed read is too short (< $minlength bp)\n";
} elsif ($pctremoved > 10) {
    print "GBStrim discarded ${pctremoved}% the reads. Less than 10% is good. Less than 20% is OK. Higher than 20% is unusual. Reads are removed if no 5' cut site is found, if sequencing adapter is removed from the 3' end but no 3' cut site is found, or if the trimmed read is too short (< $minlength bp)\n";
} else {
    print "GBStrim discarded ${pctremoved}% of the reads. Less than 10% is good.\n";
}

exit;

########################################################

# Round to neaerest hundreth
sub round100 {
    return sprintf("%.2f", $_[0]);
}

# Given an ezyme cut site overhang and padding sequences, generate a big regular expression to match it
sub adapterpattern {

    my ($overhang, $padding) = @_;
    $pattern = &expand($overhang);

    # Add padding sequences to pattern
    @paddings = split /,/, $padding;
    if ($paddings[0] eq " ") {
	$option = "?";
	shift @paddings;
    } else {
	$option = "";
    }
    $paddingpattern = join "|", @paddings;
    $pattern = "^($paddingpattern)$option($pattern)";

    return $pattern;
}

# Create regex rules to match degenerate bases
sub expand {
    my ($pat) = @_;

    $pat =~ s/N|X/./g;

    my $PURINES      = 'AG';
    my $PYRIMIDINES  = 'CT';

    ## Avoid nested situations: [ya] --/--> [[ct]a]
    ## Yet correctly deal with: sg[ya] ---> [gc]g[cta]
    if($pat =~ /\[\w*[RYSWMK]\w*\]/) {
	$pat =~ s/\[(\w*)R(\w*)\]/\[$1$PURINES$2\]/g;
	$pat =~ s/\[(\w*)Y(\w*)\]/\[$1$PYRIMIDINES$2\]/g;
	$pat =~ s/\[(\w*)S(\w*)\]/\[$1GC$2\]/g;
	$pat =~ s/\[(\w*)W(\w*)\]/\[$1AT$2\]/g;
	$pat =~ s/\[(\w*)M(\w*)\]/\[$1AC$2\]/g;
	$pat =~ s/\[(\w*)K(\w*)\]/\[$1GT$2\]/g;
	$pat =~ s/\[(\w*)V(\w*)\]/\[$1ACG$2\]/g;
	$pat =~ s/\[(\w*)H(\w*)\]/\[$1ACT$2\]/g;
	$pat =~ s/\[(\w*)D(\w*)\]/\[$1AGT$2\]/g;
	$pat =~ s/\[(\w*)B(\w*)\]/\[$1CGT$2\]/g;
	$pat =~ s/R/\[$PURINES\]/g;
	$pat =~ s/Y/\[$PYRIMIDINES\]/g;
	$pat =~ s/S/\[GC\]/g;
	$pat =~ s/W/\[AT\]/g;
	$pat =~ s/M/\[AC\]/g;
	$pat =~ s/K/\[GT\]/g;
	$pat =~ s/V/\[ACG\]/g;
	$pat =~ s/H/\[ACT\]/g;
	$pat =~ s/D/\[AGT\]/g;
	$pat =~ s/B/\[CGT\]/g;
    } else {
	$pat =~ s/R/\[$PURINES\]/g;
	$pat =~ s/Y/\[$PYRIMIDINES\]/g;
	$pat =~ s/S/\[GC\]/g;
	$pat =~ s/W/\[AT\]/g;
	$pat =~ s/M/\[AC\]/g;
	$pat =~ s/K/\[GT\]/g;
	$pat =~ s/V/\[ACG\]/g;
	$pat =~ s/H/\[ACT\]/g;
	$pat =~ s/D/\[AGT\]/g;
	$pat =~ s/B/\[CGT\]/g;
    }
    $pat =~ s/\((.)\)/$1/g;  ## Doing these last since:
    $pat =~ s/\[(.)\]/$1/g;  ## Pattern could contain [y] (for example)
    
    return $pat;
}

sub reverse_complement {
    my $dna = shift;

    # reverse the DNA sequence
    my $revcomp = reverse($dna);

    # complement the reversed DNA sequence
    $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $revcomp;
}

sub complement {
    my $dna = shift;

    # complement the reversed DNA sequence
    $dna =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $dna;
}

sub swap {
    ($a, $b) = @_;
    my $tmp = $a;
    $a = $b;
    $b = $tmp;
}    

# imported from gopher-pipelines enzymes.pm
BEGIN {
%enzymes = (
"apeki" => {
name => "ApeKI",sequence => "G|CWGC",overhang => "CWGC",},
"psti" => {
name => "PstI",sequence => "CTGCA|G",overhang => "TGCAG",},
"bstyi" => {
name => "BstYI",sequence => "R|GATCY",overhang => "GATCY",},
"bfuci" => {
name => "BfuCI",sequence => "|GATC",overhang => "GATC",},
"mboi" => {
name => "MboI",sequence => "|GATC",overhang => "GATC",},
"nsii" => {
name => "NsiI",sequence => "ATGCA|T",overhang => "TGCAT",},
"nspi" => {
name => "NspI",sequence => "RCATG|Y",overhang => "CATGY",},
"ecot22i" => {
name => "EcoT22I",sequence => "ATGCA|T",overhang => "TGCAT",},
"apai" => {
name => "ApaI",sequence => "GGGCC|C",overhang => "GGCCC",},
"taqai" => {
name => "TaqaI",sequence => "T|CGA",overhang => "CGA",},
"taqalphai" => {
name => "TaqalphaI",sequence => "T|CGA",overhang => "CGA",},
"btgi" => {
name => "BtgI",sequence => "C|CRYGG",overhang => "CRYGG",},
"taqi" => {
name => "TaqI",sequence => "T|CGA",overhang => "CGA",},
"pasi" => {
name => "PasI",sequence => "CC|CWGGG",overhang => "CWGGG",},
"bamhi" => {
name => "BamHI",sequence => "G|GATCC",overhang => "GATCC",},
"mspi" => {
name => "MspI",sequence => "C|CGG",overhang => "CGG",},
"ncoi" => {
name => "NcoI",sequence => "C|CATGG",overhang => "CATGG",},
"sbfi" => {
name => "SbfI",sequence => "CCTGCA|GG",overhang => "TGCAGG",},
"apali" => {
name => "ApaLI",sequence => "G|TGCAC",overhang => "TGCAC",},
"bglii" => {
name => "BglII",sequence => "A|GATCT",overhang => "GATCT",},
"clai" => {
name => "ClaI",sequence => "AT|CGAT",overhang => "CGAT",},
);

$tgca_r1 = ",A,TG,CAG,GCTA,ACAGA,CTCCGA,GAATCGA,AGGTCCAA,TCACGCTCA";
#$tgca_r2 = ",C,AC,GAC,AGAC,GAGAC,CGAGAC,ACGAGAC,CACGAGAC,CCACGAGAC";
$tgca_r2 = ",G,CG,TCA,ATCT,TCGAA,ACTCTA,CAGATCA,GAATCGCA,ATCGCCACA";

$gatc_r1 = ",C,TC,ATT,CATT,GTATC,CGCATC,TCGGGTT,GATAAGCT,ACAAGTCGT,GTCCAATCGT";
$gatc_r2 = ",G,AG,TCA,AAGT,ACGAA,ACTCTG,GTACGGT,TTCGACAT,CGATGTGCT";

$catg_r1 = ",C,TC,ATC,GGTC,GTATC,CGTATC,TCGGATC";
$catg_r2 = ",T,AC,TGT,AGTC,GCGTC,CTAGAC,GTCAGTC";

$cwg_r1 = "C,TC,ATC,GATC,GGATC,CGCATC,TCGTATC,AGGCTATC,CGTGTT,TTGTGTT,AAGCTATG";
$cwg_r2 = "C,TC,GTC,AGTC,ATGAC,GACGAT,GAAGTAC,ATGATCGT,ATGGT,GATGTT,TGATATT,GAATAGCGT";

$cryg_r1 = ",C,TC,ATC,GGTC,GTATC,CGTATC,TCGGATC,ATG,GGTT,GTATG,TCGGATG,GATT,GTAAC,TCAGAAA,GATC,CGTAAC";
$cryg_r2 = ",T,AC,TGT,AGTC,GCGTC,CTAGAC,GTCAGTC,AGTG,GCTGT,CTGGAC,GTCAGTT,TAC,GCGTA,GTCAGTA,AGAC,GCGAA,GTCACCA";

$cg_r1 = "T,TG,ACG,CGTA,AACTA,GTACTA,TCAGCTA,CGTCACTA,GACTAGCTA,CTGAGCACAA";
$cg_r2 = ",";

$catc_r1 = "A,TG,CAG,GCTA,ACAGA,CTCCGA,GAATCGA,AGGTCCAA,TCACGCTCA";
#$catc_r2 = ",";

%padding = (
"nsii-r1"      => $tgca_r1,
"psti-r1"      => $tgca_r1,
"sbfi-r1"      => $tgca_r1,
"apali-r1"     => $tgca_r1,
"taqi-r1"      => $tgca_r1,
"taqai-r1"     => $tgca_r1,
"taqalphai-r1" => $tgca_r1,
"ecot22i-r1"   => $tgca_r1,
"nsii-r2"      => $tgca_r2,
"psti-r2"      => $tgca_r2,
"sbfi-r2"      => $tgca_r2,
"apali-r2"     => $tgca_r2,
"taqi-r2"      => $tgca_r2,
"taqai-r2"     => $tgca_r2,
"taqalphai-r2" => $tgca_r2,
"ecot22i-r2"   => $tgca_r2,

"bamhi-r1" => $gatc_r1,
"bstyi-r1" => $gatc_r1,
"bglii-r1" => $gatc_r1,
"bfuci-r1" => $gatc_r1,
"mboi-r1"  => $gatc_r1,
"bamhi-r2" => $gatc_r2,
"bstyi-r2" => $gatc_r2,
"bglii-r2" => $gatc_r2,
"bfuci-r2" => $gatc_r2,
"mboi-r2"  => $gatc_r2,

"ncoi-r1" => $catg_r1,
"ncoi-r2" => $catg_r2,

"pasi-r1"  => $cwg_r1,
"apeki-r1" => $cwg_r1,
"pasi-r2"  => $cwg_r2,
"apeki-r2" => $cwg_r2,

"btgi-r1" => $cryg_r1,
"btgi-r2" => $cryg_r2,

"mspi-r1" => $cg_r1,
"mspi-r2" => $cg_r2,

"nspi-r1" => $catc_r1,

"clai-r2" => ",",

);
}

