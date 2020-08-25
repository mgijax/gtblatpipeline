#!/usr/bin/perl
#
# psl_gff_GT1.pl
#
# PURPOSE: 
#	Partitions a PSL file based on whether the GT sequence
#	hit once or multiple times.
#
# AUTHOR:
#	Bob Sinclair
#
#use strict;
#use warnings;

my($pslfile, $source) = @ARGV;
my $size = @ARGV;

# note lines 81 and 132 have code specific to simpify the genetrap headers

unless ($size eq 2) {
	# use of the program and expectations follows!
	print "psl_gff_GT1.pl is specific to the GEne Trap GU project - see code lines ~81 and ~132 simplifying the seqID\n";
	print "\nthis program, psl_gff_(ver)_u.pl needs two command line arguments\n\n";
	print "the first is the name of the .psl file\n";
	print "the second is a text string for the 'source' of the data (GFF col 2)\n";
	print "if you're not sure about the source\nuse your name/initials and a data origin linked with an underscore\n";
	print "keeping this unique will help if the file is needed for a separate GBrowse track\n\n";
	print "the program will produce two types of data files: input_file.gff and input_file_Gbrowse.gff\n";
	print "the first file will have just the Qname (query sequence) in column 9\n";
	print "the second will also have the gene/mrna.exon information needed by GBrowse\n\n";
	print "two versions of each of these files is produced, one with seqIDs that occur once in the input file and a second with multiples\n";
	print "it also prints some lists of seqIDs as well as a file of notes on the run\n";
	print "ask Bob about list-comparison scripts of missing seqIDs are important!\n";
	print "eg typing the following at the command line\n";
	print "\tperl psl_gff_(ver)_u.pl myfile.psl BOB_TEST\n";
	print "will produce single adn multiple insances of myfile.gff (exons only) \nand myfile_Gbrowse.gff (exons plus gene and mRNA rows).";
	print "version 7 of this prog sequentially numbers the instances of each multiply-hitting seqID\n";
	exit;
}

$pslfile =~ /(.*)\.psl/;
my $root = $1;
my $notefile = $root.'_notes.txt';
my $singleexonfile = $root.'_single.gff';
my $singlegbrowsefile = $root.'_single_Gbrowse.gff';
my $multiexonfile = $root.'_multi.gff';
my $multigbrowsefile = $root.'_multi_Gbrowse.gff';
my $listfile = $root.'_seqIDs.txt';

print "reading $pslfile and writing to $singleexonfile and $singlegbrowsefile (plus multi files)\n";

open(SINGLEEXONS, ">$singleexonfile") or die "can't open $singleexonfile for writing: $!\n";
open(SINGLEGBROWSE, ">$singlegbrowsefile") or die "can't open $singlegbrowsefile for writing: $!\n";
open(MULTIEXONS, ">$multiexonfile") or die "can't open $multiexonfile for writing: $!\n";
open(MULTIGBROWSE, ">$multigbrowsefile") or die "can't open $multigbrowsefile for writing: $!\n";

open(NOTES, ">$notefile") or die "can't open $notefile for writing: $!\n";
my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime;
my $thisday = ('Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat')[(localtime)[6]];
my $thismonth = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')[(localtime)[4]];
print NOTES "$thisday $mday $thismonth ", $year+1900, " $hour : $min : $sec\n";
print NOTES "reading $pslfile and writing to $singleexonfile and $singlegbrowsefile (plus multi files\n";

open(PSL, "<$pslfile") or die "can't open $pslfile :$!\n";
my $line;
my($chr, $start, $end, $match, $mismatch, $repmatch, $Ns, $QgapCount, $QgapBases, $TgapCount, $TgapBases, $strand, $Qname, $Qsize, $Qstart, $Qend, $Tname, $Tsize, $Tstart, $Tend, $blockCount, $blockSizes, $qStarts, $tStarts);
my $score = '.';
my $phase = '.';
my $attributes = '';
my $code_info = '';
my @blockSizes;
my @Tstarts;
#my($dummy1, $dummy2, $dummy3, $dummy4);
my $seqID = '';
my $numsinglegenes = 0;
my $numsingleexons = 0;
my $nummultigenes = 0;
my $nummultiexons = 0;
my $inputrows = 0;

my %seqIDs;


# this section increments a counter for each seqID
while ($line = <PSL>) {
	if ($line =~ /psLayout/) {
		print "Header section found, re-run BLAT with -nohead option or manually delete the 5 header lines\n";
		exit;
	}
	chomp $line;
	$seqID = '';
	($match,  $mismatch,  $repmatch,  $Ns,  $QgapCount,  $QgapBases,  $TgapCount,  $TgapBases,  $strand,  $Qname,  $Qsize,  $Qstart,  $Qend,  $Tname,  $Tsize,  $Tstart,  $Tend,  $blockCount,  $blockSizes,  $qStarts,  $tStarts) = split /\t/, $line;
	$Qname =~ /\|gb\|(.*)\..\|/;
	$seqID = $1;
	$seqIDs{$seqID} ++;
}

#while (my ($key, $value) = each %seqIDs) {
#	print LIST "$key => $value\n"; # this line for seqIDs only
	# print LIST "$key\t$value\n"; # this line for seqID tab delim to number of hits
#}

my $multilistfile = $root.'_seqIDs_multi.txt';
my $singlelistfile = $root.'_seqIDs_single.txt';

open(M_LIST, ">$multilistfile");
open(S_LIST, ">$singlelistfile");
open(LIST, ">$listfile") or die "can't open $listfile for writing: $!\n";

while (my ($key, $value) = each %seqIDs) {
	print LIST "$key\n";
	if ($value gt 1) {
		print M_LIST "$key\n";
	} else {
		print S_LIST "$key\n";
	}
}
#	print LIST "$key\n"; # this line for seqIDs only
	# print LIST "$key\t$value\n"; # this line for seqID tab delim to number of hits
#}


#open(DEBUG, ">psl_gff_06_DEBUG.txt");

close PSL;
open(PSL, "<$pslfile") or die "can't open $pslfile :$!\n";
my %count_multi_seqID;

while ($line = <PSL>) {
	if ($line =~ /psLayout/) {
		print "Header section found, re-run BLAT with -nohead option or manually delete the 5 header lines\n";
		exit;
	}
	chomp $line;
	$inputrows ++;
	$chr = '';
	$seqID = '';
	$start = '';
	$end = '';
	($match,  $mismatch,  $repmatch,  $Ns,  $QgapCount,  $QgapBases,  $TgapCount,  $TgapBases,  $strand,  $Qname,  $Qsize,  $Qstart,  $Qend,  $Tname,  $Tsize,  $Tstart,  $Tend,  $blockCount,  $blockSizes,  $qStarts,  $tStarts) = split /\t/, $line;
	$Tname =~ /chr(.+)/;
	$chr = $1;
	$start = $Tstart;
	$end = $Tend -1;
	$Qname =~ /\|gb\|(.*)\..\|/;
	$seqID = $1;
	$count_multi_seqID{$seqID} ++; # this increments the counter for the current seqID
	#print DEBUG "after increment Qname $Qname\t seqID $seqID\tcount_multi_seqID{$seqID} $count_multi_seqID{$seqID}\n"; # seems right here

	if ($seqIDs{$seqID} gt 1) {
		$seqID .= '_'.$count_multi_seqID{$seqID}; # this will append the incremental counter to the seqID, not just the total times found
		#print DEBUG "new seqID $seqID\n";
		# this prints the gene and mRNA lines for a GBrowse-compatible file
		print MULTIGBROWSE "chr$chr\t$source\tgene\t$start\t$end\t$score\t$strand\t$phase\tGene $seqID\n";
		print MULTIGBROWSE "chr$chr\t$source\tmRNA\t$start\t$end\t$score\t$strand\t$phase\tmRNA  $seqID\n";
		$nummultigenes ++;
		@blockSizes = split /,/, $blockSizes;
		@Tstarts = split /,/, $tStarts;

		# This prints exon lines for a GBrowse-compatible file and exons for a simple gff exon-only file
		for (my $i = 0; $i < $blockCount; $i ++) {
			$Tend = $Tstarts[$i] + ($blockSizes[$i] -1);
			print MULTIGBROWSE "chr$chr\t$source\texon\t$Tstarts[$i]\t$Tend\t$score\t$strand\t$phase\tmRNA $seqID\n";
			print MULTIEXONS "chr$chr\t.\texon\t$Tstarts[$i]\t$Tend\t$score\t$strand\t$phase\t$seqID\n";
			$nummultiexons ++;
		} # end for loop 


	} else {

		# this prints the gene and mRNA lines for a GBrowse-compatible file
		# print just the seqID - a lack of underscore-number indicates single hit.
		print SINGLEGBROWSE "chr$chr\t$source\tgene\t$start\t$end\t$score\t$strand\t$phase\tGene $seqID\n";
		print SINGLEGBROWSE "chr$chr\t$source\tmRNA\t$start\t$end\t$score\t$strand\t$phase\tmRNA  $seqID\n";
		$numsinglegenes ++;
		@blockSizes = split /,/, $blockSizes;
		@Tstarts = split /,/, $tStarts;

		# This prints exon lines for a GBrowse-compatible file and exons for a simple gff exon-only file
		for (my $i = 0; $i < $blockCount; $i ++) {
			$Tend = $Tstarts[$i] + ($blockSizes[$i] -1);
			print SINGLEGBROWSE "chr$chr\t$source\texon\t$Tstarts[$i]\t$Tend\t$score\t$strand\t$phase\tmRNA $seqID\n";
			print SINGLEEXONS "chr$chr\t.\texon\t$Tstarts[$i]\t$Tend\t$score\t$strand\t$phase\t$seqID\n";
			$numsingleexons ++;
		} # end for loop 
	} # end if seqID > 1 else loop

} # end while $line = <PSL>

print "there were $inputrows lines in $pslfile\n";
print "$numsinglegenes sequences occur once and contain a total of $numsingleexons exons\n";
print "$nummultigenes sequences occur >1 time and contain a total of $nummultiexons exons\n";
print NOTES "there were $inputrows lines in $pslfile\n";
print NOTES "$numsinglegenes sequences occur once and contain a total of $numsingleexons exons\n";
print NOTES "$nummultigenes sequences occur >1 time and contain a total of $nummultiexons exons\n";

close PSL;
close SINGLEEXONS;
close SINGLEGBROWSE;
close MULTIEXONS;
close MULTIGBROWSE;
close NOTES;

