#!/usr/bin/perl -w

#######################################################################################
##
##   Insider - integration & cleavage site detection tool							 
##
##   Copyright (C) 2015-2020 Ignas Bunikis, Adam Ameur							 
##																			 
##    This program is free software: you can redistribute it and/or modify    
##    it under the terms of the GNU General Public License as published by    
##    the Free Software Foundation, either version 3 of the License, or       
##    (at your option) any later version.									 
##																			 
##    This program is distributed in the hope that it will be useful,         
##    but WITHOUT ANY WARRANTY; without even the implied warranty of 		 
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 			 
##    GNU General Public License for more details.							 
##																			 
##    You should have received a copy of the GNU General Public License       
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
#######################################################################################

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my $insiderVersion = "version 1.9";

## Read command line arguments
my ($filename, $outfile, $offset, $minPeak, $minSupport, $minClipLen, $minMappingQV, $maxMappingMM, $minAlignLen, $silent, $singlePeaks, $printReads, $ots) = get_args_and_error_check();

unless($silent){
	print STDOUT "\n************ Running Insider $insiderVersion **********\n\n";
}

my %validReads = ();
my $allReadsNr = 0;
my $validReadsNr = 0;
my %uniqPosReadsPlus = ();
my %uniqPosReadsMinus = ();
my %peaks = ();
my %peakSupport = ();
my %uniqPosReadsPlusAbMin = ();
my %uniqPosReadsMinusAbMin = ();
my $printchr = "";

unless($silent){
	print STDOUT "Reading data from $filename...\n";
}

open(FILE, "<", $filename) or die "cannot open: $!";

## Filter and sort the reads
while (my $row = <FILE>) {
	chomp $row;

	if(!($row =~ /^@/)){
		my @samline = split("\t", $row);
		my $cigar = $samline[5];
		$allReadsNr++;

		if (!($samline[1] & 4)) {
			my $samlineOK = 0;
			if($ots){
				$samlineOK = 1;
			}
			else{
				if((($cigar =~ m/^(\d+)S/) && ($1 >= $minClipLen)) || (($cigar =~ m/(\d+)S$/) && ($1 >= $minClipLen))) {
					$samlineOK = 1;
				}
			}
			if($samlineOK){
				if (($samline[4] >= $minMappingQV) && (alignQV($row) <= $maxMappingMM)) {
					my $chr = $samline[2];
					my $start_pos = $samline[3]-1;
					my $end_pos = endPosition($cigar, $start_pos);
					my $strand = readStrand($samline[1]);
					$validReads{$samline[0]} = [$samline[2], $strand, $start_pos, $end_pos, $samline[5], $samline[9]];
					$validReadsNr++;

					## Create empty hash for a chromosome..
					if(!defined($uniqPosReadsPlus{$chr})){
						$uniqPosReadsPlus{$chr} = ();
					}
					if(!defined($uniqPosReadsMinus{$chr})){
						$uniqPosReadsMinus{$chr} = ();
					}

					my $alignment_len = abs $end_pos-$start_pos;
					
					if($alignment_len > $minAlignLen){
						
						if($strand eq "+"){
							push @{$uniqPosReadsPlus{$chr}{$end_pos}}, $samline[0];
							if($ots){
								push @{$uniqPosReadsPlus{$chr}{$start_pos}}, $samline[0];
							}
						}
						if($strand eq "-"){
							push @{$uniqPosReadsMinus{$chr}{$start_pos}}, $samline[0];
							if($ots){
								push @{$uniqPosReadsMinus{$chr}{$end_pos}}, $samline[0];
							}
						}
					}
				}
				
			}
		}
	}
}

close(FILE);

unless($silent){
	print STDOUT " - $validReadsNr/$allReadsNr reads met criteria for further analysis\n\n";
	if($ots){
		print STDOUT "Searching for CRISPR-Cas9 cleavage sites...\n"
	}
	else{
		print STDOUT "Searching for sites of integrated DNA...\n"
	}
}

## Filter Plus peaks by coverage
foreach my $chr (sort keys %uniqPosReadsPlus) {

	foreach my $plusPos (sort {$a<=>$b} keys %{$uniqPosReadsPlus{$chr}}){

		my $plusHeight = scalar @{$uniqPosReadsPlus{$chr}{$plusPos}};
		
		if($plusHeight > $minPeak){
			@{$uniqPosReadsPlusAbMin{$chr}{$plusPos}} = @{$uniqPosReadsPlus{$chr}{$plusPos}};
		}
	}
}

## Filter Minus peaks by coverage
foreach my $chr (sort keys %uniqPosReadsMinus) {

	foreach my $minusPos (sort {$a<=>$b} keys %{$uniqPosReadsMinus{$chr}}){

		my $minusHeight = scalar @{$uniqPosReadsMinus{$chr}{$minusPos}};

		if($minusHeight > $minPeak){
			@{$uniqPosReadsMinusAbMin{$chr}{$minusPos}} = @{$uniqPosReadsMinus{$chr}{$minusPos}};
		} 

	}
}

## Loop through clipped peaks on plus strand...
foreach my $chr (sort keys %uniqPosReadsPlus) {

	foreach my $plusPos (sort {$a<=>$b} keys %{$uniqPosReadsPlus{$chr}}){

		my $plusHeight = scalar @{$uniqPosReadsPlus{$chr}{$plusPos}};

		if($ots){
			if(defined($uniqPosReadsMinus{$chr}{$plusPos})){
				$plusHeight = scalar @{$uniqPosReadsMinus{$chr}{$plusPos}} + $plusHeight;
			}
		}
				
		if($plusHeight > $minPeak){
			
			## Check if there is some peak on minus strand in a window surrounding plus strand peak..
			for(my $coord = $plusPos-$offset; $coord<=$plusPos+$offset; $coord++){
				if(defined($uniqPosReadsMinus{$chr}{$coord})){
					my $minusHeight = scalar @{$uniqPosReadsMinus{$chr}{$coord}};

					if($ots){
						if(defined($uniqPosReadsPlus{$chr}{$coord})){
							$minusHeight = scalar @{$uniqPosReadsPlus{$chr}{$coord}} + $minusHeight;
						}
					}
										
					if($minusHeight > $minPeak){

						my $support = $plusHeight+$minusHeight; 

						## Store all coordinates with sufficient support
						if($support >= $minSupport){
							my $minPos = $plusPos;
							my $maxPos = $coord;
							my $peakDirection = "Facing";

							## Removed matched peaks
							delete $uniqPosReadsPlusAbMin{$chr}{$plusPos};
							delete $uniqPosReadsMinusAbMin{$chr}{$coord};

							if($coord < ($plusPos - 10)){
								$minPos = $coord;
								$maxPos = $plusPos;
								$peakDirection = "Opposite";
							}

							my $windowPos = int(($minPos+($maxPos-$minPos)/2)/100);

							my $peakKey = "$chr:$windowPos";
							my $support =  $plusHeight+$minusHeight;
							
							if(!defined($peakSupport{$peakKey})){
								$peaks{$peakKey} = [$chr, $minPos, $maxPos, $plusHeight, $minusHeight, $peakDirection];
								$peakSupport{$peakKey} = $support;
							}
							else{
								my $currentSupport = $peakSupport{$peakKey};
								if($support>$currentSupport){
									$peaks{$peakKey} = [$chr, $minPos, $maxPos, $plusHeight, $minusHeight, $peakDirection];
									$peakSupport{$peakKey} = $support;
								}
							}
						}
					}
				}
			}
		}
	}
}

my $totalNrpeaks = 0;

open(OUTFILE,"> $outfile") or die "Can't open file: $outfile\n";

## Print header information
resultsHeader();

## Loop through clipped peaks and print output...
foreach my $key (sort keys %peaks) {
	my @value = @{$peaks{$key}};
	my $chr = $value[0];
	my $minPos = $value[1];
	my $maxPos = $value[2];
	my $plusHeight = $value[3];
	my $minusHeight = $value[4];
	my $peakDirection = $value[5];

	if(!($chr =~ /^(chr.+|hpv.+)/)){
		$printchr = "chr$chr";
	}  else {
		$printchr = $chr;
	}
	
	if($ots){
		print OUTFILE "$printchr\t$minPos\t$minPos\t$plusHeight\n";
	}
	else{
		print OUTFILE "$printchr\t$minPos\t$maxPos\t$plusHeight\t$minusHeight\t$peakDirection\n";
	}
	
	$totalNrpeaks++;
	
	if ($printReads) {
		
		my $fileName = $chr."_".$minPos."_".$maxPos;
		
		if (-e "$fileName.fwd.fasta") {
			open(FWD,">$fileName.fwd.fasta") || die "Couldn't open file $fileName.fwd.fasta, $!";
		} else {
			open(FWD,">>$fileName.fwd.fasta") || die "Couldn't open file $fileName.fwd.fasta, $!";
		}
		if (-e "$fileName.rev.fasta") {
			open(REV,">$fileName.rev.fasta") || die "Couldn't open file $fileName.rev.fasta $!";
		} else {
			open(REV,">>$fileName.rev.fasta") || die "Couldn't open file $fileName.rev.fasta, $!";
		}
		
		for (my $i = $minPos; $i <= $maxPos; $i++) {
			
			foreach (@{$uniqPosReadsPlus{$chr}{$i}}) {
				print FWD ">$_\n", @{$validReads{$_}}[5], "\n";
			}
			
			foreach (@{$uniqPosReadsMinus{$chr}{$i}}) {
				print REV ">$_\n", @{$validReads{$_}}[5], "\n";
			}
		}
		close (FWD);
		close (REV);
	}
	
}

if($ots){
	$singlePeaks = 1;
}

if ($singlePeaks) {

	if(!$ots){
		print OUTFILE "#Single peaks on Plus strand:\n";
	}
		
	foreach my $chr (sort keys %uniqPosReadsPlusAbMin) {
			
		foreach my $plusPos (sort {$a<=>$b} keys %{$uniqPosReadsPlusAbMin{$chr}}){
			
			my $plusHeight = scalar @{$uniqPosReadsPlusAbMin{$chr}{$plusPos}};

			if(!($chr =~ /^(chr.+|hpv.+)/)){
				$printchr = "chr$chr";
			}  else {
				$printchr = $chr;
			}
			
			print OUTFILE "$printchr\t$plusPos\t$plusPos\t$plusHeight\n";
			$totalNrpeaks++;
			
			if ($printReads) {
				
				my $fileName = $chr."_".$plusPos."_SinglePeak";
				
				if (-e "$fileName.fwd.fasta") {
					open(FWD,">$fileName.fwd.fasta") || die "Couldn't open file $fileName.fwd.fasta, $!";
				} else {
					open(FWD,">>$fileName.fwd.fasta") || die "Couldn't open file $fileName.fwd.fasta, $!";
				}
				for (my $i = $plusPos-10; $i <= $plusPos+10; $i++) {
					
					foreach (@{$uniqPosReadsPlus{$chr}{$i}}) {
						print FWD ">$_\n", @{$validReads{$_}}[5], "\n";
					}
				}
				close (FWD);
			}
		}
	}

	if(!$ots){
		print OUTFILE "#Single peaks on Minus strand:\n";
	}
	
	foreach my $chr (sort keys %uniqPosReadsMinusAbMin) {

		foreach my $minusPos (sort {$a<=>$b} keys %{$uniqPosReadsMinusAbMin{$chr}}){
			
			my $minusHeight = scalar @{$uniqPosReadsMinusAbMin{$chr}{$minusPos}};

			if(!($chr =~ /^(chr.+|hpv.+)/)){
				$printchr = "chr$chr";
			} else {
				$printchr = $chr;
			}
			
			print OUTFILE "$printchr\t$minusPos\t$minusPos\t$minusHeight\n";
			$totalNrpeaks++;
			
			if ($printReads) {
				
				my $fileName = $chr."_".$minusPos."_SinglePeak";
				
				if (-e "$fileName.rev.fasta") {
					open(REV,">$fileName.rev.fasta") || die "Couldn't open file $fileName.rev.fasta $!";
				} else {
					open(REV,">>$fileName.rev.fasta") || die "Couldn't open file $fileName.rev.fasta, $!";
				}
				for (my $i = $minusPos-10; $i <= $minusPos+10; $i++) {
					
					foreach (@{$uniqPosReadsMinus{$chr}{$i}}) {
						print REV ">$_\n", @{$validReads{$_}}[5], "\n";
					}
				}
				close (REV);
			}
		}
	}
}

close(OUTFILE);

unless($silent){
	print STDOUT " - All done! $totalNrpeaks sites written to $outfile\n\n";
}

########################
##
## End of main program
##
########################

## Determine end position
sub endPosition {
	my $endPos = 0;
	if (!(scalar(@_) == 2)) {
		print "error: wrong number of arguments for endPosition"
		} else {
			$endPos = $_[1];
			$_[0] =~ s/(\d+)[NMD=X]/$endPos+=$1/eg;
		}
		return int($endPos);
	}

## Determine read strand
sub readStrand {
	my $strand;
	if (!(scalar(@_) == 1)) {
		print "error: wrong number of arguments for readStrand"
		} else {
			if ($_[0] & 16) {
				$strand = "-";
				} else {
					$strand = "+";
				}
			}	
			return $strand;
		} 

## Alignmet quality filter
sub alignQV {
	my $missPos;
	my $softclipedBases;
	my @samline2;
	my $readlen;
	my $cigar2;
	my $fraction = 0;

	if (!(scalar(@_) == 1)) {
		print "error: wrong number of arguments for alignQV"
		} else {
			$missPos = $_[0];
			@samline2 = split("\t", $_[0]);
			$cigar2 = $samline2[5];
			$readlen = length($samline2[9]);
			if ($cigar2 =~ m/^(\d+)S/) {
				$softclipedBases = $1;
				} elsif ($cigar2 =~ m/(\d+)S$/) {
					$softclipedBases = $1;
				}
				if ($missPos =~ m/MD\:Z\:[0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*/) {
					$missPos = $&;
					($missPos) =~ s/MD\:Z\://;
					my @missMatches = $missPos =~ m/[A-Z]|\^[A-Z]+/g;
					my $missMatchesNr = @missMatches;
					$fraction = ($missMatchesNr/($readlen - $softclipedBases))*100;
			# print "$readlen, $missPos, $cigar2, $softclipedBases, @missMatches, $missMatchesNr, $fraction\n";
			} else {
				$fraction = 0;
			}
		}
		return $fraction;
	}

## Header
sub resultsHeader {
	my $cmdline = qx(ps -o args $$);
	$cmdline =~ s/\n//g;
	$cmdline =~ s/COM.*insider.pl/insider.pl/;
	my $analysisDate = localtime();
	
	print OUTFILE "##\n";
	print OUTFILE "## Results were generated using Insider $insiderVersion\n";
	print OUTFILE "## The following command line was executed\:\n";
	print OUTFILE "## $cmdline\n";
	print OUTFILE "## Analysis was performed on $analysisDate\n";
	print OUTFILE "## Thank you for choosing Insider!\n";
	print OUTFILE "##\n";
}

## Argument and error handling
sub get_args_and_error_check {
	
	if (@ARGV == 0) {pod2usage(-exitval => 2, -verbose => 0);}
	
	my ($filename, $outfile, $offset, $minPeak, $minSupport, $minClipLen, $minMappingQV, $maxMappingMM, $minAlignLen, $silent,  $singlePeaks, $printReads, $ots);

	my $result = GetOptions("--help"           => sub{local *_=\$_[1];
		pod2usage(-exitval =>2, -verbose => 1)},
		"-i=s"             =>\$filename,
		"-o=s"             =>\$outfile,
		"-offset=i"        =>\$offset,
		"-minp=i"          =>\$minPeak,
		"-mins=i"          =>\$minSupport,
		"-minc=i"          =>\$minClipLen,
		"-minqv=i"         =>\$minMappingQV,
		"-maxmm=i"         =>\$maxMappingMM,
		"-minl=i"          =>\$minAlignLen,
		"-pr!"			   =>\$printReads,
		"-ps!"			   =>\$singlePeaks,
		"-ots!"            =>\$ots,					
		"-silent!"         =>\$silent)|| pod2usage(-exitval => 2, -verbose => 1);

	my $error_to_print;

	if(defined($offset)) {
		if($offset < 1){
			$error_to_print .= "\tInvalid value of 'offset' $offset\n"
		}
	}
	else{
		$offset=10; 	# Set default value for 'o'
	} 
	
	if(defined($minPeak)) {
		if($minPeak < 1){
			$error_to_print .= "\tInvalid value of 'minp' $minPeak\n"
		}
	}
	else{
		$minPeak=3;  # Set default value for 'minp'
	}

	if(defined($minSupport)) {
		if($minSupport < 1){
			$error_to_print .= "\tInvalid value of 'mins' $minSupport\n"
		}
	}
	else{
		$minSupport=10;  # Set default value for 'mins'
	}

	if(defined($minClipLen)) {
		if($minClipLen < 1){
			$error_to_print .= "\tInvalid value of 'minc' $minClipLen\n"
		}
	}
	else{
		$minClipLen=20;  # Set default value for 'minc'
	}

	if(defined($minMappingQV)) {
		if($minMappingQV < 0){
			$error_to_print .= "\tInvalid value of 'minqv' $minMappingQV\n"
		}
	}
	else{
		$minMappingQV=0;  # Set default value for 'minqv'
	}

	if(defined($maxMappingMM)) {
		if($maxMappingMM > 100){
			$error_to_print .= "\tInvalid value of 'maxmm' $maxMappingMM\n"
		}
	}
	else{
		$maxMappingMM=1;  # Set default value for 'maxmm'
	}

	if(defined($minAlignLen)) {
		if($minAlignLen < 0){
			$error_to_print .= "\tInvalid value of 'minl' $minAlignLen\n"
		}
	}
	else{
		$minAlignLen=1;  # Set default value for 'minl'
	}

	unless(defined($filename)) {
		$error_to_print .= "\tNo input file specified.\n"
	}

	unless(defined($outfile)) {
		$error_to_print .= "\tNo outfile prefix specified.\n"
	}
	
	if(defined $error_to_print) {
		my $error_msg="ERROR(s):\n$error_to_print\n";
		pod2usage(-message => $error_msg, -exitval => 2, -verbose => 0);
	}

	else{
		return ($filename, $outfile, $offset, $minPeak, $minSupport, $minClipLen, $minMappingQV, $maxMappingMM, $minAlignLen, $silent, $singlePeaks, $printReads, $ots);
	}
	
}

__END__

=head1 NAME
	
insider.pl 
	
=head1 SYNOPSIS
	
./insider.pl [options] B<--help> B<-i> B<-o> B<-ots> B<-offset> B<-minp> B<-mins> B<-minc> B<-minqv> B<-maxmm> B<-pr> B<-ps> B<-silent>

=head1 OPTIONS

=over 8

=item [REQUIRED]

=item B<-i>

Input file (aligned reads in SAM format).

=item B<-o>

Output file name.

=item [OPTIONAL]

=item B<-ots>

Set this flag for analysis of CRISPR-Cas9 cleavage sites in off-target sequencing (OTS) data.

=item B<-offset>

Offset between left-end and right-end soft clipped alignments in host genome, used when searching for insertion sites (default 10).

=item B<-minp>

Minimum number of soft clipped reads starting at the exact same position in the insertion site peak. This number of reads is required for both the left-end and right-end alignments. (default=3).

=item B<-mins>

Minimum number soft clipped reads starting at the exact same position in the insertion site peak, when taking both left-end and right-end alignments into account (default=10).

=item B<-minc>

Minimum number of soft clipped bases required for a read to be considered in the Insider analysis (default=20).

=item B<-minqv>

Minimum mapping QV value a read to be considered in the Insider analysis (default=0).

=item B<-maxmm>

Maximum fraction (%) of mismatches in aligned part of a read (default=1).

=item B<-minl>

Minimum alignment length of the reads used in the Insider analysis (default=0).

=item B<-pr>

Generate read files in FASTA format for each identified integration site or single peak.

=item B<-ps>

In the result file include peaks that passed filtering, but did not have a corresponding peak on the opposite strand.

=item B<-silent>

Do not print status to stdout.

=back

=head1 DESCRIPTION

B<This program> will read the given input file and outputs putative integration sites of DNA into a host genome.

=cut
