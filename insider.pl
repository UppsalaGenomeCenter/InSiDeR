#!/usr/bin/perl -w

#############################################################################
##
##   InSiDeR - virus integration site detection tool							 
##
##   Copyright (C) 2014 Ignas Bunikis, Adam Ameur							 
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
#############################################################################

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my $insiderVersion = "version 1.0";

# Read command line arguments
my ($filename, $outfile, $offset, $minPeak, $minSupport, $minClipLen, $silent) = get_args_and_error_check();

unless($silent){
	print STDOUT "\n************ Running InSiDeR $insiderVersion **********\n\n";
}

my %validReads = ();
my $allReadsNr = 0;
my $validReadsNr = 0;
my %uniqPosReadsPlus = ();
my %uniqPosReadsMinus = ();

unless($silent){
	print STDOUT "Reading data from $filename...\n";
}

open(FILE, "<", $filename) or die "cannot open: $!";

while (my $row = <FILE>) {
	chomp $row;

	if(!($row =~ /^@/)){
		my @samline = split("\t", $row);
		my $cigar = $samline[5];
		$allReadsNr++;
				
		if (!($samline[1] & 4)) {
			if ((($cigar =~ m/^(\d+)S/) && ($1 >= $minClipLen)) || (($cigar =~ m/(\d+)S$/) && ($1 >= $minClipLen))) {
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
				
				if($strand eq "+"){
					push @{$uniqPosReadsPlus{$chr}{$end_pos}}, $samline[0];
				}
				if($strand eq "-"){
					push @{$uniqPosReadsMinus{$chr}{$start_pos}}, $samline[0];
				}
				
			}
		}
	}
}

close(FILE);

unless($silent){
	print STDOUT " - $validReadsNr/$allReadsNr reads met criteria for further analysis\n\n";
	print STDOUT "Searching for sites of integrated DNA...\n"
}

my %peaks = ();
my %peakSupport = ();

## Loop through clipped peaks on plus strand...
foreach my $chr (sort keys %uniqPosReadsPlus) {
		
	foreach my $plusPos (sort {$a<=>$b} keys %{$uniqPosReadsPlus{$chr}}){

		my $plusHeight = scalar @{$uniqPosReadsPlus{$chr}{$plusPos}};

		if($plusHeight > $minPeak){

			my $bestSupport = 0;
			my $bestCoord = "NA";
			my $bestEnd = "NA";
			
			## Check if there is some peak on minus strand in a window surrounding plus strand peak..
			for(my $coord = $plusPos-$offset; $coord<=$plusPos; $coord++){

				if(defined($uniqPosReadsMinus{$chr}{$coord})){
					my $minusHeight = scalar @{$uniqPosReadsMinus{$chr}{$coord}};

					if($minusHeight > $minPeak){
						my $support = $plusHeight+$minusHeight; 

						## Store all coordinates with sufficient support
						if($support >= $minSupport){
							my $minPos = $plusPos;
							my $maxPos = $coord;
							if($coord < $plusPos){
								$minPos = $coord;
								$maxPos = $plusPos;
							}

							my $windowPos = int(($minPos+($maxPos-$minPos)/2)/100);

							my $peakKey = "$chr:$windowPos";
							my $support =  $plusHeight+$minusHeight;
							
							if(!defined($peakSupport{$peakKey})){
								$peaks{$peakKey} = [$chr, $minPos, $maxPos, $plusHeight, $minusHeight];
								$peakSupport{$peakKey} = $support;
							}
							else{
								my $currentSupport = $peakSupport{$peakKey};
								if($support>$currentSupport){
									$peaks{$peakKey} = [$chr, $minPos, $maxPos, $plusHeight, $minusHeight];
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

## Loop through clipped peaks and print output...
foreach my $key (sort keys %peaks) {
	my @value = @{$peaks{$key}};
	my $chr = $value[0];
	my $minPos = $value[1];
	my $maxPos = $value[2];
	my $plusHeight = $value[3];
	my $minusHeight = $value[4];

	print OUTFILE "$chr\t$minPos\t$maxPos\t$plusHeight\t$minusHeight\n";
	$totalNrpeaks++;
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


sub endPosition {
	my $endPos = 0;
	if (!(scalar(@_) == 2)) {
		print "error: wrong number of arguments for endPosition"
	} else {
	$endPos = $_[1];
	$_[0] =~ s/(\d+)[NMD]/$endPos+=$1/eg;
	}
	return int($endPos);
}

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

# Argument and error handling
sub get_args_and_error_check{
	
	if (@ARGV == 0) {pod2usage(-exitval => 2, -verbose => 0);}
	
	my ($filename, $outfile, $offset, $minPeak, $minSupport, $minClipLen, $silent);

	my $result = GetOptions("--help"           => sub{local *_=\$_[1];
							                      pod2usage(-exitval =>2, -verbose => 1)},
		                    "-i=s"             =>\$filename,
							"-o=s"             =>\$outfile,
							"-offset=i"        =>\$offset,
							"-minp=i"          =>\$minPeak,
							"-mins=i"          =>\$minSupport,
	                        "-minc=i"          =>\$minClipLen,
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
	
	unless(defined($filename)) {
		$error_to_print .= "\tNo input file specified.\n"
	}

    unless(defined($outfile)) {
		$error_to_print .= "\tNo outfile specified.\n"
	}
	
	if(defined $error_to_print) {
		my $error_msg="ERROR(s):\n$error_to_print\n";
		pod2usage(-message => $error_msg, -exitval => 2, -verbose => 0);
	}

	else{
		return ($filename, $outfile, $offset, $minPeak, $minSupport, $minClipLen, $silent);
	}
	
}

__END__

=head1 NAME
	
splitseek.pl 
	
=head1 SYNOPSIS
	
./insider.pl [options] B<--help> B<-i> B<-o> B<-offset> B<-minp> B<-mins> B<-minc> B<--silent>

=head1 OPTIONS

=over 8

=item [REQUIRED]

=item B<-i>

InSiDeR input file (SAM format).

=item B<-o>

Output file for InSiDeR results.

=item [OPTIONAL]

=item B<-offset>

Offest between left-end and right-end soft clipped alignments in host genome, used when searching for insertion sites (default 10).

=item B<-minp>

Minimum number of soft clipped reads starting at the exact same position in the insertion site peak. This number of reads is required for both the left-end and right-end alignments. (default=3).

=item B<-mins>

Minimum number soft clipped reads starting at the exact same position in the insertion site peak, when taking both left-end and right-end alignments into account (default=10).

=item B<-minc>

Minimum number of soft clipped bases required for a read to be considered in the InSiDeR analysis (default=20).

=item B<--silent>
Do not print status to stdout.

=back

=head1 DESCRIPTION

B<This program> will read the given input file and outputs putative integration sites of DNA into a host genome.

=cut
