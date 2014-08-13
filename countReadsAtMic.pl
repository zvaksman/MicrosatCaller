#!/usr/bin/perl -w
#!/usr/local/bin/perl -w
#!/bin/perl -w

######################################
# Author: Honngseok Tae
# Date: 2010
#
# countReadsAllMic.pl
# count reads on microsatellites
######################################

use strict;
use Getopt::Long qw(:config no_ignore_case);


my ($cmd, $helpFlag, $micFn, $samFn, $refFn, $flankingBases, $outFn);
$cmd = "$0 @ARGV";

$flankingBases = 5;

GetOptions(
	"h|?|help"		=> \$helpFlag,
	"mic=s"	=> \$micFn,
	"ref=s"	=> \$refFn,
	"sam=s"	=> \$samFn,
	"flanking=i"	=> \$flankingBases,
	"output=s"	=> \$outFn
) || help(1);


help(0) if defined $helpFlag;

if(!defined $samFn){ 
	$samFn = shift; 
	if($samFn && !defined $micFn && $samFn !~ /.sam/ && $samFn !~ /.bam/){
		$micFn = $samFn;
		$samFn = shift; 
	}
}
if(!defined $samFn){ print STDERR "\nNeed a sam file!!\n\n"; help(1); }
if(!defined $micFn){
	$micFn = shift;
}
if(!defined $micFn){ print STDERR "\nNeed a microsatellite list file!!\n\n"; help(1); }
if(!defined $refFn){ print STDERR "\nNeed a reference sequence file!!\n\n"; help(1); }

if(defined $outFn){
	if($micFn eq $outFn || $samFn eq $outFn  || $refFn eq $outFn ){ print STDERR " Error) Input and output files are same \n"; exit; }
}


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0  -r <ref. sequence file>  -m <microsat file>  <sam|bam file> [-o <out file>] \n";
	print STDERR "  ex) $0  -r ref.fa  -m mic.lst  GM13869_NoCot1.hg19.bam   -o GM13869_NoCot1.inMic.txt\n\n";
	exit($return);
}

my ($in, $out, @arr, $ref, $i, $j, $max, $max_i, %targetInRef, $totMicLen, $preEnd);
$totMicLen = 0;


$in = openInput($micFn);
$preEnd = 0;
while(<$in>){
	s/[\r\n]+//g;
	@arr = split /\t/;
	next if($#arr < 3 );

	$ref = $arr[0];
		
	$i = $#{$targetInRef{$ref}} + 1;

	$targetInRef{$ref}[$i]{start} = $arr[1];
	$targetInRef{$ref}[$i]{end} = $arr[1]+$arr[2]-1;
	$targetInRef{$ref}[$i]{len} = $arr[2];	
	$targetInRef{$ref}[$i]{seq} = $arr[3];
}
close($in);


$out = openOutput($outFn);
print $out "# $cmd\n";

my($mI, $mEnd, @arrTarget);

$ref = "";
my ($tStart, $cigar, $prevMI, $tPos, $totGap, $micGap, $op, $mic, @alleles, $rBuf, $rin, $refSeq);

$in = openInput($samFn);

my $prevTStart = 0;
while(<$in>){
	if(/^@/){
		next;
	}

	#####
	## 0: ID, 1: flag, 2: reference, 3: position, 4: mapping score, 5: cigar, 6/7/8: mate info, 9: seq, 10: quality, 11..: tag.
	my @arr = split /\t/;
	next if($#arr < 2);
		
	$tStart = $arr[3];
	$cigar = $arr[5];

	next if(!defined $arr[2] ||  !defined $targetInRef{$arr[2]} || ($arr[1]&0x4) != 0) ;

	if(defined $ref && $ref eq $arr[2] && $tStart < $prevTStart){
		print STDERR "[Error] Reads in the '$samFn' file is out of order. They should be sorted.\n  previous start : $prevTStart\n  $_\n";
		exit;
	}

	if(!defined $ref || $ref ne $arr[2]){		
		readRefSeq($arr[2]);
		if($ref ne ""){
			while($mI < $mEnd){
				printAlleles(\%{$arrTarget[$mI]}) if($arrTarget[$mI]{alleles});
				delete $arrTarget[$mI]{alleles};
				$mI++;
			}
		}
		$mI = 0;
		@arrTarget = sort {$a->{start} <=> $b->{start}} @{ $targetInRef{$arr[2]} };		
		$mEnd = scalar @arrTarget;		
		$ref = $arr[2];;
	}
	

	my $tEnd = $tStart -1;
	while($cigar =~ /(\d+)(\w)/g){	
		if($2 eq "D" || $2 eq "N" || $2 eq "M"){			
			$tEnd += $1;			
		}
	}
	
	while($mI < $mEnd && $arrTarget[$mI]{start} < $tStart + $flankingBases){
		printAlleles(\%{$arrTarget[$mI]}) if($arrTarget[$mI]{alleles});
		delete $arrTarget[$mI]{alleles};
		$mI++;
	}
	next if($mI == $mEnd);

	my $tI = $mI;
	
	$mic = $arrTarget[$tI];
		
	my ($qPos, $leftGap, $rightGap, $micQStart, $micQEnd, $micAlleleLen, $blockTEnd, $blockLen, $leftClipped);

	$tPos = $tStart;
	$totGap = 0;
	$micGap = 0;
	$qPos = 1;
	$leftGap = 0;
	$rightGap = 0;
	

	$micQStart = ($mic->{start} - $tStart) + 1;

	while($cigar =~ /(\d+)(\w)/g){	
		($blockLen, $op) = ($1, $2);

		if($op eq "S"){ 
			$qPos += $blockLen;
			$micQStart += $blockLen;
			$totGap += $blockLen;
		}
		elsif($op eq "D" || $op eq "N"){			
			$blockTEnd = $tPos + $blockLen - 1; # : the end position of the deletion block.
			if($blockTEnd < $mic->{start}){
				$micQStart -= $blockLen;
			}			
			elsif($tPos <= $mic->{start} && $mic->{start} <= $blockTEnd){
				$micQStart = $qPos;
			}

			
			if($mic->{start} <= $tPos){
				if($blockTEnd <= $mic->{end}) {
					$micGap -= $blockLen;
				}
				elsif($tPos <= $mic->{end}+1 && $mic->{end} < $blockTEnd) {
					$rightGap = $blockTEnd - $mic->{end};
					$micGap -= $mic->{end} - $tPos + 1;
				}            
			}
			elsif($mic->{start} <= $blockTEnd+1 && $blockTEnd <= $mic->{end}){ # $tPos < $mic->{start}
				$leftGap = $mic->{start} - $tPos;
				$micGap -= $blockTEnd - $mic->{start} + 1;            
			}
			elsif($mic->{end} < $blockTEnd){  # tPos < $mic->{start}
				$micGap -= $blockLen;
			}		
			
			$tPos += $blockLen;
			$totGap -= $blockLen;
		}
		elsif($op eq "I"){
			if($mic->{start} <= $tPos && $mic->{end} +1 >= $tPos){
				$micGap += $blockLen;
			}
			elsif($tPos < $mic->{start}){
				$micQStart += $blockLen;
			}
			$qPos += $blockLen;
			$totGap += $blockLen;
		}
		elsif($op eq "M"){
			$blockTEnd = $tPos+$blockLen-1;
			while($tI < $mEnd  && $mic->{end} < $blockTEnd){  
			
				$micQEnd = $qPos + ($mic->{end} - $tPos) + $rightGap;
				$micAlleleLen = $mic->{len} + $micGap;         
				
				if($micQEnd > $micQStart){  # if micQEnd is smaller than micQStart, it means there is a deletion including the whole microsatellite sequence.
					if($micQEnd - $micQStart + 1 == $micAlleleLen){
						addAlleleToMic($mic, \@arr, $micAlleleLen, $micQStart,  $micQEnd);
						#addAlleleToMic(mic, &listAlleles[mic->id], $micAlleleLen, mapQual, $micQStart,  $micQEnd, seqRaw, seqLen,
						#		$leftGap, $rightGap); #TODD ...baseQual,
					}
					else{
						print STRERR sprintf("[Warning] %s : (mic start %d, mic end %d in query) is not same with allele length : %d\nMic : %d-%d\n",
								$arr[0], $micQStart, $micQEnd, $micAlleleLen, $mic->{start}, $mic->{end});
						print STRERR sprintf("$tStart %d, micGap %d, mic qEnd :  qPos(%d) + (mic->end(%d) - tPos(%d)) + rightGap(%d)\n",
								$tStart, $micGap, $qPos, $mic->{end}, $tPos, $rightGap); close($out) if($out); exit(1);
					}
				}

				$tI++;
				if($tI < $mEnd){
					$mic = $arrTarget[$tI];
					$micQStart = ($mic->{start} - $tStart) + $totGap + 1;
				}
				
				$micGap = 0;
				$leftGap = 0; 
				$rightGap = 0;					
			}
			
			last if($tI >= $mEnd ||  $mic->{end} > $tEnd - $flankingBases);
			$tPos += $blockLen;
			$qPos += $blockLen;
		}
	}


	$prevTStart = $tStart;
}

if($ref ne ""){
	while($mI < $mEnd){
		printAlleles(\%{$arrTarget[$mI]}) if($arrTarget[$mI]{alleles});
		delete $arrTarget[$mI]{alleles};
		$mI++;
	}
}


close($in);


sub printAlleles
{
	my ($mic) = @_;
	my @alleles = sort {$mic->{alleles}{$b}{num} <=> $mic->{alleles}{$a}{num}} keys %{$mic->{alleles}};

	if($mic->{alleles}{$alleles[0]}{num} > 1){
		print $out "$ref\t$mic->{start}\t$mic->{len}\t" . 
			substr($refSeq, $mic->{start}-1-$flankingBases, $flankingBases) ."-".
			"$mic->{seq}-". substr($refSeq, $mic->{end},$flankingBases);
		foreach my $allele (@alleles) {
			my @lefts = sort {$mic->{alleles}{$allele}{left}{$b} <=> $mic->{alleles}{$allele}{left}{$a}} keys %{$mic->{alleles}{$allele}{left}};
			my @rights = sort {$mic->{alleles}{$allele}{right}{$b} <=> $mic->{alleles}{$allele}{right}{$a}} keys %{$mic->{alleles}{$allele}{right}};

			print $out "\t" . length($allele). "\t". ($mic->{alleles}{$allele}{num}) .
				"\t$lefts[0]-$allele-$rights[0]" 
				if($mic->{alleles}{$allele}{num} > 1);
		}
		print $out "\n";
	}
}

sub addAlleleToMic
{
	my ($mic, $arr, $micAlleleLen, $micQStart,  $micQEnd) = @_;
	
	my ($start, $len);

	my $readLen = length($arr->[9]);
	my $micSeq = substr($arr->[9], $micQStart-1, $micAlleleLen);
	$mic->{alleles}{$micSeq}{num}++;
	#if(!$mic->{alleles}{$micSeq}{left} || length($mic->{alleles}{$micSeq}{left}) < $flankingBases){
		($start, $len) = ($micQStart-1-$flankingBases, $flankingBases);
		($start, $len) = (0, $flankingBases + $start)
			if($start < 0);
		if($len == $flankingBases) { $mic->{alleles}{$micSeq}{left}{substr($arr->[9], $start, $len)}++; }
		else { $mic->{alleles}{$micSeq}{left}{substr($arr->[9], $start, $len)} = 0; }
	#}

	#if(!$mic->{alleles}{$micSeq}{right} || length($mic->{alleles}{$micSeq}{right}) < $flankingBases){
		($start, $len) = ($micQEnd, $flankingBases);
		($start, $len) = ($micQEnd, $readLen-$micQEnd)
			if($readLen -  $micQEnd < $flankingBases);
		if($len == $flankingBases) { $mic->{alleles}{$micSeq}{right}{substr($arr->[9], $start, $len)}++; }
		else{ $mic->{alleles}{$micSeq}{right}{substr($arr->[9], $start, $len)} = 0; }
	#}
}


sub readRefSeq{
	my $ref = shift;

	my ($name);
	print STDERR "reading '$ref' sequence...\n";
	$rin = openInput($refFn) if(!$rin); 
	$rBuf = <$rin> if(!$rBuf); 
	while(defined $rBuf){
		if($rBuf =~ /^>/){			
			$name = $';
			$name =~ s/^\s*//g; $name =~ s/\s.*//g;
			if($name eq $ref){
				$rBuf = <$rin>;
				last;
			}
		}
		$rBuf = <$rin>;
	}
	
	$refSeq = "";
	
	while(defined $rBuf){
		last if($rBuf =~ /^>/);
		$rBuf =~ s/[\r\n]+//g;
		$refSeq .= uc $rBuf;
		$rBuf = <$rin>;
	}

	if($refSeq eq ""){
		print STDERR "[Error] Could not read sequences for '$ref'. The order of references in '$refFn' and '$samFn' may be differ.\n";
		exit(1);
	}
}

sub openInput
{
	my ($fn) = @_;

	print "reading $fn\n";

	my ($fd);
	if($fn =~ /.bam$/){
		open($fd, "samtools view -h $fn |");
	}
	else{
		open($fd, $fn =~ /\.gz/ ? "zcat $fn|" : ($fn =~ /\.bz2/ ? "bunzip2 -c $fn|" : $fn)) || die "Could not open '$fn' : $!\n";
	}
	return $fd;
}

sub openOutput
{
	my ($fn) = @_;

	return *STDOUT unless defined $fn;

	my ($fd);
	open($fd, $fn =~ /.gz$/ ? "| gzip -c > $fn" : ($fn =~ /\.bz2/ ? "| bzip2 -c > $fn" : ">$fn")) || die "Could not write '$fn' : $!\n";
	return $fd;
}


sub rc
{
	my ($seq, $type) = @_;

	my $rc = reverse $seq;
	if(defined $type && ($type eq "rna" || $type eq "RNA")) # RNA
	{   $rc =~ y/acgtuACGTU/ugcaaUGCAA/;  }
	else ## DNA
	{   $rc =~ y/acgtuACGTU/tgcaaTGCAA/;  }

	return $rc;
}


