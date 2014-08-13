#!/usr/bin/perl 
#!/usr/local/bin/perl
#!/bin/perl 

######################################
# Author: Zalman Vaksman
# Date: 2013
#
# MSTanalysis.pl
######################################

#use strict;
#use warnings "all";
use POSIX;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case);
use Thread;
use Cwd;
use feature 'state';


my $MinDP = 10;
my $AlDP = 3;

GetOptions(
	"h|?|help"	=> \$helpFlag,

	"out=s"	=> \$outHeader,
	"DP=i"	=> \$MinDP,
	"A=i"	=> \$AlDP,
	
	
) || help(1);


sub help{
	my $return = shift;
	$return = 0 if(!defined $return);
	print STDERR " Usage: $0 -o <header of output>  [optional -DP <Minimum Total Depth> -A <minimum depth at allele>]\n\n";
	print STDERR "  ex 1) $0 -o myHeader -DP number Minimum total depth\n\n";
	print STDERR "  ex 2) $0 -o myHeader -DP 10 -A 3 \n\n";
	print STDERR "  
		------------------------------------
		If this is run as a stand-alone script, use these options
		
		Options  \n				
     -o 			 - Output and input file. The input
     					will be [-o].GATK.indel
     -DP [number]    - Cutoff per total reads 
                       counted (that passed -A (see below))	
                       Default set at 10\n
     -A [number]     - Cutoff per allele (this effects 
                       total dpeth (-DP)) 
                       Default set at 3\n
				
				\n\n";
	exit($return);
}

if(!defined $outHeader){ print STDERR "\nNeed a header of output files!!\n\n"; help(1); }
#if(!defined $refFn){ print STDERR "\nNeed a reference sequence file!!\n\n"; help(1); }
#if(! -e $outHeader){ print STDERR "\n[Error] Could not find the file '$outHeader.GATK.indel' \n"; exit(1); }



############################# data optimization and analysis ##########################################



##################### Input files

my $file = "$outHeader.GATK.indel";
my $file1 = 'pubsMarkerBand.txt';
open(my $fh, '<', $file) or die "Can't read file '$file' [$!]\n";
open(my $fh1, '<', $file1) or die "Can't read file '$file1' [$!]\n";


##################### Output files

@header = qw (Chr ChrMarker Position RefLength MSTunitL Period RefMicrSat UnqueReads GenotypeLeng Genotype TotalDP Read1Length Read1DP Read1%Tot Read1Seq Read2Length Read2SizeDiff Read2DP Read2%Tot Read2Seq Read3Length Read3SizeDiff Read3DP Read3%Tot Read3Seq Read4Length Read4SizeDiff Read4DP Read4%Tot Read4Seq Read5Length Read5SizeDiff Read5DP Read5%Tot Read5Seq Read6Length Read6SizeDiff Read6DP Read6%Tot Read6Seq Read7Length Read7SizeDiff Read7DP Read7%Tot Read7Seq Read8Length Read8SizeDiff Read8DP Read8%Tot Read8Seq Read9Length Read9SizeDiff Read9DP Read9%Tot Read9Seq Read10Length Read10SizeDiff Read10DP Read10%Tot Read10Seq Read11Length Read11SizeDiff Read11DP Read11%Tot Read11Seq Read12Length Read12SizeDiff Read12DP Read12%Tot Read12Seq Read13Length Read13SizeDiff Read13DP Read13%Tot Read13Seq Read14Length Read14SizeDiff Read14DP Read14%Tot Read14Seq Read15Length Read15SizeDiff Read15DP Read15%Tot Read15Seq Read16Length Read16SizeDiff Read16DP Read16%Tot Read16Seq Read17Length Read17SizeDiff Read17DP Read17%Tot Read17Seq Read18Length Read18SizeDiff Read18DP Read18%Tot Read18 Read19Length Read19SizeDiff Read19DP Read19%Tot Read19SeqSeq Read20Length Read20SizeDiff Read20DP Read20%Tot Read20Seq Read21Length Read21SizeDiff Read21DP Read21%Tot Read21Seq Read22Length Read22SizeDiff Read22DP Read22%Tot Read22Seq Read23Length Read23SizeDiff Read23DP Read23%Tot Read23Seq Read24Length Read24SizeDiff Read24DP Read24%Tot Read24Seq Read25Length Read25SizeDiff Read25DP Read25%Tot Read25Seq); 
$header = join ("\t", @header);
my $output1 = "$file.DP$MinDP.A$AlDP.result";
open (OUT1, ">>$output1");
#print OUT1 "Chr#, ChrMarker, Position, RefLength, MSTunitL, Period, RefMicrSat, UnqueReads, Genotype(length), Genotype, TotalDP,Read1Length, Read1DP, Read1%Tot, Read1Seq, Read2Length, Read2SizeDiff, Read2DP, Read2%Tot, Read2Seq, Read3Length, Read3SizeDiff, Read3DP, Read3%Tot, Read3Seq, Read4Length, Read4SizeDiff, Read4DP, Read4%Tot, Read4Seq, Read5Length, Read5SizeDiff, Read5DP, Read5%Tot, Read5Seq, Read6Length, Read6SizeDiff, Read6DP, Read6%Tot, Read6Seq, Read7Length, Read7SizeDiff, Read7DP, Read7%Tot, Read7Seq, Read8Length, Read8SizeDiff, Read8DP, Read8%Tot, Read8Seq, Read9Length, Read9SizeDiff, Read9DP, Read9%Tot, Read9Seq, Read10Length, Read10SizeDiff, Read10DP, Read10%Tot, Read10Seq,Read11Length, Read11SizeDiff, Read11DP, Read11%Tot, Read11Seq,Read12Length, Read12SizeDiff, Read12DP, Read12%Tot, Read12Seq,Read13Length, Read13SizeDiff, Read13DP, Read13%Tot, Read13Seq,Read14Length, Read14SizeDiff, Read14DP, Read14%Tot, Read14Seq,Read15Length, Read15SizeDiff, Read15DP, Read15%Tot, Read15Seq,Read16Length, Read16SizeDiff, Read16DP, Read16%Tot, Read16Seq,Read17Length, Read17SizeDiff, Read17DP, Read17%Tot, Read17Seq,Read18Length, Read18SizeDiff, Read18DP, Read18%Tot, Read18,Read19Length, Read19SizeDiff, Read19DP, Read19%Tot, Read19SeqSeq,Read20Length, Read20SizeDiff, Read20DP, Read20%Tot, Read20Seq,Read21Length, Read21SizeDiff, Read21DP, Read21%Tot, Read21Seq,Read22Length, Read22SizeDiff, Read22DP, Read22%Tot, Read22Seq,Read23Length, Read23SizeDiff, Read23DP, Read23%Tot, Read23Seq,Read24Length, Read24SizeDiff, Read24DP, Read24%Tot, Read24Seq,Read25Length, Read25SizeDiff, Read25DP, Read25%Tot, Read25Seq\n"; 
#print OUT1 join "\t", @header;
print OUT1 "$header\n";
if ($MinDP > 0) {
my $output2 = "$file.BelowCutOff.DP$MinDP.A$AlDP.result";
open (OUT2, ">>$output2");

my $output10 = "$file.alldata.DP$MinDP.A$AlDP.result";
open (OUT10, ">>$output10");

#print OUT2 "Chr#, ChrMarker, Position, RefLength, RefMicrSat, UnqueReads, Genotype(length), Genotype, TotalDP, Read1Length, Read1DP, Read1%Tot, Read1Seq, Read2Length, Read2SizeDiff, Read2DP, Read2%Tot, Read2Seq, Read3Length, Read3SizeDiff, Read3DP, Read3%Tot, Read3Seq, Read4Length, Read4SizeDiff, Read4DP, Read4%Tot, Read4Seq, Read5Length, Read5SizeDiff, Read5DP, Read5%Tot, Read5Seq, Read6Length, Read6SizeDiff, Read6DP, Read6%Tot, Read6Seq, Read7Length, Read7SizeDiff, Read7DP, Read7%Tot, Read7Seq, Read8Length, Read8SizeDiff, Read8DP, Read8%Tot, Read8Seq, Read9Length, Read9SizeDiff, Read9DP, Read9%Tot, Read9Seq, Read10Length, Read10SizeDiff, Read10DP, Read10%Tot, Read10Seq,Read11Length, Read11SizeDiff, Read11DP, Read11%Tot, Read11Seq,Read12Length, Read12SizeDiff, Read12DP, Read12%Tot, Read12Seq,Read13Length, Read13SizeDiff, Read13DP, Read13%Tot, Read13Seq,Read14Length, Read14SizeDiff, Read14DP, Read14%Tot, Read14Seq,Read15Length, Read15SizeDiff, Read15DP, Read15%Tot, Read15Seq,Read16Length, Read16SizeDiff, Read16DP, Read16%Tot, Read16Seq,Read17Length, Read17SizeDiff, Read17DP, Read17%Tot, Read17Seq,Read18Length, Read18SizeDiff, Read18DP, Read18%Tot, Read18,Read19Length, Read19SizeDiff, Read19DP, Read19%Tot, Read19SeqSeq,Read20Length, Read20SizeDiff, Read20DP, Read20%Tot, Read20Seq,Read21Length, Read21SizeDiff, Read21DP, Read21%Tot, Read21Seq,Read22Length, Read22SizeDiff, Read22DP, Read22%Tot, Read22Seq,Read23Length, Read23SizeDiff, Read23DP, Read23%Tot, Read23Seq,Read24Length, Read24SizeDiff, Read24DP, Read24%Tot, Read24Seq,Read25Length, Read25SizeDiff, Read25DP, Read25%Tot, Read25Seq\n"; 
#print OUT10 "Chr#, ChrMarker, Position, RefLength, RefMicrSat, UnqueReads, Genotype(length), Genotype, TotalDP, Read1Length, Read1DP, Read1%Tot, Read1Seq, Read2Length, Read2SizeDiff, Read2DP, Read2%Tot, Read2Seq, Read3Length, Read3SizeDiff, Read3DP, Read3%Tot, Read3Seq, Read4Length, Read4SizeDiff, Read4DP, Read4%Tot, Read4Seq, Read5Length, Read5SizeDiff, Read5DP, Read5%Tot, Read5Seq, Read6Length, Read6SizeDiff, Read6DP, Read6%Tot, Read6Seq, Read7Length, Read7SizeDiff, Read7DP, Read7%Tot, Read7Seq, Read8Length, Read8SizeDiff, Read8DP, Read8%Tot, Read8Seq, Read9Length, Read9SizeDiff, Read9DP, Read9%Tot, Read9Seq, Read10Length, Read10SizeDiff, Read10DP, Read10%Tot, Read10Seq,Read11Length, Read11SizeDiff, Read11DP, Read11%Tot, Read11Seq,Read12Length, Read12SizeDiff, Read12DP, Read12%Tot, Read12Seq,Read13Length, Read13SizeDiff, Read13DP, Read13%Tot, Read13Seq,Read14Length, Read14SizeDiff, Read14DP, Read14%Tot, Read14Seq,Read15Length, Read15SizeDiff, Read15DP, Read15%Tot, Read15Seq,Read16Length, Read16SizeDiff, Read16DP, Read16%Tot, Read16Seq,Read17Length, Read17SizeDiff, Read17DP, Read17%Tot, Read17Seq,Read18Length, Read18SizeDiff, Read18DP, Read18%Tot, Read18,Read19Length, Read19SizeDiff, Read19DP, Read19%Tot, Read19SeqSeq,Read20Length, Read20SizeDiff, Read20DP, Read20%Tot, Read20Seq,Read21Length, Read21SizeDiff, Read21DP, Read21%Tot, Read21Seq,Read22Length, Read22SizeDiff, Read22DP, Read22%Tot, Read22Seq,Read23Length, Read23SizeDiff, Read23DP, Read23%Tot, Read23Seq,Read24Length, Read24SizeDiff, Read24DP, Read24%Tot, Read24Seq,Read25Length, Read25SizeDiff, Read25DP, Read25%Tot, Read25Seq\n"; 

#print OUT2 join "\t", @header;
print OUT2 "$header\n";

#print OUT10 join "\t", @header;
print OUT10 "$header\n";
}


my $output3 = "$outHeader.multi.DP$MinDP.A$AlDP.result";
open (OUT3, ">>$output3");



my $output4 = "$outHeader.single.DP$MinDP.A$AlDP.result";
open (OUT4, ">>$output4");

#print OUT3 "Chr#, ChrMarker, Position, RefLength, RefMicrSat, UnqueReads, Genotype(length), Genotype, TotalDP, Read1Length, Read1DP, Read1%Tot, Read1Seq, Read2Length, Read2SizeDiff, Read2DP, Read2%Tot, Read2Seq, Read3Length, Read3SizeDiff, Read3DP, Read3%Tot, Read3Seq, Read4Length, Read4SizeDiff, Read4DP, Read4%Tot, Read4Seq, Read5Length, Read5SizeDiff, Read5DP, Read5%Tot, Read5Seq, Read6Length, Read6SizeDiff, Read6DP, Read6%Tot, Read6Seq, Read7Length, Read7SizeDiff, Read7DP, Read7%Tot, Read7Seq, Read8Length, Read8SizeDiff, Read8DP, Read8%Tot, Read8Seq, Read9Length, Read9SizeDiff, Read9DP, Read9%Tot, Read9Seq, Read10Length, Read10SizeDiff, Read10DP, Read10%Tot, Read10Seq,Read11Length, Read11SizeDiff, Read11DP, Read11%Tot, Read11Seq,Read12Length, Read12SizeDiff, Read12DP, Read12%Tot, Read12Seq,Read13Length, Read13SizeDiff, Read13DP, Read13%Tot, Read13Seq,Read14Length, Read14SizeDiff, Read14DP, Read14%Tot, Read14Seq,Read15Length, Read15SizeDiff, Read15DP, Read15%Tot, Read15Seq,Read16Length, Read16SizeDiff, Read16DP, Read16%Tot, Read16Seq,Read17Length, Read17SizeDiff, Read17DP, Read17%Tot, Read17Seq,Read18Length, Read18SizeDiff, Read18DP, Read18%Tot, Read18,Read19Length, Read19SizeDiff, Read19DP, Read19%Tot, Read19SeqSeq,Read20Length, Read20SizeDiff, Read20DP, Read20%Tot, Read20Seq,Read21Length, Read21SizeDiff, Read21DP, Read21%Tot, Read21Seq,Read22Length, Read22SizeDiff, Read22DP, Read22%Tot, Read22Seq,Read23Length, Read23SizeDiff, Read23DP, Read23%Tot, Read23Seq,Read24Length, Read24SizeDiff, Read24DP, Read24%Tot, Read24Seq,Read25Length, Read25SizeDiff, Read25DP, Read25%Tot, Read25Seq\n"; 
#print OUT4 "Chr#, ChrMarker, Position, RefLength, RefMicrSat, UnqueReads, Genotype(length), Genotype, TotalDP, Read1Length, Read1DP, Read1%Tot, Read1Seq, Read2Length, Read2SizeDiff, Read2DP, Read2%Tot, Read2Seq, Read3Length, Read3SizeDiff, Read3DP, Read3%Tot, Read3Seq, Read4Length, Read4SizeDiff, Read4DP, Read4%Tot, Read4Seq, Read5Length, Read5SizeDiff, Read5DP, Read5%Tot, Read5Seq, Read6Length, Read6SizeDiff, Read6DP, Read6%Tot, Read6Seq, Read7Length, Read7SizeDiff, Read7DP, Read7%Tot, Read7Seq, Read8Length, Read8SizeDiff, Read8DP, Read8%Tot, Read8Seq, Read9Length, Read9SizeDiff, Read9DP, Read9%Tot, Read9Seq, Read10Length, Read10SizeDiff, Read10DP, Read10%Tot, Read10Seq,Read11Length, Read11SizeDiff, Read11DP, Read11%Tot, Read11Seq,Read12Length, Read12SizeDiff, Read12DP, Read12%Tot, Read12Seq,Read13Length, Read13SizeDiff, Read13DP, Read13%Tot, Read13Seq,Read14Length, Read14SizeDiff, Read14DP, Read14%Tot, Read14Seq,Read15Length, Read15SizeDiff, Read15DP, Read15%Tot, Read15Seq,Read16Length, Read16SizeDiff, Read16DP, Read16%Tot, Read16Seq,Read17Length, Read17SizeDiff, Read17DP, Read17%Tot, Read17Seq,Read18Length, Read18SizeDiff, Read18DP, Read18%Tot, Read18,Read19Length, Read19SizeDiff, Read19DP, Read19%Tot, Read19SeqSeq,Read20Length, Read20SizeDiff, Read20DP, Read20%Tot, Read20Seq,Read21Length, Read21SizeDiff, Read21DP, Read21%Tot, Read21Seq,Read22Length, Read22SizeDiff, Read22DP, Read22%Tot, Read22Seq,Read23Length, Read23SizeDiff, Read23DP, Read23%Tot, Read23Seq,Read24Length, Read24SizeDiff, Read24DP, Read24%Tot, Read24Seq,Read25Length, Read25SizeDiff, Read25DP, Read25%Tot, Read25Seq\n"; 
#print OUT3 join \t, @header;
print OUT3 "$header\n";

#print OUT4 join \t, @header;
print OUT4 "$header\n";

############ Chrom_Region_identifier file
while (my $line1 = <$fh1>) {
	chomp $line1;
	#$line1 =~ s/	/,/g; 
	
   	push (@position, $line1);
   
   }
   
print "Opened input files and obtained the data succesfully ...........\n\n";

print "Parsing and allocating data to files ...........\n\n";
   
$firstline = <$fh>;

while ($line = <$fh>) {
	chomp $line;
	$line =~ s/ //g;
	@line = split ('	', $line);
	
	
#add chromosome location 
		if (@line[0] eq chrM) {$location = 'MpM'; #print join ("\t", @line), "\n";
		} 
		elsif (@line[0] =~ m/_/){
			$location = "extra";
			}
		else {
			$location = Chrom (@line[0], @line[1]);						
			
		}
			
			splice @line, 1, 0, $location;
			

# remove any reads (alleles) that are below the set threshhold as set by $AlDP (-AD), default is at 3
for ($z=2;$z<25;$z++) {	
		$m = 3 + ($z * 3);			
		if ($line[$m] < $AlDP) {my $cnt = 0; #$n = $m;
				@line = grep { ++$cnt < $m } @line; last;} 
		 } 

# Add number of reads and unique reads	
	my $uniq=0;
	my $totalreads=0;
	my$x;
	my$i;

	for ($x=1;$x<25;$x++) {	
		$i = 3 + ($x * 3);	
		if ($line[$i]) { $uniq++; $totalreads+=$line[$i];} else {#$uniq; $totalreads; 
		next;}		
		} 	
	my $tot = $totalreads;
	
	
#add the percent of total
	for ($x=1;$x<25;$x++) {	
		my $a = 3 + ($x * 3) + ($x - 1);	
		if ($line[$a]) { $PercTotal = ($line[$a]/$tot) * 100; splice @line, $a+1, 0, $PercTotal;} else { next;}		
		} 	
		
#9, 14, 19		
				
#add size diff of microsat	
	my $tot = $totalreads;
	for ($x=2;$x<25;$x++) {	
		my $j = ($x * 5) - 1;	
		if ($line[$j]) { my $sizeDiff = ($line[$j]-$line[5]); splice @line, $j+1, 0, $sizeDiff;} else { next;}		
		} 			

splice @line, 5, 0, $uniq, $totalreads;		

#add genotype prediction to analysis		
my ($genotype, $type);
		if (!defined @line[13]) { $genotype = "@line[7].00"; $type = 1;}
		elsif (@line[13]/@line[8] > 0.25 ) {  $genotype = "@line[7].@line[11]"; $type = 2;}
		elsif (@line[13]/@line[8] < 0.25 ) {$genotype = "@line[7].00"; $type = 1; }
		
	
	
# add MST length and MST period

# chrM	208	13	GTGTG- TTAA-GCTTG	13	53	GTGTG-TTAATTAATTAAT-GCTTG
	
	@MST = split "-", $line[4];
	$MSTlength = length $MST[1];
	$MSTperiod = @line[3]/$MSTlength;
	
	splice @line, 6, 0, $genotype, $type;	
	splice @line, 4, 0, $MSTlength, $MSTperiod;		
	

# output into all_data_file (OUT10), below_cutoff (OUT1) and above_cutoff (OUT2)
	
		print OUT10 join ("	", @line), "\n";
		
		
		if ($MinDP > 0) {
			if (@line[10] >= $MinDP) {
				print OUT1 join ("	", @line), "\n";
			} else { print OUT2 join ("	", @line), "\n"; }
		}		


# output into above_cutoff_single (OUT3), above_cutoff_multi (OUT4)

		
		if (@line[10] >= $MinDP) {
			
			$line2 = join "	", @line;	
			if ($line2 =~ m/\|/ ) {print OUT3 "$line2\n";} 
   			elsif ($line2 =~ m/^Chr#/ or $line2 =~ m/^chrM/) {;}  	
   			else { print OUT4 "$line2\n"; 
   			#print "$line2\n";
   			}
   		
   		}

		
}


#subroutine for previous part

sub Chrom  {
	($chrom, $loc) = @_;
	
	foreach $element (@position) {
		@element = split (' ', $element);
		
		if ($chrom eq $element[1] &&  $loc < $element[3] && $loc > $element[2]) {
			$return =  "$element[4]";
			#print "-- $return\n"; exit;
			return $return;
			exit;
		} else { next; }
		
	}
	
}

close OUT1;
close OUT2;
close OUT10;
close OUT3;
close OUT4;


print "Data distribution into files complete .........\n";


########################## Summary tables ####################################


open(my $fh3, '<', $output1) or die "Can't read file '$file' [$!]\n";



my $output5 = "$file.DataSum.A$AlDP.DP$MinDP.csv";
open (OUT5, ">>$output5");

$firstline = <$fh3>;
my %chromo;
my %UnqReads;
my %sizediff;
my %readLen;
my %genoty;
while (<$fh3>)  {
    chomp $_;

my ($Chr, $ChrMarker, $Position, $RefLength, $MSTunitL, $Period, $RefMicrSat, $UnqueReads, 
$GenotypeLength, $Genotype, $TotalDP, $Read1Length, $Read1DP, 
$Read1perTot, $Read1Seq, $Read2Length, $Read2SizeDiff, $Read2DP, $Read2perTot, 
$Read2Seq, $Read3Length, $Read3SizeDiff, $Read3DP, $Read3perTot, 
$Read3Seq, $Read4Length, $Read4SizeDiff, $Read4DP, $Read4perTot, $Read4Seq, $Read5Length, 
$Read5SizeDiff, $Read5DP, $Read5perTot, $Read5Seq, $Read6Length, $Read6SizeDiff, 
$Read6DP, $Read6perTot, $Read6Seq, $Read7Length, $Read7SizeDiff, $Read7DP, $Read7perTot, 
$Read7Seq, $Read8Length, $Read8SizeDiff, $Read8DP, $Read8perTot, $Read8Seq, $Read9Length, 
$Read9SizeDiff, $Read9DP, $Read9perTot, $Read9Seq, $Read10Length, $Read10SizeDiff, 
$Read10DP, $Read10perTot, $Read10Seq,$Read11Length, $Read11SizeDiff, $Read11DP, 
$Read11perTot, $Read11Seq,$Read12Length, $Read12SizeDiff, $Read12DP, 
$Read12perTot, $Read12Seq,$Read13Length, $Read13SizeDiff, $Read13DP, $Read13perTot, 
$Read13Seq,$Read14Length, $Read14SizeDiff, $Read14DP, $Read14perTot, 
$Read14Seq,$Read15Length, $Read15SizeDiff, $Read15DP, $Read15perTot, 
$Read15Seq,$Read16Length, $Read16SizeDiff, $Read16DP, $Read16perTot, 
$Read16Seq,$Read17Length,$Read17SizeDiff, $Read17DP, $Read17perTot, 
$Read17Seq,$Read18Length, $Read18SizeDiff, $Read18DP, $Read18perTot, 
$Read18,$Read19Length, $Read19SizeDiff, $Read19DP, $Read19perTot, $Read19SeqSeq,
$Read20Length, $Read20SizeDiff, $Read20DP, $Read20perTot, $Read20Seq,$Read21Length, 
$Read21SizeDiff, $Read21DP, $Read21perTot, $Read21Seq,$Read22Length, $Read22SizeDiff, $Read22DP, 
$Read22perTot, $Read22Seq,$Read23Length, $Read23SizeDiff, $Read23DP, $Read23perTot, 
$Read23Seq,$Read24Length, $Read24SizeDiff, $Read24DP, $Read24perTot, $Read24Seq,$Read25Length, 
$Read25SizeDiff, $Read25DP, $Read25perTot, $Read25Seq) = split "	", $_; 


# Total counts and reads per chromosome
 
$chromo{$Chr}{count}++;
$chromo{$Chr}{DP} += $TotalDP;

# counts per unique reads
$UnqReads{$UnqueReads}{count}++;
$UnqReads{$UnqueReads}{DP} += $TotalDP;

# counts for size differences between microsats
$sizediff{$Read2SizeDiff}{count}++;
$sizediff{$Read3SizeDiff}{count}++;
$sizediff{$Read4SizeDiff}{count}++;
$sizediff{$Read5SizeDiff}{count}++;
$sizediff{$Read6SizeDiff}{count}++;
$sizediff{$Read7SizeDiff}{count}++;
$sizediff{$Read8SizeDiff}{count}++;
$sizediff{$Read9SizeDiff}{count}++;
$sizediff{$Read10SizeDiff}{count}++;
$sizediff{$Read11SizeDiff}{count}++;
$sizediff{$Read12SizeDiff}{count}++;
$sizediff{$Read13SizeDiff}{count}++;
$sizediff{$Read14SizeDiff}{count}++;
$sizediff{$Read16SizeDiff}{count}++;
$sizediff{$Read17SizeDiff}{count}++;
$sizediff{$Read18SizeDiff}{count}++;
$sizediff{$Read19SizeDiff}{count}++;
$sizediff{$Read20SizeDiff}{count}++;
$sizediff{$Read21SizeDiff}{count}++;
$sizediff{$Read22SizeDiff}{count}++;
$sizediff{$Read23SizeDiff}{count}++;
$sizediff{$Read24SizeDiff}{count}++;
$sizediff{$Read25SizeDiff}{count}++;

# counts per unique reads
$readLen{$Read1Length}{count}++;
$readLen{$Read1Length}{DP} += $TotalDP;

# hetero or homo-zygotic loci
$genoty{$Genotype}{count}++;
$genoty{$Genotype}{$UnqReads}++;

}

print OUT5 "The input file was $output1 
The settings were (Allele Depth) -A $AlDP and -DP $MinDP";

print  OUT5 "Total counts and reads per allele per chromosome\n";
print OUT5 "Chrom#,Counts,MeanDP\n";
hash_counts (%chromo);


print OUT5 "Total counts of nmber of locus \n";
print OUT5 "Allele #,Counts,averageDP/locus\n";
hash_counts (%UnqReads);



$counts=0;
$depth=0;
print OUT5 "Counts of size differences from major alleles \n";
print OUT5 "Allele #	sizeDiff\n";
foreach $occurs (sort keys %sizediff) {
	if ($occurs =~ m/\d/) {
    	print OUT5 "$occurs	$sizediff{$occurs}{count}\n";
    	$counts += $sizediff{$occurs}{count};
	}
}

print  OUT5 "Total	$counts\n\n";


print OUT5 "Total counts per read length \n";
print OUT5 "ReadLength,Counts,averageDP/locus\n";
hash_counts (%readLen);








sub hash_counts {
	(%hash) = @_;
	$counts=0;
	$depth=0;	
	foreach $occurs (sort keys %hash) {
    print OUT5 $occurs,"	",$hash{$occurs}{count},"	",$hash{$occurs}{DP}/$hash{$occurs}{count},"\n";
    $counts += $hash{$occurs}{count};
    $depth += $hash{$occurs}{DP};
    
}

print OUT5 "Total	$counts	",$depth/$counts,"\n\n";
}



