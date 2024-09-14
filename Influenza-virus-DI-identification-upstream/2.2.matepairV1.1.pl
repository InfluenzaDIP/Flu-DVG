#!/usr/bin/perl -w
#Author: Jia Xinmiao
#Email: jiaxm@big.ac.cn
#Date: 2015-07-21
#Description: reads配对
#reads ID @ST-E00144:122:HF2YVCCXX:4:1101:7303:1801 1:N:0:=
my $version=1.00;

use strict;
use Getopt::Long;

my %opts;
GetOptions(\%opts,"f=s","r=s","om=s","os=s","h");
if (!(defined $opts{f} and defined $opts{r} and defined $opts{om} and defined $opts{os}) || defined $opts{h}) { #necessary arguments
&usage;
}

my $filef=$opts{'f'};
my $filer=$opts{'r'};

open IN,"<$filef";
my %hash;
while (my $aline=<IN>) {
	chomp $aline;
	if ($aline=~/^(@\S+)/) {$hash{$1}=1;}
}
close IN;

my @name;
my %mateID;
open I,"<$filer";#input file  
while (my $aline=<I>) {
	chomp $aline;
	if ($aline=~/^(@\S+)/) {
		if (exists $hash{$1}) {
			$mateID{$1}="mate";#将一致的ID存入哈希表
			#push @name,$1;#存入数组
		}
	}
}
close I;

my $outf3=$opts{'om'}."_1_mt.fastq";
my $outr3=$opts{'om'}."_2_mt.fastq";
my $outos=$opts{'os'}."_single.fastq";
open F,">$outf3"; #output file  
open R,">$outr3";
open S,">$outos";

open IN ,"<$filef";
open I,"<$filer";

my ($nameF,$nameR,$seqF,$seqR,$nF,$nR,$qvF,$qvR);

my $numf=0;
my $numr=0;
my $nums1=0;
my $nums2=0;

while (!eof(IN)) {
	$nameF=<IN>;
	$seqF=<IN>;
	$nF=<IN>;
	$qvF=<IN>;
	chomp $nameF;
	$nameF=~/^(@\S+)/;
	if (exists $mateID{$1}) {
		print F "$nameF\n$seqF$nF$qvF";
		$numf += 1;
	}
	else{
		print S "$nameF\n$seqF$nF$qvF";
		$nums1 += 1;
	}

}

while (!eof(I)) {
	$nameR=<I>;
	$seqR=<I>;
	$nR=<I>;
	$qvR=<I>;
	chomp $nameR;
	$nameR=~/^(@\S+)/;
	if (exists $mateID{$1}) {
		print R "$nameR\n$seqR$nR$qvR";
		$numr += 1;
	}
	else{
		print S "$nameF\n$seqF$nF$qvF";
		$nums2 += 1;
	}
}

print "Mate pair reads number: $numf \n";

print "Not paired: reads1 $nums1 ; reads2 $nums2 \n";

close F;
close R;
close IN;
close I;
sub usage{
print <<"USAGE";
Version $version
Usage:
$0 -f -r  -om -os
options:
-f input fastq1 file
-r input fastq2 file
-om output file prefix (paired)
-os output file prefix (notpaired)
-h help
USAGE
exit(1);
}

