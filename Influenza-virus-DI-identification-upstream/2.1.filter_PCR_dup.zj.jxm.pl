#!/usr/bin/perl
#Author: Jun Zhong (钟君) <joshchung@foxmail.com>
#Copyright: The script is for personal communication only and please do not leak out without onwer's permission.
#版权说明：此脚本仅限于作者个人交流，如没有取得作者允许请不要扩散。
#Institute: Beijing Institute of Genomics, Chinese Academy of Sciences
#Description:
# This program .
#History:
# 2011/09/20 V 1.0  
#Xinmiao Jia modified 
#2015/07/17
use strict;
use warnings;

#########modified context
use Getopt::Long;
my %opts;
GetOptions(\%opts,"fq1:s","fq2:s","outprefix:s");
my $usage= <<"USAGE";
	Program: $0
	INPUT:
		-fq1        Input filter fastq file1
		-fq2        Input filter fastq file2
	OUTPUT:
		-outprefix        Output filter fastq result file prefix
USAGE
die $usage unless(defined $opts{fq1} and $opts{fq2} and $opts{outprefix});

my $forwardFQ=$opts{fq1};
my $reverseFQ=$opts{fq2};
my $outputP=$opts{outprefix};
########modified context

my ($name, %name_of, %define_name, $num_1, $num_2,$name_f, $name_r);

my $cut_size = 25;#截取的长度

open F1, "$forwardFQ" || die "$!";;
open F2, "$reverseFQ" || die "$!";;

while (!eof(F1) && !eof(F2)) {
	chomp(my $first_line_f = <F1>);
	chomp(my $first_line_r = <F2>);
	if ($first_line_f =~ /^(@\S+)/) {#modified according to the reads name
		$name_f = $1;
	}
	else {
		print "First line of forward reads file not the name!\n";
		die;
	}
	if ($first_line_r =~ /^(@\S+)/) {#modified according to the reads name
		$name_r = $1;
	}
	else {
		print "First line of reverse reads file is not the name!\n";
		die;
	}
	if ($name_f ne $name_r) {
		print "The name of forward read is different from reverse!\n";
		print "forward is /$first_line_f backward is /$first_line_r \n";
		die;
	}
	$name = $name_f;

	chomp(my $seq_f = <F1>);
	chomp(my $seq_r = <F2>);
	my $pre_f = substr($seq_f, 0, $cut_size);
	my $pre_r = substr($seq_r, 0, $cut_size);
	my $new = $pre_f . $pre_r;	
	if (length($new) != ($cut_size*2)) {
		print "\$new is not ($cut_size*2)!!!\n";
		die;
		
	}
	$name_of{$new} = $name;

	<F1>; <F2>;
	<F1>; <F2>;

	$num_1++;
}

foreach  (values %name_of) {
	$define_name{$_} = 1;
}
close F1;
close F2;


open F1, "$forwardFQ" || die "$!";
open F2, "$reverseFQ" || die "$!";

#####modified output name
my $outputP1="$outputP"."_1";
my $outputP2="$outputP"."_2";
####

open F1_FILTER, ">$outputP1.rmdup.fq";
open F2_FILTER, ">$outputP2.rmdup.fq";

while (!eof(F1) && !eof(F2)) {
	my $first_line_f = <F1>;
	my $first_line_r = <F2>;
	my $seq_f = <F1>;
	my $seq_r = <F2>;
	my $plus_f = <F1>;
	my $plus_r = <F2>;
	my $qual_f = <F1>;
	my $qual_r = <F2>;
	if ($first_line_f =~ /^(@\S+)/) {#modified according to the reads name
		if (exists $define_name{$1} ) {
			$num_2++;
			print F1_FILTER "$first_line_f$seq_f$plus_f$qual_f";
		}
	}
	if ($first_line_r =~ /^(@\S+)/) {#modified according to the reads name
		if (exists $define_name{$1} ) {
			print F2_FILTER "$first_line_r$seq_r$plus_r$qual_r";
		}
	}

}
open F3, ">$outputP.dup.count";
print F3 "Num of total reads: $num_1\n";
my $num_dup = ($num_1 - $num_2);
my $account_for = sprintf "%.2f%%", ($num_dup/$num_1)*100;
print F3 "Num of dup reads: $num_dup\nDup rate: $account_for";

#added 
close F1;
close F2;
close F1_FILTER;
close F2_FILTER;
close F3;

