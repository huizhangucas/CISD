#!/usr/bin/perl -w

use strict;
use warnings;

my $loopfile1=$ARGV[0];#find which loops in loop1 are supported by loop2 ,loop1 may be predicted loop, loop2 may be chia-pet loop
my $loopfile2=$ARGV[1];#
my $outputfile=$ARGV[2];#
my $tmpdir="./tmp_loop_overlap_loop";
my $loop1;
my $loop2a;
my $loop2b;

system("mkdir -p $tmpdir");
system("rm -rf $tmpdir/tmp*");
system("intersectBed -a $loopfile1 -b $loopfile2 -wa -wb > $tmpdir/tmp0");
system("cut -f 4-6 $tmpdir/tmp0 > $tmpdir/tmp046");
system("paste $tmpdir/tmp046 $tmpdir/tmp0 > $tmpdir/tmp1");
system("awk 'BEGIN {OFS=\"\\t\"}{print \$4,\$5,\$6,\$1,\$2,\$3}' $loopfile2 > $tmpdir/tmp_loop2reverse ");
system("intersectBed -a $tmpdir/tmp1 -b $tmpdir/tmp_loop2reverse -wa -wb | cut -f 4-21 > $tmpdir/tmp2");

open(FILE,"<$tmpdir/tmp2");
open(LOG,">$tmpdir/tmp_result");
while(<FILE>){
	chomp;
	my @array=split;
	$loop1="$array[0]\t$array[1]\t$array[2]\t$array[3]\t$array[4]\t$array[5]";
	$loop2a="$array[6]\t$array[7]\t$array[8]\t$array[9]\t$array[10]\t$array[11]";
	$loop2b="$array[15]\t$array[16]\t$array[17]\t$array[12]\t$array[13]\t$array[14]";
	if ($loop2a eq $loop2b) {
		print LOG "$loop2a\t$loop1\n";#### first output loop2£¬then out put loop1
	}
}
close(LOG);
system("cut -f 7-12 $tmpdir/tmp_result | sort -u | sort -k2n -o $outputfile");
#system("rm -rf $tmpdir");
