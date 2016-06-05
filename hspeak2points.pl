#!/usr/bin/perl -w

use strict;
use warnings;

$inputpeak=$ARGV[0];
$outputsite=$ARGV[1];
open(PEAK,"< $inputpeak");
open(LOG,"> $outputsite");

#open(PEAK,"<hspeaks0.95_chrAllnew.bed");
#open(LOG,">candidate_sites");

print "converting high score peaks to candidate sites\n";

$k=0;
while (<PEAK>){
$k++;
$line=$_;
chomp($line);

my ($chr,$start,$end)=split(/\t/,$line);
if ($end-$start >2001){
	$start=$start+1000;
	$end=$end-1000;
my $center=($start+$end)/2;
$start2=$start+($center-$start)%100;
$end2=$end-($center-$start)%100;
for ($i=$start2;$i<=$end2;$i=$i+100){
	$a=$i-100;
	$b=$i+100;	
	print LOG "$chr\t$a\t$b\n";
}
}else{
	my $center=int(($start+$end)/2);
	$a=$center-100;
	$b=$center+100;
	print LOG "$chr\t$a\t$b\n";
#	print LOG "$line\n";
}

}
print "candidate sites are saved to $outputsite\n";
