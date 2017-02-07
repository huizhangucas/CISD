#!/usr/bin/perl -w

use strict;
use warnings;

my $inputpeak=$ARGV[0];
my $outputsite=$ARGV[1];
open(PEAK,"< $inputpeak");
open(LOG,"> $outputsite");

#open(PEAK,"<hspeaks0.95_chrAllnew.bed");
#open(LOG,">candidate_sites");

print "converting high score peaks to candidate sites\n";

my $k=0;
while (<PEAK>){
$k++;
my $line=$_;
chomp($line);

my ($chr,$start,$end)=split(/\t/,$line);
if ($end-$start >2001){
	$start=$start+1000;
	$end=$end-1000;
my $center=($start+$end)/2;
my $start2=$start+($center-$start)%100;
my $end2=$end-($center-$start)%100;
my $i;
for ($i=$start2;$i<=$end2;$i=$i+100){
	my $a=$i-100;
	my $b=$i+100;	
	print LOG "$chr\t$a\t$b\n";
}
}else{
	my $center=int(($start+$end)/2);
	my $a=$center-100;
	my $b=$center+100;
	print LOG "$chr\t$a\t$b\n";
#	print LOG "$line\n";
}

}
print "candidate sites are saved to $outputsite\n";
