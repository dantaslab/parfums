#!/usr/bin/perl -w

# Written by Alejandro Reyes

use strict;

if (@ARGV != 1) {
	die "\nUsage: GetAssembly-Illumina.pl <InputFile>\n";
}

my $infile = shift;
open (IN, "<$infile") or die ("Couldn't open file: $infile\n");

my $name1="";
my %seqs=();
while (my $line=<IN>){
  chomp $line;
  if ($line =~ /^>(.*)$/){
    $name1=$1;
    $name1=~ s/\s*$//;
  }else{
    $seqs{$name1}.=$line;
  }
}
close IN;

while (my $li = <STDIN>){
  my @arr = split /\s+/, $li;
  for (my $i=$arr[6]-1; $i < $arr[7]; $i++){
    substr $seqs{$arr[0]}, $i, 1, "N";
  }
}

foreach my $k (keys(%seqs)){
  my $count=1; 
  my @sub=split /N+/, $seqs{$k};
  foreach my $l (@sub){
    print ">$k\_$count\n$l\n" unless (length($l) < 50);
  }
}
