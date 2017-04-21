#!/usr/bin/perl -w
use strict;

die ("Usage: FR-Hit_clean.pl <Output_FR-Hit> <ToGet_Seqs>\n") unless (scalar(@ARGV)==2);

my $file = shift @ARGV;
my $file2 = shift @ARGV;

open (IN, "<$file") or die ("Couldn't open file: $file\n");
open (OUT, ">$file2") or die ("Couldn't open file: $file2\n");

my %used=();

while (my $l = <IN>){
  chomp ($l);
  my @array = split /\t/, $l;
  $array[0] =~ s/\_\d$//;
  next if ($used{$array[0]});
  next if ($array[6] =~/\+/ && $array[8] =~/start/);
  next if ($array[6] =~/\-/ && $array[8] =~/stop/);
  next if ($array[8] =~/stop/ && $array[-1] < 90);
  next if ($array[8] =~/start/ && $array[-2] > 10);
  $used{$array[0]}=1;
  print OUT "$array[0]_0\n$array[0]_1\n";
}
close IN;
close OUT;
