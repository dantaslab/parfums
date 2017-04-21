#!/usr/bin/perl -w

# Written by Alejandro Reyes

use strict;

if (@ARGV != 1) {
	die "\nUsage: GetAssembly-Illumina.pl <InputFile>\n";
}

my %adapter=();
my %total=();
my $infile = shift;
open (IN, "<$infile") or die ("Couldn't open file: $infile\n");
$infile=~ s/.fna//;

open (OUT, ">$infile.noVector") or die ("Couldn't open file: $infile.clean\n");
my $pass=0;

my $switch = 0;
while (my $li = <STDIN>){
  $switch =1 if ($li =~ /Num. pairs/);
  if ($switch == 1 && $li =~ /^\s+\d+/){
    my @arr = split /\s+/, $li;
    $arr[8]=~ /\((\d+)\)/;
    $pass=0;
    $pass=1 if ($arr[6] <= 5);
    my $left = $1;
    $pass=2 if ($left <= 5);
    $pass=3 if (abs($arr[7]-$arr[6])/($arr[7]+$left) > 0.85);
    $pass=4 if ($total{$arr[5]});
    $total{$arr[5]}=1;

    next unless $pass > 0;
    $arr[5]=~ s/\_\d$//;
    $adapter{$arr[5]}=1;
  }
  $switch=0 if ($li =~ /^\d+\s+matching entries/);
}


my $seq1;
my $seq2;
my $name1;
my $name2;

while (my $line=<IN>){
  chomp $line;
  
  if (($line =~ /^>(.*)\_0\s*$/)){
    my $n = $1;
    $name1=$line;
    $seq1=<IN>;
    $name2=<IN>;
    $seq2=<IN>;
    chomp $seq1;
    chomp $name2;
    chomp $seq2;
    
    die ("Second sequence doesn't match name\n") unless ($name2 =~ /$n/);
 
    print OUT "$name1\n$seq1\n$name2\n$seq2\n" unless ($adapter{$n});
  }
}
close IN;
close OUT;
