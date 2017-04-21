#!/usr/bin/perl -w

# Written by Alejandro Reyes

use strict;

if (@ARGV != 1) {
	die "\nUsage: GetAssembly-Illumina.pl <InputFile>\n";
}

my %seqs=();
my %adapter=();
my %vector=();
my $infile = shift;
open (IN, "<$infile") or die ("Couldn't open file: $infile\n");

$infile=~ s/.fna//;

open (OUT, ">$infile.clean") or die ("Couldn't open file: $infile.clean\n");

my $switch = 0;
while (my $li = <STDIN>){
  $switch =1 if ($li =~ /Num. pairs/);
  if ($switch == 1 && $li =~ /^\s+\d+/){
    my @arr = split /\s+/, $li;
    $arr[8]=~ /\((\d+)\)/;
    my $n=$1;
    if ($li =~ /Adapter/){
      $arr[5]=~ s/\_\d$//;
      next unless ($arr[-3] !~ /\(/ && $arr[-3] < 6);
      if ($arr[6] < 35){
	$seqs{$arr[5]}=1;
      }else{
	$adapter{$arr[5]}=$arr[6]-1;
      }
    }
    if ($li =~ /Vector/){
      if ($arr[6] <= 5){
	$vector{$arr[5]}[0]=0;
	$vector{$arr[5]}[1]=$arr[7];
      }
      if ($n <= 5){
	$vector{$arr[5]}[0]=$arr[6]-1;
	$vector{$arr[5]}[1]=$arr[7]+$n;	
      }else{
	$arr[5]=~ s/\_\d$//;
	$seqs{$arr[5]}=1;
      }
    }
  }
  $switch=0 if ($li =~ /^\d+\s+matching entries/);
}


my $seq1;
my $seq2;
my $name1;
my $name2;
my $sevN=0;
my $numN=0;

while (my $line=<IN>){
  chomp $line;
  
  if (($line =~ /^>(.*)\_0\s*$/)){
    $numN=0;
    $sevN=0;
    my $n = $1;
    $name1=$line;
    $name1=~ s/^>//;
    $seq1=<IN>;
    $name2=<IN>;
    $name2=~ s/^>//;
    $seq2=<IN>;
    chomp $seq1;
    chomp $name2;
    chomp $seq2;
        
    die ("Second sequence doesn't match name\n") unless ($name2 =~ /$n/);
    $sevN = 1 if ($seq1 =~ /NNNNN+/ or $seq2 =~ /NNNN+/);
    $numN = $seq1 =~ s/N/N/g;
    $numN += $seq2 =~ s/N/N/g;
    next if ($seqs{$n} or $numN > 6 or $sevN);

    if ($vector{$name1}){
      for (my $i=$vector{$name1}[0]; $i<$vector{$name1}[1]; $i++){
	substr $seq1, $i, 1, "X";
      }
    }if ($vector{$name2}){
      for (my $i=$vector{$name2}[0]; $i<$vector{$name2}[1]; $i++){
	substr $seq2, $i, 1, "X";
      }
    }if ($adapter{$n}){
      $seq1 = substr $seq1, 0, $adapter{$n};
      $seq2 = substr $seq2, 0, $adapter{$n};
    }
    $seq1 =~ s/X+//;
    $seq2 =~ s/X+//;
    print OUT ">$name1\n$seq1\n>$name2\n$seq2\n";
  }
}
close IN;
close OUT;
