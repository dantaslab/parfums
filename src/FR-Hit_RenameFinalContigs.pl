#!/usr/bin/perl -w
use strict;



die ("Usage: FR-Hit_cleanChimera.pl <Output_FR-Hit> <Contigs> > <output>\n") unless (scalar(@ARGV) == 2);



my $file = shift @ARGV;
my $fna = shift @ARGV;
my $prefix = shift @ARGV;
my @seqs=();

open (IN, "<$file") or die ("Couldn't open file: $file\n");
open (FNA, "<$fna") or die ("Couldn't open file: $fna\n");

my %length=();
my %seqs=();
my $name="";

while (my $j = <FNA>){
  chomp($j);
  if ($j =~ /^>/){
    my @nameTable = split /\s/, $j;
    $name=$nameTable[0];
    $name =~ s/\s+//;
    $name =~ s/^>//;
  }else{
    $seqs{$name}.=$j;
    $length{$name}=length($seqs{$name});
  }
}

close FNA;


my %genomes=();
my $old="";
my $n="";

while (my $l = <IN>){
  chomp ($l);
  my @array = split /\t/, $l;
  $n = $array[0];
  $n =~ s/\_\d$//;

  $array[7] =~ s/\%$//;
  $array[1] =~ s/nt$//;
  next if ($array[7] < 85 or  ($array[3]/$array[1]) < 0.9);
  for (my $i=$array[-2]; $i<=$array[-1]; $i++){
    $genomes{$array[8]}[$i]++;                                 # Gives the coverage of the genome
  }
}


my $count=1;
foreach my $k (keys(%genomes)){
  die ("Unknown genome >$k<\n") unless ($length{$k});
  for (my $i=0; $i<=$length{$k}; $i++){
	$genomes{$k}[$i]=0 unless ($genomes{$k}[$i]);
  }
  my @sort= sort(@{$genomes{$k}});
  my $mean = $sort[int($length{$k}/2)];

  print ">Contig\_$count\_Mean:$mean\_Len:$length{$k}\n$seqs{$k}\n";
  $count++;
}



