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
  if ($j =~ /^>(.*)$/){
    $name=$1;
    $name =~ s/\s+//;
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
  next if ($array[7] < 95);                                  # Next if the percent id is less than 95%
  next unless ((($array[3]/$array[1]) >= 0.95) or ($array[-2]< 5) or ($array[-1] > ($length{$array[8]}-5)));                     # Next unless the percent covered is over 95% of the read or is in the first or last 30bp
  for (my $i=$array[-2]; $i<=$array[-1]; $i++){
    $genomes{$array[8]}[$i-1][0]++;          # Marks the coverage of all sequences
    if ($i > ($array[-2]+10) && $i < ($array[-1]-10)){
       $genomes{$array[8]}[$i-1][1]++;
    }elsif (($i <= 12) && $i < ($array[-1]-10)){
       $genomes{$array[8]}[$i-1][1]++;
    }elsif (($i > ($array[-2]+10)) && ($i > ($length{$array[8]}-15))){
       $genomes{$array[8]}[$i-1][1]++;
    }
  }
}
close IN;

foreach my $k (keys(%genomes)){
  die ("Unknown genome >$k<\n") unless ($length{$k});
 # print "\n$k\n";
  #next unless ($k =~ /NODE_1_length_2354_cov_40.158878/);
  for (my $i=0; $i<$length{$k}; $i++){
    $genomes{$k}[$i][0]=0 unless ($genomes{$k}[$i][0]);
    $genomes{$k}[$i][1]=1 unless ($genomes{$k}[$i][1]);
    $genomes{$k}[$i][2]=$genomes{$k}[$i][0]/$genomes{$k}[$i][1];
   # print "$i\t$genomes{$k}[$i][0]\t$genomes{$k}[$i][1]\t$genomes{$k}[$i][2]\n";
    substr $seqs{$k}, $i, 1, "N" if ((($genomes{$k}[$i][1] < 20) && ($genomes{$k}[$i][2] > 10)) or ($genomes{$k}[$i][2] > 100)); 
  }
  
  my @split = split /N+/, $seqs{$k};
  for (my $h =0; $h<@split; $h++){
    print ">$k\_$h\n$split[$h]\n" if (length($split[$h]) > 50);
  }
}



