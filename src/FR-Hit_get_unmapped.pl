#!/usr/bin/perl -w
use strict;



die ("Usage: FR-Hit_get_unmapped.pl <Output_FR-Hit> <FASTA_FILE> <prefix>\n") unless (scalar(@ARGV) == 3);

use constant WAIT_TIME => '60' ; # 60 seconds = 1 minute
use constant NO_CHANGE_IN_FILE_FOR_X_SECONDS => '30';



my $file = shift @ARGV;
my $fastaFile = shift @ARGV;
my $prefix = shift @ARGV;
my $x = 0;
my $dir=`pwd`;
chomp $dir;
my $nd= "$dir/$prefix";

my @seqs=();

open (IN, "<$file") or die ("Couldn't open file: $file\n");
open (OUT, ">$prefix\_ToGet_missing.txt") or die ("Couldn't open file: $prefix\_ToGet_missing.txt\n");

my $old="";
my %used=();
my %names = ();

while (my $l = <IN>){
  chomp ($l);
  my @array = split /\t/, $l;

  $array[7] =~ s/\%$//;
  $array[1] =~ s/nt$//;

  if ($array[0] ne $old){
    if ($old ne ""){
      &process_array(\@seqs);
    }
    $old = $array[0];
    @seqs=();
  }
  push (@seqs, [@array]);
}
close IN;
&process_array(\@seqs);

foreach my $l (keys(%used)){
  my $read_name1 = "$l\_0";
  my $read_name2 = "$l\_1";
  print OUT "$read_name1\n$read_name2\n";
  $names{$read_name1} = 1;
  $names{$read_name2} = 1;
}
close OUT;

#Reads the Fasta Sequences

open (IN, "<$fastaFile") or die ("Couldn't open file: $fastaFile\n");

my $counter = 0;
while (my $line = <IN>) {
  if ($line =~ /^>/){
    $counter = 0;
    chomp $line;
    $line =~ s/^>//;			#Only if the name does not contain >, otherwise comment this line
    my @temp = split (/\s+/, $line);

    unless ($names{$temp[0]}){
      print ">$line\n";
      $counter = 1;
    }
  }

  elsif (($line =~ /^\w/) && ($counter == 1)){
    print $line;
  }

}
close IN;

#system("/home/gdlab/shared/old-shared/scripts/PARFuMS_scripts/GetSequences_inverted.pl $prefix\_ToGet_missing.txt $mapFile > $prefix\_MissingFromFirstPass.fna");


########### subroutines #######

sub process_array{
  my @arr = @{$_[0]};

  my $sum=0;
  my @seq=();
  for (my $j=0; $j<@arr; $j++){
    my @coords= sort ($arr[$j][4], $arr[$j][5]);
    for (my $i=$coords[0]; $i<=$coords[1]; $i++){
      $seq[$i]++;
    }
  }
  for (my $i=1; $i<=$arr[0][1]; $i++){
    $sum++ if ($seq[$i]);
  }
  my $name = $arr[0][0];
  $name=~ s/\_\d$//;
  $used{$name}=1 if ($sum > ($arr[0][1]*0.8));
#  print "Sum is $sum while size is $arr[0][1]\n";
}
