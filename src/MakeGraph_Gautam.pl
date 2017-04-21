#!/usr/bin/perl -w

use strict;

die ("Usage: MakeGraph_Gautam.pl <Contig_File> <Blast_File>  > GraphOutfile\n") unless scalar(@ARGV)==2;


# global variables for the different output file formats


my $genome_file = shift;

open (IN, $genome_file) or die ("Couldn't open file $genome_file\n");

my %genomeSize=();
my %genomeVector=();
my $name="";
my $seq="";

while (my $line =<IN>){
  chomp $line;
  
  if ($line =~ /^>/){
    $line =~ s/>//;
    if ($name ne ""){
      $genomeSize{$name}=length($seq);
      &lookVector($seq, $name);
    }
    $name = $line;
    $name =~ s/\s+$//;
    $seq="";
  }else{
    $seq.=$line;
  }
}
$genomeSize{$name}=length($seq);
&lookVector($seq, $name);
close IN;



my $blast_file = shift;

open (IN2, $blast_file) or die ("Couldn't open file $blast_file\n");

my %genomeCounts;
my %linkCounts;
my %linkLabels;
my %prevLine=();
my %b=();
my $prev_id="blah";

while (<IN2>) {
	chomp;

	#my ($score, $startG, $stopG, $startQ, $stopQ, $dir, $name, $target)

	my @line = parse_result_line($_);	
	my $id = $line[-2];
	$id =~ s/\_[0,1]$//;
	
	if ($id eq $prev_id){
	  $b{$line[-2]}=$line[0] unless $b{$line[-2]};
	  push (@{$prevLine{$line[-2]}}, [@line]) if ($line[0]>0.99*$b{$line[-2]});
	}else{
	  &processline(\%prevLine) if (scalar(keys (%prevLine)) == 2);
	  %prevLine=();
	  %b=();
	  $b{$line[-2]}=$line[0];
	  push (@{$prevLine{$line[-2]}}, [@line]);
	  $prev_id=$id;
	}
}
close IN2;



print "digraph G \{\n\trankdir=LR\;\n\tsize=\"8,5\"\n\tnode \[shape = circle\]\;\n";


my %contigMedian=();

#Foreach contig calculate the meadian coverage per bp
foreach my $k (keys (%genomeCounts)){
  my @count_array=();
  for (my $i=0; $i<$genomeSize{$k}; $i++) {
    if ($genomeCounts{$k}[$i]){
      push @count_array, $genomeCounts{$k}[$i];
    }else{
      push @count_array,"0";
    }
  }
  $contigMedian{$k}= &median(\@count_array);
  print "\t$k \[label=\"$k\\nCov:$contigMedian{$k} Len:$genomeSize{$k}\\n$genomeVector{$k}\"\]\;\n";
}


foreach my $k (keys (%genomeCounts)){
  my $cum=0;
  for (my $i=0; $i<$genomeSize{$k}; $i++) {
    $cum += $genomeCounts{$k}[$i] if ($genomeCounts{$k}[$i]);		
  }
  my $avg=int($cum/$genomeSize{$k});
}

foreach my $s (keys (%linkCounts)){
  foreach my $t (keys (%{$linkCounts{$s}})){
    print "\t$s -> $t \[label=\"$linkLabels{$s}{$t}\\n$linkCounts{$s}{$t}\" \] \;\n" if ($linkCounts{$s}{$t} > 10);
  }
}print "\}\n";


##############SubRoutines###########

sub median{
  my ($array_ref) = @_;
  my $count = scalar @$array_ref;
  my @array = sort { $a <=> $b } @$array_ref;
  if ($count % 2) {
    return $array[int($count/2)];
  } else {
    return ($array[$count/2] + $array[$count/2 - 1]) / 2;
  }
} 

sub processline{
  my %lines = %{$_[0]};
  my @l=keys(%lines);
  
  die ("No tiene 2 keys\n") if (scalar(@l) != 2);

  foreach my $k (@{$lines{$l[0]}}){
    foreach my $t (@{$lines{$l[1]}}){
      my @a=@{$k};
      my @b=@{$t};
      next if ($a[-1] eq $b[-1]);

      #Here I need to see if the total length between the links is less than 300 make it;
      my $sizeA = 0;
      my $sizeB = 0;
      my $dirA="";
      my $dirB="";
      die ("No existe tamanio para genoma A >$a[-1]<\n") unless $genomeSize{$a[-1]};
      die ("No existe tamanio para genoma B >$b[-1]<\n") unless $genomeSize{$b[-1]};

      if ($a[-3] == 0){
	$sizeA=$genomeSize{$a[-1]}-$a[1];
	$dirA="5'";
      }else{
	$sizeA=$a[2];
	$dirA="3'";
      }
      if ($b[-3] == 0){
	$sizeB=$genomeSize{$b[-1]}-$b[1];
	$dirB="5'";
      }else{
	$sizeB=$b[2];
	$dirB="3'";
      }      
      next unless ($sizeA + $sizeB) < 300;


      if ($linkCounts{$a[-1]}{$b[-1]}){
	$linkCounts{$a[-1]}{$b[-1]}++;
	$linkLabels{$a[-1]}{$b[-1]}="$dirA\_$a[-1]\_to_$dirB\_$b[-1]" unless $linkLabels{$a[-1]}{$b[-1]};
      }else{
	$linkCounts{$b[-1]}{$a[-1]}++;
	$linkLabels{$b[-1]}{$a[-1]}="$dirB\_$b[-1]\_to_$dirA\_$a[-1]" unless $linkLabels{$b[-1]}{$a[-1]};
      }
    }
  }
  
  for (my $j=0; $j<@{$lines{$l[0]}}; $j++){
    for my $t ($lines{$l[0]}[$j][1] .. $lines{$l[0]}[$j][2]){
      $genomeCounts{$lines{$l[0]}[$j][-1]}[$t]++;
    }
  }
  for (my $j=0; $j<@{$lines{$l[1]}}; $j++){
    for my $t ($lines{$l[1]}[$j][1] .. $lines{$l[1]}[$j][2]){
      $genomeCounts{$lines{$l[1]}[$j][-1]}[$t]++;
    }
  }

}	  

sub parse_result_line {
  my ($line) = @_;
  my @pieces = split /\t/;
  my $dir=0;
  
  #	return ($score, $startG, $stopG, $startQ, $stopQ, $name, $target)
  
  if ($pieces[8] > $pieces[9]) { # infer direction from the orientation of the match positions and put smallest first
    $dir=1;
    my $temp = $pieces[8];
    $pieces[8] = $pieces[9];
	  $pieces[9] = $temp;
  }
  return $pieces[11], $pieces[8], $pieces[9], $pieces[6], $pieces[7], $dir, $pieces[0], $pieces[1];
}


sub lookVector{
  my ($seq, $name)=@_;
  
  my $return="";
  my $res="";

  $res=`echo $seq | tre-agrep -2 -i ATCAAGCTTATCGATACCGTC`;
  $res.=`echo $seq | tre-agrep -2 -i GACGGTATCGATAAGCTTGAT`;
  if ($res ne ""){
    $return.=" Vector_3";
    $res="";
  }
  $res=`echo $seq | tre-agrep -2 -i GACCTCGAGGGGGGG`;
  $res.=`echo $seq | tre-agrep -2 -i CCCCCCCTCGAGGTC`;
  if ($res ne ""){
    $return.=" Vector_5";
    $res="";
  }

  $genomeVector{$name}=$return;

}
