#!/usr/bin/perl -w

use strict;
use POSIX qw(log10);


die ("Usage: Stich_Contigs_Gautam.pl <Contig_File> <Graph_File> <AnnotationFile>  > GraphOutfile\n") unless scalar(@ARGV)==3;


# global variables for the different output file formats


my $genome_file = shift;

open (IN, $genome_file) or die ("Couldn't open file $genome_file\n");

my %contigs=();
my %conLen=();
my $name="";

while (my $line =<IN>){
  chomp $line;
  
  if ($line =~ /^>/){
    $line =~ s/>//;
    $name = $line;
    $name =~ s/\s+$//;
  }else{
    $contigs{$name}.=$line;
    $conLen{$name}=length($contigs{$name});
  }
}
close IN;



my $blast_file = shift;

open (IN2, $blast_file) or die ("Couldn't open file $blast_file\n");

my %ContigInfo=();
my %linkCounts=();

while (my $l=<IN2>) {
  next unless ($l =~ /label/);
  chomp $l;
  $l =~ s/\[label\s*=\s*\"/\t/;
  $l =~ s/\\n/\t/;
  $l =~ s/Cov:/\t/;
  $l =~ s/Len:/\t/;
  $l =~ s/\"\]\;/\t/;
  
  
  
  my @line = split /\s+/, $l;
  shift(@line);
  if ($l =~ /->/){
    my $c=$line[0];
    my $d=$line[2];
    my $dirA=0;
    if ($line[3] =~ /^5/){
      $dirA=1;
    }else{
      $dirA=-1;
    }
    my $dirB=0;
    if ($line[3] =~ /to_5/){
      $dirB=-1;
    }else{
      $dirB=1;
    }
    @{$linkCounts{$c}{$d}}=($dirA, $dirB);
    $dirA=$dirA*-1;
    $dirB=$dirB*-1;
    @{$linkCounts{$d}{$c}}=($dirB, $dirA);
  }else{
    my $n = shift @line;
    shift @line;
    @{$ContigInfo{$n}}=@line;
  }
}
close IN2;


my $infile3=shift;
open (IN3, $infile3) or die ("Couldn't open file: $infile3\n");

my %annotation=();
my $name1="";
my @coord="";
my %seen=();
my %stich=();

while (my $li=<IN3>){
  next unless ($li =~ /\w+|^\//);
  chomp $li;
  if ($li =~ /^\w+|^\//){
    $name1=$li;
  }elsif ($li =~ /^\t\d+/){
    my @t = split /\t/, $li;
    shift @t;
    @coord=@t;
  }else{
    if ($li !~ /StartPosQ/) {
    	my @t = split /\t+/, $li;
    	#my @Acc = split /\_AccNo_/, $t[1];
    	#@{$annotation{$Acc[1]}{$name1}}=(@coord, $t[3], $t[4], $t[5]);
    	#New Part of Code introduced############
   	my $ind = grep {$t[$_] =~ /\_AccNo_/} (0..scalar(@t)-1);
    	#Helpful for Debug
	#if (not $ind) {
	#	print "$li\n";
    	#}	 
    	my @Acc = split /\_AccNo_/, $t[$ind];
    	@{$annotation{$Acc[1]}{$name1}}=(@coord, $t[$ind+2], $t[$ind+3], $t[$ind+4]);
    	########New Code End###############
	my @k = keys (%{$annotation{$Acc[1]}});
    	if (scalar(@k) > 1){
    	  foreach my $n (@k){
    	    if ($linkCounts{$n}){
    	      foreach my $t (@k){
    	        if ($linkCounts{$n}{$t}){
    	          my @a=($n, $t, $Acc[1], $Acc[0], @{$annotation{$Acc[1]}{$n}}, @{$annotation{$Acc[1]}{$t}}, $linkCounts{$n}{$t}[0], $linkCounts{$n}{$t}[1], $conLen{$n}, $conLen{$t});
    	          push @{$stich{$Acc[1]}}, [@a];
    	          #print "Contigs $n and $t have the same annotation\n\t$Acc[1], $Acc[0]\n";
    	          #print "\tCoordinates are $annotation{$Acc[1]}{$n}[0] - $annotation{$Acc[1]}{$n}[1]";
    	          #print "\twith DB: hits from $annotation{$Acc[1]}{$n}[3] - $annotation{$Acc[1]}{$n}[4] and PerId $annotation{$Acc[1]}{$n}[2]\n";
    	          #print "\tCoordinates are $annotation{$Acc[1]}{$t}[0] - $annotation{$Acc[1]}{$t}[1]";
    	          #print "\twith DB: hits from $annotation{$Acc[1]}{$t}[3] - $annotation{$Acc[1]}{$t}[4] and PerId $annotation{$Acc[1]}{$t}[2]\n";
    	          #print "\tDir join is $linkCounts{$n}{$t}[0] and $linkCounts{$n}{$t}[1]\n";
    	          #print "\tThe length is $conLen{$n} and $conLen{$t}\n";
    	          delete $linkCounts{$n}{$t};
    	          delete $linkCounts{$t}{$n};
    	        }
    	      }
    	    }
    	  }
    	}
     }
  }
}

foreach my $k (keys (%stich)){
  #print "$k\n";
  &proc_comp(\@{$stich{$k}});
  #foreach my $l (@{$stich{$k}}){
  #  my $pl = join ",", @{$l};
  #  print "\t$pl\n";
  #}
}

foreach my $h (keys (%contigs)){
  print ">$h\n$contigs{$h}\n" unless $seen{$h};
}

############### Subroutines ##########

sub proc_comp{
  #### Assumption, for any given genbank entry / hit only 1 protein present in the sample.
  my $start=0;
  my $stop=0;
  my $len=0;
  my $contigs="";
  my $hcov=shift(@_);
  

  foreach my $t (@{$hcov}){
    my $cut=0;
    # for each entry with hit to a db, search if the hit overlap if the same, get the largest annotation, if overlap make an assembly. Assign that
    # If already exists, check if overlap if does, assemble if contained, ignore.
    # In @s, ContigA, ContigB, FromDB_A, ToDB_A, FromDB_B, ToDB_B,     DirA,     DirB
    #           0        1        2         3         4       5          6         7
    my @s = ($t->[0], $t->[1], $t->[7], $t->[8], $t->[12], $t->[13], $t->[14], $t->[15]);
    my ($c_a, $c_b, $len_f, $con) = &get_assembly($s[2], $s[3], $s[4], $s[5], $s[0], $s[1]);
    ## if worth make the assembly, if contained discard the contig to %seen;
    
    my $new_con="";
    
    if ($con eq $s[0]){
      #print "Contig is $s[0], discarding $s[1]\n";
      $seen{$s[1]}=1;
      $new_con=$contigs{$s[0]};
    }elsif ($con eq $s[1]){
      #print "Contig is $s[1], discarding $s[0]\n";
      $seen{$s[0]}=1;
      $new_con=$contigs{$s[1]};
    }else{
      #print "Assembling $s[0] and $s[1], the direction is $s[6] and $s[7], the database pos are $s[2]-$s[3] and $s[4]-$s[5] can I just stich?\n";
      my $seqA = $contigs{$s[0]};
      my $seqB = $contigs{$s[1]};
      my ($min_s, $max_s) = sort {$a <=> $b} ($s[2], $s[4]);
      my ($min_z, $max_z) = sort {$a <=> $b} ($s[3], $s[5]);

      $cut = $min_z-$max_s;
      
      if ($s[6] == -1){
	my $rev = reverse($seqA);
	$rev =~ tr/ATCGatcg/TAGCtagc/;
	$seqA=$rev;
      }
      if ($s[7] == -1){
	my $rev = reverse($seqB);
	$rev =~ tr/ATCGatcg/TAGCtagc/;
	$seqB=$rev;
      }
      if ($cut > 0){
	my $lc = $cut*3;
	if ($max_s == $s[2]){
	  $seqA = substr ($seqA, $lc);
	}elsif ($max_s == $s[4]){
	  $seqB = substr ($seqB, $lc);
	}
      }
      $new_con=$seqA.$seqB;
    }
	

    

    if (($start==0) && ($stop ==0)){
      $start=$c_a;
      $stop=$c_b;
      $len=$len_f;
      $contigs=$con;
      #save the assembly
      print ">StichContigs_$con\n$new_con\n";
      #print ">StichContigs_$con\n";

      #print "Going to assemble $con\n";
      #print "Assembling $s[0] and $s[1], the direction is $s[6] and $s[7], the database pos are $s[2]-$s[3] and $s[4]-$s[5] and I am taking $cut aa out of the second seq?\n" if (($con =~/\s+/) && ($cut>0));
    }else{
      my ($start_new, $stop_new, $len_new, $contigs_new) = &get_assembly($start, $stop, $c_a, $c_b, $contigs, $con);
      ## again if worth print the assembly, else discard the contigs to %seen;
      if ($len_new > $len){
	#print "Is worth assembling $con\n";
	print ">StichContigs_$con\n$new_con\n";
	#print ">StichContigs_$con\n";
	#print "Assembling $s[0] and $s[1], the direction is $s[6] and $s[7], the database pos are $s[2]-$s[3] and $s[4]-$s[5] and I am taking $cut aa out of the second seq?\n" if (($con =~/\s+/) && ($cut>0));
      }else{
	#print "I am discarding $con\n";
	my @p = split /\s+/, $con;
	foreach my $o (@p){
	  $seen{$o}=1;
	}
      }
      ($start, $stop, $len, $contigs)=($start_new, $stop_new, $len_new, $contigs_new)
      
    }
  }
  #print "Assembly goes from $start, $stop, and is made of $contigs\n";
}
      
    
    

sub get_assembly{
  my ($start_a, $stop_a, $start_b, $stop_b, $con_a, $con_b) = @_;
  my ($min_s, $max_s) = sort {$a <=> $b} ($start_a, $start_b);
  my ($min_z, $max_z) = sort {$a <=> $b} ($stop_a, $stop_b);
  my $lenA = $stop_a - $start_a;
  my $lenB = $stop_b - $start_b;
  my $len_t= $max_z - $min_s;
  my $c_a;
  my $c_b;
  my $len_f;
  my $con="";
  if (($lenA / $len_t) > 0.8){
    $c_a=$start_a;
    $c_b=$stop_a;
    $len_f=$lenA;
    $con=$con_a;
  }elsif (($lenB / $len_t) > 0.8){
    $c_a=$start_b;
    $c_b=$stop_b;
    $len_f=$lenB;
    $con=$con_b;
  }else{
    $c_a=$min_s;
    $c_b=$max_z;
    $len_f=$len_t;
    $con="$con_a $con_b";
  }
  return ($c_a, $c_b, $len_f, $con);
}
  



