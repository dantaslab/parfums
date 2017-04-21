#!/usr/bin/perl -w

# Written by Alejandro Reyes

use strict;

if (@ARGV != 3) {
	die "\nUsage: Parse-Cog_Blast.pl [Cog_DB_Names] [Contig File] [Blast Outfile]  > <outfile>\n\n";
}

#Reads the Database Names, Gets the length and the taxonomical classification.

my $INfile = shift @ARGV;

open (IN, "<$INfile") or die ("Couldn't open file: $INfile\n");

my %DBnames = ();

#print "Begining to process DB\n";
while (my $lines=<IN>){			### NEED a tab delimited file: Accession \t Anything \t COG Categories;
	chomp $lines;
	if ($lines =~ /^\S+/){
		my @elements = split (/\t/, $lines);
		$DBnames{$elements[0]}=$elements[-1];
	}
}
close IN;

my $genome_file = shift;

open (IN, $genome_file) or die ("Couldn't open file $genome_file\n");

my %genomeSize=();
my $name="";
my $seq="";

while (my $line =<IN>){
  chomp $line;
  
  if ($line =~ /^>/){
    $line =~ s/>//;
    $line =~ s/\s+$//;
    if ($name ne ""){
      $genomeSize{$name}=length($seq);
    }
    $name = $line;
    $seq="";
  }else{
    $seq.=$line;
  }
}
$genomeSize{$name}=length($seq);
my $c = scalar(keys %genomeSize);
print "$c\n";
close IN;






#Reads the Blast Output in Table format

my $file = shift @ARGV;

open (IN, "<$file") or die ("Couldn't open file: $file\n");

# Loops in each line

print "Query\tDB_Name\tE-val\tScore\tPercId\tStartPosQ\tEndPosQ\tStartDB\tStopDB\tPercent_query\tPercent_DB\tAnnotation\n"; 

my @hits=();


while (my $lines=<IN>){
	chomp $lines;
	#print "Parsing Blast Output file\n";	
	if ($lines !~ /^$/){
		#Each element from the file is taken:
		# 0) Read Name
		# 1) DB Accession
		# 2) Percent ID
		# 3) Length Alignment
		# 4) Mismatches
		# 5) Gaps
		# 6) start coordinate in the query sequence
		# 7) end coordinate in the query sequence
		# 8) start coordinate in the db sequence
		# 9) end coordinate in the db sequence
		# 10) E-value
		# 11) Score
		
		my @element = split (/\t/,$lines);
		if (scalar(@element) != 12){
			die ("La entrada \n$lines\n, no tiene 12 elementos\n");
		}
		print "@element\n";
		my @k = &get_dir(\@element);
		print "@k\n";
		@element = @k;

		if(($hits[-1]) && ($element[0] ne $hits[-1][0])){
			&process_hit(\@hits);
			@hits=();
		}
		push @hits, [@element];
	}
}
&process_hit(\@hits);

################ Subroutines ##########

sub get_dir{
  my $array=shift(@_);
  
  my $dir="";
  if ((($array->[6] > $array->[7]) && ($array->[8] > $array->[9])) or (($array->[6] < $array->[7]) && ($array->[8] < $array->[9]))){
    $dir=1;
  }else{
    $dir=-1;
  }
  push (@{$array}, $dir);
  ($array->[6], $array->[7])=sort {$a <=> $b} ($array->[6], $array->[7]);
  ($array->[8], $array->[9])=sort {$a <=> $b} ($array->[8], $array->[9]);
  return @{$array};
}


sub process_hit {
  my $hits=shift(@_);
  my $last="";
  my @same=();
  my @new=();
    
  foreach my $k(@{$hits}){
    if ($k->[1] ne $last){
      if ($last ne ""){
	if (scalar(@same) > 1){
	  my $s=&proc_pair(\@same);
	  push (@new, $s);
	}else{
	  push (@new, \@{$same[0]});
	}
      }
      $last=$k->[1];
      @same=();
    }
    push (@same, $k);
  }
  &proc_comp(\@new);
}

sub proc_pair{
  my $array2D = shift(@_);
  
  # In this case both the query and subject are the same.
  my @q=();
  my @db=();
  my @return=();
  my $score=0;
  
  for (my $i=0; $i<@{$array2D}; $i++){
    my ($qstart, $qstop) = ($array2D->[$i][6], $array2D->[$i][7]);
    my ($dbstart, $dbstop)= ($array2D->[$i][8], $array2D->[$i][9]);
    
    # If the hits are overlaping in the query and overlaping or not inclusive on the subject, treat it as a big hit.
    # If the hits are not overlaping but gap is <10% of covered region, join them
    # Else, pick the largest of the hits and return it;
   my $overlap=0;
   for (my $m=$qstart; $m<=$qstop; $m++){
	$overlap++ if ($q[$m]);
   }
   
    next if ($overlap > (($qstop - $qstart)*0.6));
    
    for (my $m=$qstart; $m<=$qstop; $m++){
      $q[$m]=$array2D->[$i][2]/100 unless $q[$m];
    }
    for (my $m=$dbstart; $m<=$dbstop; $m++){
      $db[$m]=1;
    }
    $score+=$array2D->[$i][11];
  }
  
  my $q_len=0;
  my $q_st=0;
  my $db_len=0;
  my $db_st=0;
  my $q_g=0;
  my $switch=0;
  my $perId=0;
  
  for (my $i=0; $i<$#q; $i++){
    if ($q[$i]){
      $q_len++;
      $perId+=$q[$i];
      $q_st=$i if ($switch ==0);
      $switch=1;
    }else{
      $q_g++ if ($switch == 1);
      $q_len++ if ($switch == 1);
    }
  }
  
  $switch=0;
  for (my $i=0; $i<$#db; $i++){
    if ($db[$i]){
      $db_len++;
      $db_st=$i if ($switch ==0);
      $switch=1;
    }else{
      $db_len++ if ($switch == 1);
    }
  }

  if (($q_g/$q_len) < 0.1){
    my $pid=($perId/$q_len)*100;
    $array2D->[0][2]=sprintf ("%.2f", $pid);              #Calculate the PerID
    $array2D->[0][6]=$q_st;
    $array2D->[0][7]=$#q;
    $array2D->[0][3]=$#q-$q_st+1;
    $array2D->[0][8]=$db_st;
    $array2D->[0][9]=$#db;
    $array2D->[0][11]=$score;
  }
  return $array2D->[0];  
  
  # We are interested in columns 6, 7, 8, 9.
  # I need PercID, E-value, score.
  # For Score: Sum the scores
  # For PerID: PerID * lengthID on each and sum the divide the total length
  # For e-value, take the lowest
}

sub proc_comp{
  my $TwoDarray=shift(@_);
  
  foreach my $array (@{$TwoDarray}){
    die ("There is no genome size for $array->[0]\n") unless ($genomeSize{$array->[0]});
    my $lenq = $genomeSize{$array->[0]};	      
    my $perq = sprintf ("%.2f", (abs($array->[6] - $array->[7])+1)*100/$lenq);
	
    #my ($name, $dlen) = split /-/, $array->[1];
    my @annArray = split /-/, $array->[1];
    my $dlen = pop(@annArray);
    my $name = join ('-', @annArray);
    next unless ($DBnames{$name});
    
    my $posL="";

    ## Pos is 3' if the 3' of the sequence could be contained within the contig: i.e if the hit to DB starts in pos 20, and the hit in query starts in pos 40, means that even if the blast hit is not on the 3' end of the db, the region may as well be contained.
    if ($array->[-1] == 1){
      $posL.= "3'" if ($array->[6] >= $array->[8]*3);
      $posL.= "5'" if (($lenq-$array->[7]) >= ($dlen-$array->[9])*3);
    }else {
      $posL.= "3'" if (($lenq-$array->[7]) >= $array->[8]*3);
      $posL.= "5'" if ($array->[6] >= ($dlen-$array->[9])*3);
    }

    my $perd = sprintf ("%.2f", ($array->[9] - $array->[8] + 1) *100 / $dlen);

    $posL = "3'5'" if ($perd >= 90);
    $posL="internal" if $posL eq "";
    
    print "$array->[0]\t$array->[1]\t$array->[10]\t$array->[11]\t$array->[2]\t$array->[6]\t$array->[7]\t$array->[8]\t$array->[9]\t$perq\t$perd\t$posL\_$DBnames{$name}\_AccNo_$array->[1]\n";
  }
}


	
	
