#!/usr/bin/perl -w

# Written by Alejandro Reyes

use strict;

if (@ARGV != 1) {
	die "\nUsage: Posparse-Cog_Blast.pl [Parsed Blast Outfile] > <outfile>\n\n";
}


my $file = shift @ARGV;

open (IN, "<$file") or die ("Couldn't open file: $file\n");

my @hits=();

my $trash=<IN>;
my $less_than_12 = 0;
while (my $lines=<IN>){
	chomp $lines;	
	if ($lines !~ /^$/){
		#Each element from the file is taken:
		# 0) Read Name
		# 1) DB Accession
		# 2) E-value
		# 3) Score
		# 4) Percent ID
		# 5) Start POs
		# 6) Stop Pos
	        # 7) Start DB
	        # 8) Stop DB
		# 9) PercId Query
		# 10) PercId Database
		# 11) Functional Annotation
		
		
		my @element = split (/\t/,$lines);
		if (scalar(@element) != 12){
			$less_than_12++;
			#if ($less_than_12 > 3) {
			#	die ("La entrada \n$lines\n, no tiene 12 elementos\n");
			#}
			next;
		}
		
		if(($hits[-1]) && ($element[0] ne $hits[-1][0])){
		  my @nhits =  sort {$b->[3] <=> $a->[3]} @hits;
		  &proc_comp(\@nhits);
		  @hits=();
		}
		push @hits, [@element];
	}
}
my @s_hits =  sort {$b->[3] <=> $a->[3]} @hits;
&proc_comp(\@s_hits);
#print "$less_than_12\n";



sub proc_comp{
  my %list=();
  my $hcov=shift(@_);

  foreach my $t (@{$hcov}){
    my @s = ($t->[-1], $t->[3], $t->[4], $t->[7], $t->[8], $t->[5], $t->[6]);
    my $switch=0;
    foreach my $l (keys (%list)){
      my ($m, $n, $score) = split /\_/, $l;
      next if (($t->[6] < $m) or ($t->[5] > $n));
      my $ov=0;
      $ov += ($m-$t->[5]) if (($m - $t->[5]) > 0);
      $ov += ($t->[6] - $n) if (($t->[6]-$n) > 0);
      my $overlap = $ov / ($t->[6] - $t->[5]+1);
      
      my ($min_s, $max_s) = sort {$a <=> $b} ($m, $t->[5]);
      my ($min_z, $max_z) = sort {$a <=> $b} ($n, $t->[6]);

     # print "Comparing query $t->[5] to $t->[6] while list is $m to $n and overlap is $overlap with $ov bp of overlap\n";

      if ($overlap <= 0.4){
	$switch=1;
	last if  ($t->[3] < ($score*0.9));
	#print "Score is $score and current score is $t->[3]\n";
	#my ($bs, $tr) = sort {$b <=> $a} ($score, $t->[3]);
	my $x = "$min_s\_$max_z\_$score";
	$list{$x}=$list{$l};
	push (@{$list{$x}}, [@s]);
	delete($list{$l}) unless ($x eq $l);
	last;
      }
    }
    my $x = "$t->[5]\_$t->[6]\_$t->[3]";
    push (@{$list{$x}}, [@s]) if $switch == 0;
  }

  my @k = keys (%list);
  my @t=();
  for (my $l=0; $l<@k; $l++){
    my @pl = split /\_/, $k[$l];
    push (@t, [$pl[0], $l]);
  }
  my @sorted = sort {$a->[0] <=> $b->[0]} @t;
  
  print "$hcov->[0][0]\n";
  for (my $i=0; $i<@sorted; $i++){
    my @ann = &parse_list(\@{$list{$k[$sorted[$i][1]]}});
    my $ann = join "\n\t\t", @ann;
    my ($start, $stop) = split /\_/, $k[$sorted[$i][1]];
    print "\t$start\t$stop\n\t\t$ann\n";
    
    
  }
}


sub parse_list{
	
  my $list=shift(@_);
#  my @sorted=  sort{$b->[1] <=> $a->[1]} @{$list};
  my @sorted= @{$list}; 
 
  my @new=();
  my %seen=();
  
  foreach (my $i=0; $i<@sorted; $i++){
    # print "Is looking at this $sorted[$i][0]\n";
    next if ($seen{$sorted[$i][0]});
 #   last if ($sorted[$i][1] < ($sorted[0][1]*0.9));
    my $pl = join "\t", @{$sorted[$i]};
    push @new, $pl;
    
    $seen{$sorted[$i][0]}=1;
  }
  return @new;
}
