#!/usr/bin/perl -w
use strict;

use FindBin qw($Bin);
use File::Basename;
use File::Spec::Functions qw(catdir);

use lib catdir(dirname($Bin), '/lib');
use Getopt::Long;
use Capture::Tiny ':all';
use IO::Compress::Gzip qw(gzip $GzipError);
use Log::Log4perl;
use Cwd;


die ("Usage: FR-Hit_parse.pl <Output_FR-Hit> <ContigFile> <prefix> <OUTFILE>\n") unless (scalar(@ARGV) == 4);

use constant WAIT_TIME => '60' ; # 60 seconds = 1 minute
use constant NO_CHANGE_IN_FILE_FOR_X_SECONDS => '30';



my $file = shift @ARGV;
my $fna = shift @ARGV;
my $prefix = shift @ARGV;
my $output_file = shift @ARGV;
my $dirname = dirname($file);
my @seqs=();
my %links=();
my $linkCount = 0;
my $linkThreshold = 5;
my %linkNum=();

open (IN, "<$file") or die ("Couldn't open file: $file\n");
open (FNA, "<$fna") or die ("Couldn't open file: $fna\n");

my %len=();
my %refs=();
my $name="";
while (my $k = <FNA>){
  chomp($k);
  if ($k =~ /^>(\S+)\s*$/){
    print "Name is >$name< and size is $len{$name}\n";
    $name=$1;
  }else{
    $refs{$name}.=$k;
    $len{$name}=length($refs{$name});
  }
}


my $old="";
my $n="";

while (my $l = <IN>){
  chomp ($l);
  my @array = split /\t/, $l;
  $n = $array[0];
  $n =~ s/\_\d$//;

  $array[7] =~ s/\%$//;
  $array[1] =~ s/nt$//;
  next if ($array[7] < 95);

  if ($n ne $old){
    if ($old ne ""){
      &process_array(\@seqs, $n);
    }
    $old = $n;
    @seqs=();
  }
  push (@seqs, [@array]);
}
close IN;
&process_array(\@seqs, $n);

my %contigs=();
my %used=();
my %skip=();
my $counter=1;
my $count2=1;

#die ("End of Test\n");

foreach my $l (keys(%links)){
  foreach my $t (keys(%{$links{$l}})){
    if ($links{$l}{$t} >= $linkThreshold){
    	$linkCount++;
    	$linkNum{$l}{$t} = $linkCount;
    }
  }
}
if ($linkCount > 39){
   $linkThreshold = 10;
   $linkCount = 0;
   foreach my $l (keys(%links)){
   	foreach my $t (keys(%{$links{$l}})){
    	if ($links{$l}{$t} >= $linkThreshold){
    		$linkCount++;
    		$linkNum{$l}{$t} = $linkCount;
    	}
   	}
   }
   if ($linkCount > 39){
   		$linkThreshold = 20;
   		$linkCount = 0;
   		foreach my $l (keys(%links)){
   			foreach my $t (keys(%{$links{$l}})){
    			if ($links{$l}{$t} >= $linkThreshold){
    				$linkCount++;
    				$linkNum{$l}{$t} = $linkCount;
    			}
   			}
   		}
   			if ($linkCount > 39){
   				$linkThreshold = 40;
   				$linkCount = 0;
   				foreach my $l (keys(%links)){
   					foreach my $t (keys(%{$links{$l}})){
    					if ($links{$l}{$t} >= $linkThreshold){
    						$linkCount++;
    						$linkNum{$l}{$t} = $linkCount;
    					}
   					}
   				}
   				if ($linkCount > 39){
   					$linkThreshold = 60;
   					$linkCount = 0;
   					foreach my $l (keys(%links)){
   						foreach my $t (keys(%{$links{$l}})){
    						if ($links{$l}{$t} >= $linkThreshold){
    						$linkCount++;
    						$linkNum{$l}{$t} = $linkCount;
    						}
   						}
   					}
   					if ($linkCount > 39){
   						$linkThreshold = 80;
   						$linkCount = 0;
   						foreach my $l (keys(%links)){
   							foreach my $t (keys(%{$links{$l}})){
    							if ($links{$l}{$t} >= $linkThreshold){
    							$linkCount++;
    							$linkNum{$l}{$t} = $linkCount;
    							}
   							}
   						}
   						if ($linkCount > 39){
   							$linkThreshold = 100;
   							$linkCount = 0;
   							foreach my $l (keys(%links)){
   								foreach my $t (keys(%{$links{$l}})){
    								if ($links{$l}{$t} >= $linkThreshold){
    								$linkCount++;
    								$linkNum{$l}{$t} = $linkCount;
    								}
   								}
   							}
   						}
   					}
   				}
   			}
   }
}
foreach my $l (keys(%links)){
  foreach my $t (keys(%{$links{$l}})){
    if ($links{$l}{$t} >= $linkThreshold){
      print "Evaluating link between $l and $t that has $links{$l}{$t}\n";
      print "Link number $linkNum{$l}{$t} has a count of $links{$l}{$t} reads. A total of $linkCount links are being evaluated and the threshold for evaluation is $linkThreshold reads\n";
      my @gens = split /\!/, $l;
      my @pos = split /\!/, $t;
      my $sub1="";
      my $sub2="";
      my $return_file="";

      $sub1 = substr $refs{$gens[0]}, 0, 100 if ($pos[0] =~/start/);
      $sub1 = substr $refs{$gens[0]}, -100 if ($pos[0] =~/stop/);
      $sub2 = substr $refs{$gens[1]}, 0, 100 if ($pos[1] =~/start/);
      $sub2 = substr $refs{$gens[1]}, -100 if ($pos[1] =~/stop/);
      $used{"$gens[0].$counter"}=$refs{$gens[0]};
      $used{"$gens[1].$counter"}=$refs{$gens[1]};
      $skip{$gens[0]}=1;
      $skip{$gens[1]}=1;

      open (OUT, ">Link_ToMap.fna") or die ("Couldn't open file Link_ToMap.fna\n");
      print OUT ">$pos[0]\_$gens[0]\n$sub1\n>$pos[1]\_$gens[1]\n$sub2\n";
      close OUT;
      system("fr-hit -d Link_ToMap.fna -a $dirname/$prefix.fasta.clean -o Link_Map.txt -g 1 -q 50 -c 95");
      system("perl /opt/apps/labs/gdlab/software/parfums/1.0/src/FR-Hit_cleanLink.pl Link_Map.txt Link_ToGet_Seqs.txt");
      system("perl /opt/apps/labs/gdlab/software/parfums/1.0/src/GetSequences.pl Link_ToGet_Seqs.txt $dirname/$prefix.fasta.clean >> Link_ToMap.fna");
      system("cd-hit-est -i Link_ToMap.fna -o Link_OutCdHit.fna -G 0 -aS 0.99 -g 1 -r 1 -c 0.9");
      my $len=60;
      system("phrap -minmatch 10 -maxmatch 30 -bandwidth 0 -minscore 15 Link_OutCdHit.fna");
      system("fr-hit -d Link_OutCdHit.fna.contigs -o Link_Map.txt -a $dirname/$prefix.fasta.clean -m 30");
      system("perl /opt/apps/labs/gdlab/software/parfums/1.0/src/FR-Hit_cleanChimera.pl Link_Map.txt Link_OutCdHit.fna.contigs > Link_NoChimera.fna");
      my @contig =`cat Link_NoChimera.fna`;
      my $na = "";
      my $contigsRound=0;
      my $notPass="";

      foreach my $s (@contig){
      chomp($s);
      if ($s =~ /^>/){
	  if ($na ne ""){
	    if (length($contigs{$na}) < $len){
	      $notPass.=">$na\n$contigs{$na}\n";
	      delete $contigs{$na};
	    }else{
	      $contigsRound++;
	    }
	    #print "Contig is \n$na\n$contigs{$na}\n";
	  }
	  $na=$s;
	  $na.=".$counter";
	  $counter++;
	}else{
	  $contigs{$na}.=$s;
	}
      }
      if (length($contigs{$na}) < $len){
	$notPass.=">$na\n$contigs{$na}\n";
	delete $contigs{$na};
      }else{
	$contigsRound++;
      }
      #print "Contig is \n$na\n$contigs{$na}\n";

      print "Testing link No $count2 of $linkCount total links\n";
      $count2++;


      #Let's do some cleaning

      #die "No Contigs pass the threshold, len was $len, contigs are:\n$notPass\n" unless ($contigsRound > 0);
      #system("rm -r Link* SGE");
    }
  }
}

open (OUT, ">$output_file") or die ("Couldn't open file $output_file\n");

foreach my $o (keys(%refs)){
  print OUT ">$o\n$refs{$o}\n" unless ($skip{$o});
}
foreach my $p (keys(%used)){
  print OUT ">$p\n$used{$p}\n" unless length($used{$p}) < 30;
}
foreach my $u (keys(%contigs)){
  print OUT "$u\n$contigs{$u}\n" if ($contigs{$u} && length($contigs{$u}) > 30);
}
close OUT;


########### subroutines #######

sub process_array{
  my @arr = @{$_[0]};
  my $name = $_[1];
  my %seqs=();
  my %gens=();

  for (my $i=0; $i<@arr; $i++){
    push (@{$seqs{$arr[$i][8]}{$arr[$i][0]}}, [@{$arr[$i]}]);
    $gens{$arr[$i][0]}{$arr[$i][8]}++;
  }
  my @ke = keys(%gens);
  return unless scalar(@ke ==2);
  my @geA = keys (%{$gens{$ke[0]}});
  my @geB = keys (%{$gens{$ke[1]}});
  return unless (scalar(@geA)==1 && scalar(@geB)==1 && $geA[0] ne $geB[0]);


  my $switch1=0;
  my $l=0;

  my @a1=@{${$seqs{$geA[0]}{$ke[0]}}[0]};
  my @a2=@{${$seqs{$geB[0]}{$ke[1]}}[0]};

  die ("No existe len de >$geA[0]<\n") unless ($len{$geA[0]});
  die ("No existe len de >$geB[0]<\n") unless ($len{$geB[0]});


  my $link="";
  if ($a1[6] eq "+"){
    if($a1[-1] > ($len{$geA[0]}-100)){
      $link="stop_1";
    }else{
      $link="weird1";
    }
  }else{
    if ($a1[-2] < 100){
      $link="start_1";
    }else{
      $link="weird_1";
    }
  }
  if ($a2[6] eq "+"){
    if($a2[-1] > ($len{$geB[0]}-100)){
      $link.="!stop_2";
    }else{
      $link.="!weird2";
    }
  }else{
    if ($a2[-2] < 100){
      $link.="!start_2";
    }else{
      $link.="!weird_2";
    }
  }

  #print "@a1\n@a2\n$link\nWhere\n" if ($geA[0] eq "Clean_seq_964_8" && $geB[1] eq "Clean_seq_169_11" && $link =~ /stop_1-start_2/);

  if ($link !~ /weird/){
    if ($link eq "start_1!start_2" or $link eq "stop_1!stop_2" ){
      if ($links{"$geB[0]!$geA[0]"}{$link}){
	$links{"$geB[0]!$geA[0]"}{$link}++;
      }else{
	$links{"$geA[0]!$geB[0]"}{$link}++;
      }
    }elsif ($link eq "start_1!stop_2"){
      if ($links{"$geB[0]!$geA[0]"}{"stop_1!start_2"}){
	$links{"$geB[0]!$geA[0]"}{"stop_1!start_2"}++;
      }else{
	$links{"$geA[0]!$geB[0]"}{$link}++;
      }
    }elsif ($link eq "stop_1!start_2"){
      if ($links{"$geB[0]!$geA[0]"}{"start_1!stop_2"}){
	$links{"$geB[0]!$geA[0]"}{"start_1!stop_2"}++;
      }else{
	$links{"$geA[0]!$geB[0]"}{$link}++;
      }
    }
    #print "$ke[0]\t$geA[0]\t$seqs{$geA[0]}{$ke[0]}[0][-2]\t$seqs{$geB[0]}{$ke[1]}[0][-1]\n";
    #print "Link is $link\n";
  }

  return;
}


################### Subroutines #####################

sub wait_for_file{
  my($fileName) = @_;

  # Using constants WAIT_TIME and NO_CHANGE_IN_FILE_FOR_X_SECONDS

  print "Wait_for_file subroutine invoked at ".localtime().", waiting for file $fileName\n";

  # Check for the finished file every WAIT_TIME minutes
  until (-e "$fileName") {
    print "$fileName does not exist yet at ".localtime()."\n";
    sleep WAIT_TIME;
  }

  print "File $fileName found at ".localtime()."\n";

  # make sure the file is no longer being modified (changing size)
  my$time = 0;
  my$size = -s $fileName;

  until ($time == NO_CHANGE_IN_FILE_FOR_X_SECONDS) {
    if ($size == -s $fileName) {
      $time++;
      sleep 1;
      if ($time%5 == 0) {
	print "No change in file size for $time seconds\n";
      }
    } else {
      $time = 0;
      $size = -s $fileName;
      print "The file size changed, sleeping for 1 minute\t";
      sleep 60;
      print"Waking up, try again\n";
    }
  }
  print "File $fileName exists and hasn't been modifed for at least 1 minute at time ".localtime()."\n\n";
  return $fileName;
}
