#!/usr/bin/perl -w
############################
#   Author: Manish Boolchandani
#   Date: 2015/02/19
#   Purpose: Parse through cross match and remove/trim
#   the reads with adapters and vectors
###########################

use strict;
use FindBin qw($Bin);
use File::Basename;
use File::Spec::Functions qw(catdir);

use lib catdir(dirname($Bin), '/lib');
use Getopt::Long;
use Log::Log4perl;

#Global variables for arguments
my ($bcFile, $inputDir, $sackedIDs);
my (%seqDel, %adapters, %vectors);

#START OF THE PROGRAM
#Subroutine to parse the command line
&parse_command_line();

#setup log file configuration path
my $dirname = dirname($Bin);
my $log_conf = "$dirname/config_files/log4perl.conf";
Log::Log4perl::init_once($log_conf);
my $logger = Log::Log4perl->get_logger("Adapter_Trim");
$logger -> info("Trimming of Adapter sequences started");

my @uncheckedIDs;
if ($sackedIDs) {
    @uncheckedIDs = split(/\,/, $sackedIDs);
}
#Subroutine to parse through cross-match output.
$logger -> debug("Parsing cross_match output");
open (my $bcIN, "< $bcFile") or die "$!";
chomp(my @BC = <$bcIN>);
foreach my $barcode (@BC) {
    my @seqID = split(/\s+/, $barcode);
    if (not grep {$_ eq $seqID[1]} @uncheckedIDs) {
        my $cmFile = "$inputDir/$seqID[1]/$seqID[1].cmoutput";
        my $fastaFile = "$inputDir/$seqID[1]/$seqID[1].fasta";
        if ((-e $cmFile) && (-e $fastaFile)) {
		my $cmd = "perl $dirname/src/CleanAdapter-Illumina_PE_mod.pl $fastaFile $cmFile";
		$logger -> debug("Clean Adaptors CMD: $cmd");
		my $out = `$cmd`;
        #    &read_CrossMatch($cmFile);
        #    my $seqkey = scalar(keys %seqDel);
        #    my $adapKey = scalar(keys %adapters);
        #    my $vecKey = scalar(keys %vectors);
        #    #$logger -> debug("SeqDel: $seqkey\nAdapters:$adapKey\nVectors:$vecKey\n");
        #    #Subroutine to trim adapters and vectors
        #    $logger -> debug("Function call to trim $seqID[1].fasta");
        #    &trim_Sequences($fastaFile);
        #    undef %seqDel; undef %adapters; undef %vectors;
	        if (-e "$fastaFile.clean") {
	       		my $seqC = `grep '^>' $fastaFile.clean | wc -l`;
	       		$logger -> debug("Read Count in $fastaFile.clean --> $seqC");
	      	 }
        }
        else {
            $logger -> logdie("Couldn't locate fasta file or CM output for $seqID[1]. Please check: $!");
        }
    }
}
exit 0;


#END OF THE PROGRAM

sub trim_Sequences {
    my $fastaFile = shift;
    open (my $fastaIN, "< $fastaFile") or $logger -> logdie("$!");
    open (my $OUT, "> $fastaFile.clean") or $logger -> logdie("$!");
    my ($seq1, $id1, $seq2, $id2, $comID);
    while(my $line = <$fastaIN>) {
        chomp($line);
        if ($line =~ /^>/) {
            $id1 = $line; $id1 =~ s/^>//;
            $comID = $id1;
            $comID =~ s/\_\d$//;
            chomp($seq1 = <$fastaIN>);
            chomp($id2 = <$fastaIN>); $id2 =~ s/^>//;
            chomp($seq2 = <$fastaIN>);
            $logger -> logdie("Forward and Reverse IDs didn't match: $comID") if ($id2 !~ /$comID/);
            my $numN = $seq1 =~ s/N/N/g;
            $numN += $seq2 =~ s/N/N/g;
            if ((exists $seqDel{$comID}) or ($numN > 6)) {
                next;
            }

            if (exists $vectors{$id1}) {
                my @tmpArr = @{$vectors{$id1}};
                $seq1 = substr($seq1, $tmpArr[0], $tmpArr[1] - $tmpArr[0], 'XXXX');
            }
            if (exists $vectors{$id2}) {
                my @tmpArr = @{$vectors{$id2}};
                $seq2 = substr($seq2, $tmpArr[0], $tmpArr[1] - $tmpArr[0], 'XXXX');
            }
            if (exists $adapters{$comID}) {
               $seq1 = substr($seq1, 0, $adapters{$comID});
               $seq2 = substr($seq2, 0, $adapters{$comID});
            }
            $seq1 =~ s/X+//; $seq2 =~ s/X+//;
            print $OUT ">$id1\n$seq1\n>$id2\n$seq2\n";
        }
    }
    my $seqID = basename($fastaFile);
    $logger -> debug("File Trimmed: $seqID.fasta");
    close($fastaIN); close($OUT);
}

sub read_CrossMatch {
    my $cmFile = shift;
    my (@CMInfo, @data);
    open (my $cmIN, "< $cmFile") or $logger -> logdie("$!");
    while(my $line = <$cmIN>) {
        if ($line =~ /^Maximal single base matches/) {
            while((defined($line = <$cmIN>)) && ($line !~ /matching entries/)) {
                if ($line =~ /^\s+\d+/) {
                    chomp($line);
                    push(@data, $line);
                }
            }
            push(@CMInfo, @data);
            @data=();
        }
    }
    my $count = scalar(@CMInfo);
    print "$count\n";
    foreach my $cmline (@CMInfo) {
        my @vals = split(/\s+/, $cmline);
        my $subID = $vals[5];
        $subID =~ s/\_\d$//;
        my $flankLen = $1 if ($vals[8] =~/\((\d+)\)/);
        if ($cmline =~ /Adapter/) {
            next unless (($vals[-3] !~ /\(/) && ($vals[-3] < 6));
            if ($vals[6] < 35) {
                $seqDel{$subID} = 1;
            }
            else {
                $adapters{$subID} = $vals[6] - 1;
            }
        }
        if ($cmline =~ /Vector/) {
            if ($vals[6] < 6) {
                $vectors{$vals[5]} = [ 0,$vals[7] ];
            }
            if ($flankLen <= 5) {
                $vectors{$vals[5]} = [ $vals[6]-1,$vals[7]+$flankLen ];
            }
            else {
                $seqDel{$subID} = 1;
            }
        }
    }
}

sub get_Logfile_Name {
    return "$inputDir/parfums.log";
}

sub parse_command_line {
    my $help;
    my $result = GetOptions ( "bcfile=s" => \$bcFile,
		  "dir=s"  => \$inputDir,
		  "help" => \$help,
		  "unchecked=s" => \$sackedIDs
        ) ;
          &usage() if ($help);
}
