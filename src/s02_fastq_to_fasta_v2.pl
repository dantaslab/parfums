#!/usr/bin/perl -w

#
# Created by: Manish Boolchandani
# Creation Date: 2016-05-19
# Version: 1.0
# Description: Converts fastq to fasta format
#
#Introduce random sampling for reads more than 4,000,000
#
#

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

my $MAX_READ = 4000000;
my $MIN_R = 100; my $MAX_R = 900;
#Global variables to store command line input
my ($inputDir, $bcfile, $sackedIDs, $help);
my @uncheckedIDs;
#START OF THE PROGRAM
&parse_command_line();

if ($sackedIDs) {
	@uncheckedIDs = split(/\,/, $sackedIDs);
}
#Initialize log files
my $dirname = dirname($Bin);
my $log_conf = "$dirname/config_files/log4perl.conf";
Log::Log4perl::init_once($log_conf);
my $logger = Log::Log4perl->get_logger("Fastq_To_Fasta");
$logger -> info("Converting fastq to fasta format");

my $currentDir = getcwd;

#$inputDir = "$currentDir/parfum_run2";
my $fastqDir = "$inputDir/PreprocessedFiles";

    if(-d $fastqDir){
        #my $fasta_dir = "$inputDir/FastaFiles";
        #if(!-d $fasta_dir){
        #    mkdir "$inputDir/FastaFiles" or die("Error Creating direcotry: $!");
        #}
        open (my $bcIN, "< $bcfile") or $logger -> logdie("$!");
	$logger -> info("Maximum number of reads allowed: $MAX_READ");
	#$logger -> info("We will randomly subsample, if > $MAX_READ");

        while (my $line = <$bcIN>) {
            chomp($line);
            my @vals = split(/\s+/, $line);
            my $sampleID = $vals[1];
            my $out_filePath = "$inputDir/$sampleID";
            if (not -d $out_filePath) {
                mkdir "$out_filePath" or $logger -> logdie("Error creating directory: $!");
            }
            my $fastaFile = "$out_filePath/$sampleID.fasta";
	    my $fw_fastqFile = "$fastqDir/$sampleID\_FW.fq";
	    my $rc_fastqFile = "$fastqDir/$sampleID\_RV.fq";
	    my $mergeFile = "$fastqDir/$sampleID.fastq";
 
	    #$logger -> info("Checking Reads in sample: $sampleID");
	    my $flag = &check_readCount($sampleID, $fw_fastqFile, $rc_fastqFile, $mergeFile);
	    #$logger -> debug("Reading file: $sampleID.fastq");
	    if ($flag) {
            	my $fq_to_fasta_cmd = "$dirname/lib/seqtk-master/seqtk seq -a $mergeFile > $fastaFile";

            	my ($std_out,$std_err,$exit_status)= capture {
                	system("$fq_to_fasta_cmd");
            	};
            	if($exit_status){
                	$logger -> logdie("Error in executing seq-tk program:$std_err");
            	}
            	else {
                	$logger -> debug("$sampleID.fasta created");
                	gzip $fastaFile => "$fastaFile.gz" or $logger -> logdie("gzip failed:$GzipError");
            	}
	    }
        }
    }
    else {
        print "Error: Directory doesnt exists\n";
        exit;
    }

    $logger -> debug("Format Conversion done successfully");

    #END OF THE PROGRAM

sub check_readCount {
	my $sampleid = $_[0];
	my $fwFile = $_[1];
	my $rcFile = $_[2];
	my $outFile = $_[3];
	my $chk_fw_cmd = "awk \'NR%4 == 1 {print \$0;}\' $fwFile | wc -l";
	my $chk_rc_cmd = "awk \'NR%4 == 1 {print \$0;}\' $rcFile | wc -l";
	my $fCount = `$chk_fw_cmd`; chomp($fCount);
	chomp(my $rCount = `$chk_rc_cmd`);
	$logger ->info("SampleID => $sampleid\tFW_reads => $fCount\tRV_reads => $rCount\n");
	
	if (($fCount == 0) && ($rCount == 0)) {
		$logger ->info("No reads in $sampleid file: $fCount and $rCount");
		return 0;
	} 
	elsif ($fCount == $rCount) {
		$fCount = $fCount + 0;
		if ($fCount > $MAX_READ) {
			my $fwOut = "$fwFile.sample";
			my $rcOut = "$rcFile.sample";
			my $percent = sprintf ("%.2f", ($MAX_READ/$fCount));
			my $random = $MIN_R + int(rand($MAX_R - $MIN_R));
			my $fw_cmd = "$dirname/lib/seqtk-master/seqtk  sample -s$random $fwFile $percent > $fwOut";
			my $rc_cmd = "$dirname/lib/seqtk-master/seqtk  sample -s$random $rcFile $percent > $rcOut";
			$logger -> info("There exists $fCount reads in the given fasta file. Thus, approx. $percent percent of reads will be used in the remainder of assembly pipeline");
			$logger -> debug("Seqtk run for $fwFile");
			#$logger ->debug("$fw_cmd");
			#$logger ->debug("$rc_cmd");
			my ($fw_out,$fw_err,$fw_exit)= capture { system("$fw_cmd"); };
			if ($fw_exit) {
				$logger -> logdie("Could not complete FW file sampling: $fw_err");
			}
			my ($rc_out,$rc_err,$rc_exit)= capture { system("$rc_cmd"); };
			if ($rc_exit) {
				$logger -> logdie("Could not complete RV file sampling: $rc_err");
			}
			if ((-e $fwOut) && (-e $rcOut)) {
				&merge_fastq_files($fwOut, $rcOut, $outFile);
			}
			return 1;
		}
		else {
			&merge_fastq_files($fwFile, $rcFile, $outFile);
			return 1;
		}
	}
	else {
		$logger ->logdie("FW and RV file has unequal entries: $fwFile: $fCount and $rcFile: $rCount");
		return 0;
	}
}		

sub merge_fastq_files {
	my $fwFile = $_[0];
	my $rcFile = $_[1];
	my $outFile = $_[2];
	open (my $fwIn, "< $fwFile") or die "Cannot open $fwFile $!";
	open (my $rcIn, "< $rcFile") or die "Cannot open $rcFile $!";
	open (my $out, "> $outFile") or die "Cannot open $outFile $!";
	while(defined(my $fwline = <$fwIn>) && defined(my $rcline = <$rcIn>)) {
		chomp($fwline); chomp($rcline);
		my $fwID = $fwline; my $rcID = $rcline;
		if (($fwline =~ /^@/) && ($rcline =~ /^@/)) {
			my @seqID1 = split(/\#/, $fwline);
			my @seqID2 = split(/\#/, $rcline);
			if ($seqID1[0] eq $seqID2[0]) {
				#print "$fwline\t$rcline\n";	
				chomp(my $fwSeq = <$fwIn>); chomp(my $rcSeq = <$rcIn>);
				$fwline = <$fwIn>; $rcline = <$rcIn>;
				chomp(my $fwQual = <$fwIn>); chomp(my $rcQual = <$rcIn>);
				print $out "$fwID\n$fwSeq\n+\n$fwQual\n";
				print $out "$rcID\n$rcSeq\n+\n$rcQual\n";
			}
			else {
				$logger ->logdie ("Incorrect format of $fwFile and $rcFile; cannot match IDs $seqID1[0] and $seqID2[0]");
			}
		}		
	} 
}

sub get_Logfile_Name {
    return "$inputDir/parfums.log";
}

sub parse_command_line {
    my $result = undef;
    $result = GetOptions ( "dir=s" => \$inputDir,
		  "bcfile=s"  => \$bcfile,
		  "maxread=i" => \$MAX_READ,
		  "unchecked=s" => \$sackedIDs,
		  "help" => \$help
		  );
    &usage() if ($help);
}
