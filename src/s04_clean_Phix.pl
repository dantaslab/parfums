#!/usr/bin/perl -w
use strict;

#####################################
#   Author: Manish Boolchandani
#   Date: 2015/02/27
#   Purpose: Clean PhiX
#   Output: ID.noPhiX.fasta files
#####################################

use strict;
use FindBin qw($Bin);
use File::Basename;
use File::Spec::Functions qw(catdir);

use lib catdir(dirname($Bin), '/lib');
use Getopt::Long;
use Log::Log4perl;
use PARFuMs_Subs;
use File::Temp qw(tempdir tempfile);

#Global variables
my ($bcFile, $inputDir, $phixPath, $sackedIDs);
my %frhit_out;
my @uncheckedIDs;
#START OF THE PROGRAM
&parse_command_line();

if ($sackedIDs) {
	@uncheckedIDs = split(/\,/, $sackedIDs);
}
#Setup the log configuration path
my $dirname = dirname($Bin);
my $log_conf = "$dirname/config_files/log4perl.conf";
Log::Log4perl::init_once($log_conf);
my $logger = Log::Log4perl->get_logger("PhiX_Clean");
$logger ->info("STEP 3: REMOVING PhiX and VECTOR SEQUENCES");

#Temp Directory for frhit run
my $template = "phiX_frhitRun_XXXX";
my $tempDir = tempdir($template, DIR => "$inputDir/tmp/", CLEANUP => 0);

#Change the permissions of temporary directory in order to make it accessible by everyone. 
if( -d $tempDir ) {
        $logger -> debug("$tempDir exists");
        chmod(0755, $tempDir) or $logger -> die("Unable to change the permissions of dir $tempDir");
} else {
        $logger -> debug("$tempDir does not exist");
}

#Subroutine to run fr-hit on all files
my $frhit_script = &generate_FRhit_Script();
PARFuMs_Subs::submit_jobArray($frhit_script, "phiX_frhit", $tempDir, $logger);

#Subroutine to clean the file based on fr-hit output
$logger -> info("frhit for PhiX completed");
my $idCount = scalar(keys %frhit_out);
$logger -> debug("$idCount IDs in frhit_out");
my $refH = PARFuMs_Subs::frhit_clean(\%frhit_out, $logger);
my %matchReads = %{$refH};

#Subroutine to remove match reads from fastafile
&clean_FastaFile(\%matchReads);

exit 0;

sub clean_FastaFile {
    my %matchR = %{$_[0]};
    foreach my $k (keys %matchR) {
        my $fastaFile = "$inputDir/$k/$k.fasta.clean";
        my $outFile = "$inputDir/$k/$k.noPhiX.fasta";
        if (not exists $matchR{$k}{"No_Match_Reads"}) {
            if (-e $fastaFile) {
                $logger -> debug("Removing PhiX vector instances: $k");
                open (my $in, "< $fastaFile") or die "$!";
                open (my $out, "> $outFile") or die "$!";
                while(my $line = <$in>) {
                    chomp($line);
                    if ($line =~ /^>/) {
                        my $id = $line;
                        $id =~ s/>//;
                        $id =~ s/\_\d$//;
                        if (not exists $matchR{$k}{$id}) {
                            chomp(my $seq = <$in>);
                            print $out ">$line\n$seq\n";
                        }
                    }
                }
		close($out)
            }
        }
        else {
            system("cp", "$fastaFile", "$outFile");
        }
	my $seqC = `grep '^>' $outFile | wc -l`;
	$logger -> debug("Sequence Count in $outFile --> $seqC");
    }
}

sub generate_FRhit_Script {
    open (my $bcIN, "< $bcFile") or die "$!";
    my $scriptfile = "$inputDir/tmp/cleanPhiX_FRhit_script.sh";
    open (my $out, "> $scriptfile") or die "$!";
    chomp(my @filenames = <$bcIN>);
    foreach my $line (@filenames) {
        my @seqID = split(/\s+/, $line);
	if (not grep {$seqID[1] eq $_}@uncheckedIDs) {
        	my $cleanfile = "$inputDir/$seqID[1]/$seqID[1].fasta.clean";
       		if (-e $cleanfile) {
            		my $outFile = "$inputDir/$seqID[1]/$seqID[1]_Map_PhiX.txt";
            		my $frhit_cmd = &get_FRhit_cmd($cleanfile,$phixPath,$outFile);
            		print $out "$frhit_cmd\n";
            		$frhit_out{$seqID[1]} = $outFile;
        	}
        	else {
            		$logger->info ("Could not locate $cleanfile: Please check error log");
        	}
	}
	else {
		$logger ->info("$seqID[1] is not included in analysis");
	}
    }
    return $scriptfile;

}

sub get_FRhit_cmd {
    my @param = @_;
    print "@param\n";
    my $frhit = "fr-hit -a FASTA_FILE_1 -d FASTA_FILE_2 -o OUTPUT_FILE -c 70 -m 40 -r 0";
    $frhit =~ s/FASTA_FILE_1/$param[0]/;
    $frhit =~ s/FASTA_FILE_2/$param[1]/;
    $frhit =~ s/OUTPUT_FILE/$param[2]/;
    return $frhit;
}

sub frhit_qsub {
    $logger -> info("Running fr-hit on cluster using qsub");
    my $filename = shift;
    my $jobName = "parfum_frhit";
    my $workDir = "tmp_$jobName";
    if (not(-d $workDir)) {
        mkdir $workDir;
    }
    chdir($workDir);
    my $cmd = "nq $filename | qsub -N $jobName -sync y";
    my $job_return = `$cmd`;
    #$logger -> info("$job_return");
    my $jobID = $1 if ($job_return =~ /^Your job-array (\d+)\./);
    #print "$jobID\n";
    $logger -> info("job $jobID completed");
}

sub get_Logfile_Name {
    return "$inputDir/parfums.log";
}

sub parse_command_line {
    my $help;
    my $result = GetOptions (
        "bcfile=s" => \$bcFile,
		"dir=s"  => \$inputDir,
		"phix=s" => \$phixPath,
		"unchecked=s" => \$sackedIDs,
        "help" => \$help
    ) ;
    &usage() if ($help);
}
