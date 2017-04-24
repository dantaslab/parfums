#!/usr/bin/perl -w

#########################################################
#
#   Author: Manish Boolchandani
#   Date: 2016/05/20
#   Description: PARFuMs Wrapper program
#   Version: 1.1.0
#   Command Use: perl parfumWrapper.pl dir OutputFiles --step (1:3|[1-4]|0)
#
#########################################################

use strict;
use FindBin qw($Bin);
use File::Basename;
use File::Spec::Functions qw(catdir);

use lib "$Bin/lib";
use Getopt::Long;
use Log::Log4perl;
use Capture::Tiny ':all';
use PARFuMs_Subs;

#Global variables for command line
my $VERSION = "1.1.0";
my ($outputDir,$stepIn,$configPath);
my %config;
my @procStep;
my $MAXSTEP = 6;
my $IDs_Not_Passed = 0;
my $dirname = dirname(__FILE__);
#print $dirname;

#BOOLEAN VARIABLES
my $unmatch_stat = 0;

#START OF THE PROGRAM

#Subroutine to read command line inputs (config file path & step information)
&parse_command_line();

#Subroutine to read config file
&read_config_file();

#Initialize directories and files
&makeOutputDir();

#Initiating a Log file to track program's progress.
my $log_conf = "$Bin/config_files/log4perl.conf";
Log::Log4perl::init($log_conf);
my $logger = Log::Log4perl->get_logger("PARFuM_Wrapper");
$logger -> info("STARTED PARFUMS PIPELINE $VERSION");

print("PARFUMs started\n");

#Run Programs
foreach my $proc (@procStep) {
    if ($proc == 1) {
        my $processFlag = 0;
        my $cmd = "perl $dirname/src/s01_preprocessing_v5.pl --bcfile $config{'BC_file'} --fwd_file $config{'FW_file'} --rev_file $config{'RC_file'} --$config{'barcode_pos'} --dir $config{'DIR'}";
        if ($unmatch_stat) {
            $cmd = $cmd . " --unmatch";
        }
        my ($stdout,$stderr,$exit) = capture { system("$cmd"); };
        print "$stdout\n$stderr\n$exit\n";
	if ($exit) {
            &error_Notify($proc,$stderr);
        } else {
            $processFlag = 1;
        }

        if ($processFlag) {
            my $cmd1 = "perl $dirname/src/s02_fastq_to_fasta_v2.pl --dir $config{'DIR'} --bcfile $config{'BC_file'} --maxread $config{'MAX_READ'} --unchecked $IDs_Not_Passed";
	    $logger -> debug("$cmd1");
            my ($stdout1,$stderr1,$exit1) = capture { system("$cmd1"); };
            ($exit1) ? &error_Notify($proc,$stderr1) : $logger -> info("Fastq to Fasta conversion complete");
        }
        $logger -> info("STEP 1: PREPROCESSING COMPLETED");
    }
    elsif ($proc == 2) {
        my $processFlag = 0;
	my $suffix = "fasta";
	&check_output_files($suffix);
	$logger -> info("IDs not included at Step$proc: $IDs_Not_Passed");
        my $cmd = "perl $dirname/src/s03_make_adapter_run_CM.pl --bcfile $config{'BC_file'} --dir $config{'DIR'} --adapA $config{'adapA'} --adapB $config{'adapB'} --endVec1 $config{'endVec1'} --endVec1R $config{'endVec1R'} --endVec2 $config{'endVec2'} --endVec2R $config{'endVec2R'} --unchecked $IDs_Not_Passed";
        $logger -> debug("Step 2: $cmd");
	my ($stdout,$stderr,$exit) = capture { system("$cmd"); };
        #if ($exit) {
        #    &error_Notify($proc,$stderr);
        #} else {
        #    $processFlag = 1;
        #}
        #$logger ->info("IDs not Passed: $IDs_Not_Passed");
        #if ($processFlag) {
        #    $logger -> info("CrossMatch Run complete. Trimming Fasta files");
	#    my $cmd = "perl $dirname/src/trim_Adaptors_Vectors.pl --bcfile $config{'BC_file'} --dir $config{'DIR'} --unchecked $IDs_Not_Passed";
        #    my ($stdout1,$stderr1,$exit1) = capture { system("$cmd"); };
        #    #print "$stdout1\n";
        #    ($exit1) ? &error_Notify($proc,$stderr1) : $logger -> info("Adapters trimmed successfully");
        #}
	$logger -> info("STEP 2: ADAPTER TRIMMING COMPLETED");
    }
    elsif ($proc == 3) {
        my $processFlag = 1;
        my $suffix = "fasta.clean";
        &check_output_files($suffix);
	$logger -> info("IDs not included at Step$proc: $IDs_Not_Passed");
	my $cmd = "perl $dirname/src/s04_clean_Phix.pl --bcfile $config{'BC_file'} --dir $config{'DIR'} --phix $config{'phixDIR'} --unchecked $IDs_Not_Passed";
        
	$logger -> debug("cmd = $cmd");
	
	my ($stdout,$stderr,$exit) = capture {
            system("$cmd");
        };
        if ($exit) {
            &error_Notify($proc,$stderr);
        } else {
            $processFlag = 1;
        }
        if ($processFlag) {
            $logger -> info("PhiX cleaning complete. Going to RevVector");
            my $cmd = "perl $dirname/src/s05_remove_Vectors.pl --bcfile $config{'BC_file'} --dir $config{'DIR'} --vecfile $config{'vecfile'}";
            my ($stdout1,$stderr1,$exit1) = capture { system("$cmd"); };
            ($exit1) ? &error_Notify($proc,$stderr1) : $logger -> info("RevVector reads cleaned");
        }
	$logger -> info("STEP 3: PhiX and VECTOR SEQUENCES REMOVED");
    }
    elsif ($proc == 4) {
        my $suffix = "noVector.fasta";
        &check_output_files($suffix);
	$logger ->info("IDs not included at Step$proc: $IDs_Not_Passed");
        my $cmd = "perl $dirname/src/s06_velvet_assembly_v1.pl --bcfile $config{'BC_file'} --dir $config{'DIR'} --dblink $config{'DB_link'} --unchecked $IDs_Not_Passed";
        $logger -> info("$cmd");
	my ($stdout,$stderr,$exit) = capture {
            system("$cmd");
        };
        #print "$stdout\n";
        if ($exit) {
            &error_Notify($proc,$stderr);
        } else {
            $logger -> info("STEP 4: VELVET ASSEMBLY COMPLETED");
        }

    }
    elsif ($proc == 5) {
        my $suffix = "ForPhrap1.fasta";
        &check_output_files($suffix);
	$logger ->info("IDs not included at Step$proc: $IDs_Not_Passed");
        #$logger -> debug("IDs not Passed: $IDs_Not_Passed");
        my $cmd = "perl $dirname/src/s07_phrap_assembly_v2.pl --bcfile $config{'BC_file'} --dir $config{'DIR'} --dblink $config{'DB_link'} --unchecked $IDs_Not_Passed";
        $logger ->debug("$cmd");
	my ($stdout,$stderr,$exit) = capture {
            system("$cmd");
        };
        #print "$stdout\n";
        if ($exit) {
            &error_Notify($proc,$stderr);
        } else {
            $logger -> info("STEP 5: PHRAP ASSEMBLY COMPLETED");
        }

	#my $suffix = "last-contigs.fasta";
        #&check_output_files($suffix);
    }
    elsif ($proc == 6) {
	my $suffix = "last-contigs.fasta";
	&check_output_files($suffix);
	$logger ->info("IDs not included at Step$proc: $IDs_Not_Passed");
        #$logger -> debug("IDs not Passed: $IDs_Not_Passed");
        my $cmd;
	if ($config{'resANNOTATE'}) {
		$cmd = "perl $dirname/src/s08_annotate_v1.pl --bcfile $config{'BC_file'} --dir $config{'DIR'} --dblink $config{'DB_link'} --unchecked $IDs_Not_Passed --resannotate $config{'resANNOTATE'}";
	}
	else { 
		$cmd = "perl $dirname/src/s08_annotate_v1.pl --bcfile $config{'BC_file'} --dir $config{'DIR'} --dblink $config{'DB_link'} --unchecked $IDs_Not_Passed";
        }
	#print "$cmd\n";
	my ($stdout,$stderr,$exit) = capture {
            system("$cmd");
        };
        #print "$stdout\n";
        if ($exit) {
            &error_Notify($proc,$stderr);
        } else {
            $logger -> info("STEP 6: ANNOTATION COMPLETED");
        }
    }
    else {
        print "Step Count doesn't make sense in current program\n";
        $logger -> logdie("Step Count doesn't make sense in current program");
    }
}

$logger -> info("END OF PARFUMS PIPELINE");

exit 0;
#END OF THE PROGRAM

sub get_Logfile_Name {
    my $outputdir = $config{'DIR'};
    return "$outputdir/parfums.log";
}

sub check_output_files {
    my $suffix = shift;
    my $bcfile = $config{'BC_file'};
    my $dir = $config{'DIR'};
    open (my $in, "< $bcfile") or die "$!";
    while(my $line = <$in>) {
        chomp($line);
        my @seqID = split(/\s+/, $line);
        my $outfile = "$dir/$seqID[1]/$seqID[1].$suffix";
        if (!-e $outfile || -z $outfile) {
            $logger ->info("$seqID[1].$suffix is empty/not found. Halting its processing here");
            if ($IDs_Not_Passed) {
		my @currIDs = split(/\,/, $IDs_Not_Passed);
		if (not(grep(/^$seqID[1]$/, @currIDs))) {
                	$IDs_Not_Passed = join(',', ($IDs_Not_Passed, $seqID[1]));
		}
	    }
            else {
                $IDs_Not_Passed = $seqID[1];
            }
        }
    }
}

sub error_Notify {
    my $step = shift;
    my $stderr = shift;
    $logger ->logdie("Error occured at step $step: $stderr");
}

sub parse_command_line {
    my ($help, $ver);
    my $result = GetOptions (
            "step=s" => \$stepIn,
            "config=s" => \$configPath,
            "help" => \$help,
            "version" => \$ver
    );
    usage() if ($help);
    if ($ver) {
        print "PARFuMS version: $VERSION\n";
        exit;
    }
    #die "Please provide step and config file path" if (!$stepIn);
    if ($stepIn =~ /:/) {
        my @num = split(/\:/, $stepIn);
        for(my $i=$num[0]; $i <= $num[1]; $i++) {
            push(@procStep, $i);
        }
    }
    elsif ($stepIn == 0) {
        for(my $i=1; $i <= $MAXSTEP; $i++) {
            push(@procStep, $i);
        }
    }
    else {
        push(@procStep, $stepIn);
    }
    #$logger->info("Input Step Information: @procStep");
}

sub read_config_file {
    die "Couldn't locate config file at $configPath\n"  if not(-e $configPath);
    #$logger->info("Input Config File: $configPath");
    open (my $configIN, "< $configPath") or die "$!";
    chomp(my @allData = <$configIN>);
    foreach my $line (@allData) {
        if (($line !~ /^#/) && ($line !~ /^\s*$/)) {
            my @vals = split(/\s+\:\s+/, $line);
            $config{$vals[0]} = $vals[1];
        }
    }
    $unmatch_stat = 1 if ($config{'unmatch_disp'} eq 'yes');
}

sub makeOutputDir {
    if (not(-d $config{'DIR'})) {
		system("mkdir", "$config{'DIR'}");
        system("mkdir", "$config{'DIR'}/PreprocessedFiles");
        system("mkdir", "$config{'DIR'}/tmp");
        system("mkdir", "$config{'DIR'}/tmp/qsub_scripts");
        #$logger->info("New Directory: $config{'DIR'} and its subdirectories made");
	}
	else {
        #$logger->info("$config{'DIR'} exists: Overwriting files");
        if (not(-d "$config{'DIR'}/PreprocessedFiles")) {
            system("mkdir", "$config{'DIR'}/PreprocessedFiles");
        }
        if (not(-d "$config{'DIR'}/tmp")) {
            system("mkdir", "$config{'DIR'}/tmp");
            system("mkdir", "$config{'DIR'}/tmp/qsub_scripts");
        }
        else {
            if (not(-d "$config{'DIR'}/tmp/qsub_scripts")) {
                system("mkdir", "$config{'DIR'}/tmp/qsub_scripts");
            }
        }
    }
}

sub usage {
print<<EOF;

PARFUMS version $VERSION, by Manish Boolchandani (manish\@wustl.edu),
This program is a wrapper to run all steps of the PARFUMS pipeline.

Usage: perl parfum_wrapper.pl --config [CONFIG_FILE] --step [1:5|0|N1:N2];

Please find format of CONFIG_FILE in README.

Each step corresponds to a specific stage on PARFUMS pipeline.
Step information is entered to
    run a complete pipeline (enter 0),
    run a specific step (any number 1 to 5),
    or specific set of steps (2:4 to run Step 2,3,4).

Description of each step:

Step 1: Preprocessing input FastQ files.
    Reads FW- and RC- FASTQ files and splits it into smaller files based
    on barcode matching. It further converts each small FASTQ into FASTA file
    format.

Step 2: Adapter Cleaning
    Creates adapter file for each specific ID and trim sequences based on
    adapter matching.

Step 3: Clean PhiX and Rev Vector reads.
    Identify reads associated with PhiX and RevVector and removes them from
    FASTA file. (Output filename: ID.noVector.fasta)

Step 4: Velvet Assembly

Step 5: Phrap Assembly

Step 6: Annotation
EOF
exit 1;

}
