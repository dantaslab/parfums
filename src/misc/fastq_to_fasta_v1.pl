#!/usr/bin/perl -w

#
# Created by: Manish Boolchandani
# Creation Date: 2015-02-15
# Version: 1.0
# Description: Converts fastq to fasta format
#
#
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

#Global variables to store command line input
my ($inputDir, $bcfile, $help);

#START OF THE PROGRAM
&parse_command_line();

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
        while (my $line = <$bcIN>) {
            chomp($line);
            my @vals = split(/\s+/, $line);
            my $filename = $vals[1];
            my $file_path = "$inputDir/$filename";
            if (not -d $file_path) {
                mkdir "$file_path" or $logger -> logdie("Error creating directory: $!");
            }
            my $fasta_file = "$file_path/$filename.fasta";
            $logger -> debug("Reading file: $filename.txt");

            print "$dirname/lib/seqtk-master/seqtk seq -a $fastqDir/$filename.txt > $fasta_file\n";

            my ($std_out,$std_err,$exit_status)= capture {
                system("$dirname/lib/seqtk-master/seqtk seq -a $fastqDir/$filename.txt > $fasta_file");
            };
            if($exit_status){
                $logger -> logdie("Error in executing seq-tk program:$std_err");
            }
            else {
                $logger -> debug("$filename.fasta created");
                gzip $fasta_file => "$fasta_file.gz" or $logger -> logdie("gzip failed:$GzipError");
            }
        }
    }
    else {
        print "Error: Directory doesnt exists\n";
        exit;
    }

    #$logger -> info("Format Conversion done successfully");

    #END OF THE PROGRAM

sub get_Logfile_Name {
    return "$inputDir/parfums.log";
}

sub parse_command_line {
    my $result = undef;
    $result = GetOptions ( "dir=s" => \$inputDir,
		  "bcfile=s"  => \$bcfile,
		  "help" => \$help
		  );
    &usage() if ($help);
}
