#!/usr/bin/perl -w
#
################################
#   Code: Remove vector sequences from fasta reads
#   Author: Manish Boolchandani
#   Date: 2015/03/01
###############################

use strict;
use FindBin qw($Bin);
use File::Basename;
use File::Spec::Functions qw(catdir);

use lib catdir(dirname($Bin), '/lib');
use Getopt::Long;
use Log::Log4perl;
use PARFuMs_Subs;
use File::Temp qw(tempdir tempfile);


#Global Variables
my $MAXSIZE = 200000;
my ($bc_file, $inputDir, $vec_file);
my (%spFasta, %totalReads, %cm_outfile);

#Subroutine to parse command line
&parse_command_line();

#setup log file configuration path
my $dirname = dirname($Bin);
my $log_conf = "$dirname/config_files/log4perl.conf";
Log::Log4perl::init($log_conf);
my $logger = Log::Log4perl->get_logger("Remove_Vec");
$logger -> info("Removing vector associated reads");

my $template = "vec_CMRun_XXXX";
my $tempDir = tempdir($template, DIR => "$inputDir/tmp/", CLEANUP => 0);

#Change the permissions of the temporary directory in order to make it accessible by everyone
if( -d $tempDir ) {
        $logger -> debug("$tempDir exists");
        chmod(0755, $tempDir) or $logger -> die("Unable to change the permissions of dir $tempDir");
} else {
        $logger -> debug("$tempDir does not exist");
}

#Subroutine to split fasta files
&read_FastaFile();

#Subroutine to create jobscript of crossmatch
my $cm_jobfile = &make_crossmatch_jobscript();

#Run CrossMatch and make final output file
PARFuMs_Subs::submit_jobArray($cm_jobfile, "vec_CM_Run", $tempDir, $logger);

foreach my $id (keys %cm_outfile) {
	my $clean_filename = "$inputDir/$id/$id.noVector.fasta";
	if (-e $clean_filename) {
		unlink $clean_filename or $logger ->die("Unable to remove pre-existing $clean_filename");
	}
	$logger ->debug("$id\t@{$cm_outfile{$id}}");
	foreach my $filename (@{$cm_outfile{$id}}) {
		system("cat $filename >> $clean_filename");
	}
	$logger ->debug("$clean_filename is formed");
}
#Process CrossMatch output
#&remove_Vector();

#$logger -> info("STEP 3: Completed!!");

exit 0;

sub remove_Vector {
    my (%processedIDs, %vectorReads);
    foreach my $id (keys %cm_outfile) {
	my @sub_outfiles;
        my $outfile = "$inputDir/$id/$id.noVector.fasta";
        my @cmfiles = @{$cm_outfile{$id}};
        my @fastafiles = @{$spFasta{$id}};
        for(my $i=0; $i<scalar(@cmfiles); $i++) {
            my $cmName = basename($cmfiles[$i], ".cmoutput");
            my $fastaName = basename($fastafiles[$i], ".fasta");
            if ($cmName eq $fastaName) {
		my $cmd = "perl $dirname/src/RemoveVector_Gautam.pl $fastafiles[$i] $cmfiles[$i]";
		$logger -> debug("Remove Vector CMD: $cmd");
		my $outfile = `$cmd`;
		push(@sub_outfiles, $outfile);
                #my %fastaseq = &fastafile_parser($fastafiles[$i]);
                #my @CMout = &cm_parser($cmfiles[$i]);
                #foreach my $line (@CMout) {
                #    my @vals = split(/\s+/, $line);
                #    my $flank = $1 if ($vals[8] =~ /\((\d+)\)/);
                #    my $pass = 0;
                #    $pass = 1 if ($vals[6] <= 5);
                #    $pass = 2 if ($flank <= 5);
                #    my $perOverlap = abs(($vals[7] - $vals[6])/($vals[7] + $flank));
                #    $pass = 3 if ($perOverlap > 0.85);
                #    $pass = 4 if ($processedIDs{$vals[5]});
                #    $processedIDs{$vals[5]} = 1;
                #    next unless $pass > 0;
                #    $vals[5] =~ s/\_\d$//;
                #    $vectorReads{$vals[5]} = 1;
                #}
                #foreach my $id (keys %fastaseq) {
                #    if (not exists $vectorReads{$id}) {
                #        print $OUT ">$id\_0\n$fastaseq{$id}[0]\n>$id\_1\n$fastaseq{$id}[1]\n";
                #    }
                #}
            }
            else {
                $logger->info("ERROR: $fastafiles[$i] not processed as filenames didn't match ($cmName $fastaName)");
            }
        }
	if (-e $outfile) {
		system("rm -f $outfile");
	}
	foreach my $subfile (@sub_outfiles) {
		chomp($subfile);
		my $out = `cat $subfile >> $outfile`;
	}
	my $seqC = `grep -c '^>' $outfile`;
	$logger ->debug("$outfile sequence count: $seqC");
    }
}

sub cm_parser {
    my $filename = shift;
    my @CMInfo; my @data;
    open (my $cmIN, "< $filename") or die "$!";
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
    return @CMInfo;
}

sub fastafile_parser {
    my $filename = shift;
    my %fastaseq;
    if (-e $filename) {
        open (my $in, "< $filename") or die "$!";
        while(my $line = <$in>) {
            chomp($line);
            if ($line =~ /^>/) {
                my $id1 = $line;
                $id1 =~ s/>//; $id1 =~ s/\_\d$//;
                chomp(my $fwseq = <$in>);
                chomp(my $id2 = <$in>);
                $id2 =~ s/>//; $id2 =~ s/\_\d$//;
                chomp(my $rcseq = <$in>);
                if ($id1 eq $id2) {
                    $fastaseq{$id1} = [ $fwseq,$rcseq ];
                }
                else {
                    $logger->info("ERROR: remove_Vectors.pl: FW- and RC-seq IDs didn't match ($id1 $id2)");
                }
            }
        }
    }
    return %fastaseq;
}

sub make_crossmatch_jobscript {
    my $cm_jobfile = "$tempDir/vec_CMrun_script.sh";
    $logger -> info ("Generating a script file to run cross_match: vec_CMrun_script.sh");
    open (my $OUT, "> $cm_jobfile") or $logger -> logdie("$!");
    my $codepath = "$dirname/src/RemoveVector_Gautam.pl";
    my $CM_cmd = PARFuMs_Subs::getCM_cmd();
    foreach my $id (sort {$totalReads{$a} <=> $totalReads{$b}} keys %totalReads) {
        if (exists $spFasta{$id}) {
            my @files = @{$spFasta{$id}};
            foreach my $filename (@files) {
                my $name = fileparse($filename, ".fasta");
                my $outFile = "$filename.noVector";
                my $tmp_cmd = $CM_cmd;
                $tmp_cmd =~ s/FASTA_FILE_1/$filename/g;
                $tmp_cmd =~ s/FASTA_FILE_2/$vec_file/g;
                $tmp_cmd =~ s/CODE_TO_REMOVE/$codepath/g;
                #my $cm_cmd = "cross_match $filename $inputDir/$id/Adapters_$id.fna" . $cm_param . "> $filename.cmoutput";
                print $OUT "$tmp_cmd\n";
                if (exists $cm_outfile{$id}) {
                    push(@{$cm_outfile{$id}}, $outFile);
                }
                else {
                    $cm_outfile{$id} = [ $outFile ];
                }
            }
        }
    }
    close($OUT);
    return $cm_jobfile;
}

sub read_FastaFile {
    my %processFiles;
    my $splitDir = "$tempDir/splitFastaFiles";
    if (not -d $splitDir) {
        mkdir $splitDir;
    }
    open (my $in, "< $bc_file") or die "Could not find barcode file at $bc_file: $!";
    chomp(my @BC = <$in>);
    foreach my $line (@BC) {
        my @seqID = split(/\s+/, $line);
        my $fastaFile = "$inputDir/$seqID[1]/$seqID[1].noPhiX.fasta";
        if (-e $fastaFile) {
            chomp(my $readCount = `grep '^>' $fastaFile | wc -l`);
            $totalReads{$seqID[1]} = $readCount;
            if ($readCount > $MAXSIZE) {
                $processFiles{$seqID[1]} = $fastaFile;
            }
            else {
                system("cp $fastaFile $splitDir/$seqID[1].fasta");
                $spFasta{$seqID[1]} = [ "$splitDir/$seqID[1].fasta" ];
            }
        }
    }
    my %tmpFasta = PARFuMs_Subs::split_FastaFiles(\%processFiles, $splitDir, $MAXSIZE, $logger);
    @spFasta{keys %tmpFasta} = values %tmpFasta;
}

sub get_Logfile_Name {
    return "$inputDir/parfums.log";
}

sub parse_command_line {
    my $help;
    my $result = GetOptions (
        "bcfile=s" => \$bc_file,
        "dir=s" => \$inputDir,
        "vecfile=s" => \$vec_file,
        "help" => \$help
    );
}


