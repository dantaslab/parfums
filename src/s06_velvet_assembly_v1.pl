#!/usr/bin/perl -w

####################################################################################################
# Created by:           Manish Boolchandani
# Creation Date:        2015/03/02
# Version:              1.0
# Description:          Short Read assembly using multiple rounds of Velvet
#
# Modify By:
# Modification Date:
# Remarks:
#
####################################################################################################

use strict;
use FindBin qw($Bin);
use File::Basename;
use File::Temp qw(tempdir tempfile);
use File::Spec;
use File::Spec::Functions qw(catdir);

use lib catdir(dirname($Bin), 'lib');
use Getopt::Long;
use Capture::Tiny;
use Log::Log4perl;
use Cwd;
use PARFuMs_Subs qw(split_FastaFiles submit_jobArray get_FRhit_cmd get_CDhit_cmd);

# Variable declaration
my ($sackedIDs, $bc_file, $inputDir, $paramVelveth1, $paramVelvetg1, $runType, $jobFile, $dbLink);
my (%inputSeqFile, %totalReads, %mapFile, %cdhitFile, %frhitFile, %noChimeraFile);
my (@params, @uncheckedIDs);
my $counter = 0;

#------------------------------------------------------------------------#
# START OF MAIN PROGRAM
#------------------------------------------------------------------------#

#Subroutine to parse command line arguments
&parse_command_line();

#setup log file configuration path
my $dirname = dirname($Bin);
my $log_conf = "$dirname/config_files/log4perl.conf";
Log::Log4perl::init_once($log_conf);
my $logger = Log::Log4perl -> get_logger("Velvet_Assembly");
$logger -> info("STEP 4: VELVET ASSEMBLY STARTED");

#my $MAXSIZE = 100000
my $MAXSIZE = 10000; # maximum number of reads to keep in one file
my $runDir = getcwd;
my $template = "VelvetRun1_XXXX";
my $tempDir = tempdir( $template, DIR => "$inputDir/tmp", CLEANUP => 0 );

#Change the permissions of the temporary directory in order to make it accessible by everyone
if( -d $tempDir ) {
        $logger -> debug("$tempDir exists");
        chmod(0755, $tempDir) or $logger -> die("Unable to change the permissions of dir $tempDir");
} else {
        $logger -> debug("$tempDir does not exist");
}

if ($sackedIDs) {
    @uncheckedIDs = split(/\,/, $sackedIDs);
}
#------------------------------------------------------------------------#
# Velvet Round 1                                                         #
#------------------------------------------------------------------------#

# Read barcode file and input sequence files
$logger -> info("IDs not included: @uncheckedIDs");
&read_inputSequences();

#Default parameters for velvet,
#Old_Parfum CMD: velveth assemble_040104_AX.noVector.fasta 31 -shortPaired 040104_AX.noVector.fasta; velvetg assemble_040104_AX.noVector.fasta -cov_cutoff 10 -ins_length 100 -min_contig_lgth 100
$paramVelveth1 = "31 -shortPaired";
$paramVelvetg1 = "-cov_cutoff 10 -ins_length 100 -min_contig_lgth 100";
push(@params, ($paramVelveth1, $paramVelvetg1));

#Subroutine to split _NoVector.fasta  file into $MAXSIZE sequences for velvet run
my %splitFile = &split_fastafile();

$logger -> info("Velvet Assembly Round-1 Started");
$jobFile = &generate_script("velvet",\%splitFile); # generate script file to run velvet assembly on cluster
PARFuMs_Subs::submit_jobArray($jobFile, "VelvetRun1", $tempDir, $logger); # Submit array of jobs on SGE cluster using qsub

$logger -> debug("Round-1 completed");
my $contigRef = &merge_VelvetOutput(\%splitFile, "velvet1_contigs.fasta" ,$tempDir);  #Merge velvet output files and return list of files in an hash reference
my %contigFiles = %$contigRef; #hash of contig files for each ID

$jobFile = &generate_script("cd-hit-est",\%contigFiles, "cd-hit1.fasta"); #generate script file to run cd-hit on cluster
PARFuMs_Subs::submit_jobArray($jobFile, "CD-hit1","$tempDir", $logger); #submit array of jobs

$jobFile = &generate_script("fr-hit",\%mapFile, "Map1.txt");
PARFuMs_Subs::submit_jobArray($jobFile, "FR-hit1","$tempDir",$logger); #submit array of jobs

$jobFile = &generate_script("remove-chimera",\%mapFile, "NoChimera.fasta");
PARFuMs_Subs::submit_jobArray($jobFile, "Remove-Chimera1","$tempDir",$logger); #submit array of jobs

$jobFile = &generate_script("cd-hit-est", \%noChimeraFile, "cd-hit2.fasta" );
PARFuMs_Subs::submit_jobArray($jobFile, "CD-hit2","$tempDir",$logger); #submit array of jobs

$jobFile = &generate_script("fr-hit", \%mapFile,"Map2.txt");
PARFuMs_Subs::submit_jobArray($jobFile,"FR-hit2","$tempDir",$logger);

$jobFile = &generate_script("unmapped-reads",\%frhitFile,"Missing1stPass.fasta");
PARFuMs_Subs::submit_jobArray($jobFile, "UnmappedReads1", "$tempDir",$logger);

print "Finshed First Round\n";

#-----------------------------------------------------------------#
# Velvet Round 2                                                  #
#-----------------------------------------------------------------#

#Velvet Parameters for round2
$template = "VelvetRun2_XXXX";
$tempDir = tempdir( $template, DIR => "$inputDir/tmp", CLEANUP => 0 );

#Change the permissions of the temporary directory in order to make it accessible by everyone
if( -d $tempDir ) {
        $logger -> debug("$tempDir exists");
        chmod(0755, $tempDir) or $logger -> die("Unable to change the permissions of dir $tempDir");
} else {
        $logger -> debug("$tempDir does not exist");
}

$MAXSIZE = 1500;
#$MAXSIZE = 2000;
&split_fastafile("Missing1stPass.fasta");

my $paramVelveth2 = "31 -shortPaired";
my $paramVelvetg2 = "-cov_cutoff 7 -ins_length 80 -min_contig_lgth 100";
$params[0] = $paramVelveth2;
$params[1] = $paramVelvetg2;

$logger -> info("Velvet Assembly Round-2 Started");
$jobFile = &generate_script("velvet",\%splitFile); # generate script file to run velvet assembly on cluster
PARFuMs_Subs::submit_jobArray($jobFile, "VelvetRun2", $tempDir,$logger); # Submit array of jobs on SGE cluster using qsub
$contigRef = &merge_VelvetOutput(\%splitFile, "velvet2_contigs.fasta" , $tempDir);  #Merge velvet output files and return list of files in an hash reference
%contigFiles = %$contigRef; #hash of contig files for each ID

$jobFile = &generate_script("combine-contigs", \%inputSeqFile, "ForCD-hit.fasta");
PARFuMs_Subs::submit_jobArray($jobFile, "CombineContigs1", $tempDir, $logger);

$jobFile = &generate_script("cd-hit-est", \%contigFiles, "cd-hit3.fasta");
PARFuMs_Subs::submit_jobArray($jobFile,"CD-hit3", "$tempDir",$logger);

$jobFile = &generate_script("fr-hit",\%mapFile, "Map3.txt");
PARFuMs_Subs::submit_jobArray($jobFile, "FR-hit3","$tempDir",$logger); #submit array of jobs

$jobFile = &generate_script("unmapped-reads",\%frhitFile,"Missing2ndPass.fasta");
PARFuMs_Subs::submit_jobArray($jobFile, "UnmappedReads2", "$tempDir",$logger);

print "\nFinished 2nd round of Velvet\n";

#-----------------------------------------------------------------#
# Velvet Round 3                                                  #
#-----------------------------------------------------------------#

#Velvet Parameters for round 3
$template = "VelvetRun3_XXXX";
$tempDir = tempdir( $template, DIR => "$inputDir/tmp", CLEANUP => 0 );

#Change the permissions of the temporary directory in order to make it accessible by everyone
if( -d $tempDir ) {
        $logger -> debug("$tempDir exists");
        chmod(0755, $tempDir) or $logger -> die("Unable to change the permissions of dir $tempDir");
} else {
        $logger -> debug("$tempDir does not exist");
}

my $paramVelveth3 = "31 -shortPaired";
my $paramVelvetg3 = "-cov_cutoff 10 -ins_length 80 -min_contig_lgth 100";
$params[0] = $paramVelveth3;
$params[1] = $paramVelvetg3;

$logger -> info("Velvet Assembly Round-3 Started");
$jobFile = &generate_script("velvet",\%splitFile); # generate script file to run velvet assembly on cluster
PARFuMs_Subs::submit_jobArray($jobFile, "VelvetRun3", $tempDir, $logger); # Submit array of jobs on SGE cluster using qsub
$contigRef = &merge_VelvetOutput(\%splitFile, "velvet3_contigs.fasta" , $tempDir);  #Merge velvet output files and return list of files in an hash reference
%contigFiles = %$contigRef; #hash of contig files for each ID

$jobFile = &generate_script("combine-contigs", \%inputSeqFile, "ForCD-hit2.fasta");
PARFuMs_Subs::submit_jobArray($jobFile, "CombineContigs2", $tempDir, $logger);

$jobFile = &generate_script("blast-adapters", \%contigFiles, "ForCD-hit3.fasta");
PARFuMs_Subs::submit_jobArray($jobFile, "blastAdapter", $tempDir, $logger);

$jobFile = &generate_script("cd-hit-est", \%contigFiles, "cd-hit4.fasta");
PARFuMs_Subs::submit_jobArray($jobFile,"CD-hit4", "$tempDir", $logger);

$jobFile = &generate_script("fr-hit",\%mapFile, "Map4.txt");
PARFuMs_Subs::submit_jobArray($jobFile, "FR-hit4","$tempDir", $logger); #submit array of jobs

$jobFile = &generate_script("remove-chimera",\%mapFile, "ForPhrap1.fasta");
PARFuMs_Subs::submit_jobArray($jobFile, "Remove-Chimera2","$tempDir", $logger); #submit array of jobs

#print "\n Finished 3rd round of Velvet assembly\n";
#print "contigFiles => \n";
#print %contigFiles;
#print "\n No Chimera CD-hit files cdhitFile =>\n";
#print %cdhitFile;
#print "\n Map file mapFile=> \n";
#print %mapFile;
#print "\n splitFile => \n";
#print %splitFile;
#print "\n frhitFile =>\n";
#print %frhitFile;
#print "\n inputSeqFile =>\n";
#print %inputSeqFile;


$logger -> info("End of Velvet Assembly program");
exit(0);

#END OF MAIN PROGRAM

####################################################################################################

# subroutine to read input sequences from barcode file

sub read_inputSequences{
    my $suffix = "noVector.fasta";
    open (my $IN,"<$bc_file") or $logger -> logdie("$!");
    chomp(my @BC = <$IN>);
    my $rec = join(":", @uncheckedIDs);
    $logger ->info("IDs not included: $rec");
    foreach my $line (@BC){
        my @seqID = split(/\s+/,$line);
        my $filename = "$seqID[1].$suffix";
        if (not grep {$_ eq $seqID[1]} @uncheckedIDs) {
            my $filepath = "$inputDir/$seqID[1]/$filename";
            if( -e $filepath){
                $inputSeqFile{$seqID[1]} = $filepath;
            }
            else{
                $logger -> logdie("ERROR: File doesnt exists: $filepath");
            }
        }
    }
    $logger -> debug("barcode file stored");
}

#Subroutine to split *_NoVector.fna fasta file into $MAXSIZE sequences for velvet run
sub split_fastafile{
    my %processFiles; my %splitFile;
    my $outDir = "$tempDir/splitFastaFiles";
    if (not -d $outDir) {
        mkdir $outDir;
    }
    foreach my $sampleID (keys %inputSeqFile){
        my $fastaFile = $inputSeqFile{$sampleID};
        if(-e $fastaFile){
            $mapFile{$sampleID} = $fastaFile;
            chomp (my $readCount = `grep -c '^>' $fastaFile`);
            $totalReads{$sampleID} = $readCount;
            if ($readCount > $MAXSIZE){
                $processFiles{$sampleID} = $fastaFile;
            }
            else{
                system("cp $fastaFile $outDir/$sampleID.fasta");
                $splitFile{$sampleID} = [ "$outDir/$sampleID.fasta" ];
            }
        }
        else{
            $logger -> logdie("ERROR: File doesn't exists: $fastaFile");
        }
    }
    $logger -> debug("Split Fasta Files Process Started: MAXSIZE: $MAXSIZE");
    my %tmpFasta = PARFuMs_Subs::split_FastaFiles(\%processFiles, $outDir, $MAXSIZE, $logger);
#   Merging tmpFasta with splitFile array
    @splitFile{keys %tmpFasta} = values %tmpFasta;
    $logger -> debug("Split Fasta Files complete");
    return %splitFile;
}

sub generate_script {
    my $runType = $_[0];
    my %splitFile = %{$_[1]};
    my $suffix="";
    my $jobfile = "$inputDir/tmp/qsub_scripts/${runType}_script_v$counter.sh";
    my $jobfilename = basename($jobfile);
    my @cmds;
    if(scalar(@_) == 3){
        $suffix =$_[2];
    }

    open (my $OUT, "> $jobfile") or $logger -> logdie($!);
    $logger -> info("Generating qsub script file to run $runType:$jobfilename");

    if($runType =~ "velvet"){
        print "Writing Velvet script file\n";
        @cmds = get_VelvetCmd(\@params, \%splitFile);
    }
    elsif($runType =~ "cd-hit-est") {
        $logger -> info("Writing CD-HIT script file");
	foreach my $k (keys %splitFile) { $logger -> debug("CD-HIT file: $k\t$splitFile{$k}"); }
        @cmds = get_CDhitCmd(\%splitFile, $suffix);
    }
    elsif($runType =~ "fr-hit" ){
	$logger -> info("Writing FR-HIT script file");
        print "Writing FR-Hit script file\n";
        @cmds = get_FRhitCmd(\%inputSeqFile, \%cdhitFile, $suffix);
    }
    elsif($runType =~ "remove-chimera"){
        $logger -> info("Writing Remove-Chimera script file");
        my $cmdRef = get_RmChimeraCmd(\%cdhitFile,\%frhitFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType =~ "unmapped-reads"){
        $logger -> info("Writing script file to get unmapped reads");
        my $cmdRef = get_UnmappedReadsCmd(\%mapFile,\%frhitFile,$suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType =~ "combine-contigs"){
        $logger -> info("Writing script to combine contig files");
        my $cmdRef = get_combineContigsCmd(\%splitFile, \%cdhitFile, \%contigFiles, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType =~ "blast-adapters"){
        $logger -> info("Writing script file to blast contigs againt adapters database");
        my $cmdRef = get_BlastAdaptCmd(\%contigFiles, $suffix);
        @cmds = @$cmdRef;
    }
    else {
        $logger -> logdie("Invalid runType: $runType");
    }
    #$logger->info("$runType CMD: @cmds");
    if ( @cmds > 0 ){
        foreach my $cmd (@cmds){
            print $OUT "$cmd\n";
        }
    } else {
        $logger -> logdie("Nothing in the command list to submit on cluster!");
    }
    $counter++;
    close($OUT);
    return $jobfile;
}

sub get_VelvetCmd {
    my @parameters = @{$_[0]};
    my $paramVelveth = $parameters[0];
    my $paramVelvetg = $parameters[1];
    my %splitFile = %{$_[1]};
    my @cmds;
    my $velvet_cmd = "velveth assemble_INCLUDE_DIR $paramVelveth INCLUDE_FILE ; velvetg assemble_INCLUDE_DIR $paramVelvetg";
    foreach my $id ( keys %splitFile){
        my @files = @{$splitFile{$id}};
        foreach my $file (@files){
            my $filename = basename($file, ".fasta");
            my $tmp_cmd = $velvet_cmd;
            $tmp_cmd =~ s/INCLUDE_DIR/$filename/g;
            $tmp_cmd =~ s/INCLUDE_FILE/$file/g;
            push(@cmds, $tmp_cmd);
        }
    }
    return @cmds;
}

sub get_FRhitCmd {
    my %mapFile = %{$_[0]};
    my %cdhitFile = %{$_[1]};
    my $suffix = $_[2];
    my @cmds;
    foreach my $sampleId ( keys %mapFile) {
        my $file = $mapFile{$sampleId};
        my $dirname = dirname($file);
        my $INFILE = $cdhitFile{$sampleId};
        my $ALIGNFILE = $file;
        my $OUTFILE = "$dirname/${sampleId}.$suffix";
        my $frhit_cmd = "fr-hit -d $INFILE -a $ALIGNFILE -m 30 -o $OUTFILE";
        $frhitFile{$sampleId} = $OUTFILE;
        push(@cmds, $frhit_cmd);
    }
    return @cmds;
}

sub get_CDhitCmd {
    my %splitFile = %{$_[0]};
    my $suffix = $_[1];
    my @cmds;
    foreach my $sampleId ( keys %splitFile){
        my $file = $splitFile{$sampleId};
        my $filename = basename($file,".fasta");
        my $dirname = dirname($file);
        my $INFILE = $file;
        my $OUTFILE = "$dirname/${sampleId}.$suffix";
	$logger -> debug("Writing CD-HIT script file");
	#Previous cmd line didn't include option of -c 0.9 which enforce >90% seq identity to remove redundant contigs (It might be by default option)
        #my $cdhit_cmd = "cd-hit-est -i $INFILE -o $OUTFILE -g 1 -r 1";
        my $cdhit_cmd = "cd-hit-est -i $INFILE -o $OUTFILE -g 1 -r 1 -c 0.9";
        $cdhitFile{$sampleId} = $OUTFILE;
        push(@cmds, $cdhit_cmd);
    }
    $logger -> debug("@cmds");
    return @cmds;
}

sub get_RmChimeraCmd {
    my %cdhitFile = %{$_[0]};
    my %frhitFile = %{$_[1]};
    my $suffix = $_[2];
    my @cmds;
    foreach my $sampleId ( keys %totalReads){
        if(exists $cdhitFile{$sampleId}){
            my $file = $cdhitFile{$sampleId};
            my $dirname = dirname($file);
            my $MAPFILE = $frhitFile{$sampleId};
            my $CDHITFILE = $file;
            my $OUTFILE = "$dirname/$sampleId.$suffix";
            my $rmchimera_cmd = "perl $Bin/FR-Hit_cleanChimera.pl $MAPFILE $CDHITFILE > $OUTFILE";
               $noChimeraFile{$sampleId} = $OUTFILE;
            push(@cmds, $rmchimera_cmd);
        }
    }
    return \@cmds;
}

sub merge_VelvetOutput {
    my %velvetFiles = %{$_[0]};
    my $suffix = $_[1];
    my $velvet_OutDir = $_[2];
    my %contigFiles;
    foreach my $seqID (keys %velvetFiles){
        my @allfiles = @{$velvetFiles{$seqID}};
        my $combinedFile = "$seqID.$suffix";
        my $outFile = "$inputDir/$seqID/$combinedFile";
	if (-e $outFile) {
		$logger -> debug("$outFile exists!! will be deleted");
		system("rm -f $outFile");
	}
        foreach my $file (@allfiles){
            my $filename = basename($file, ".fasta");
            my $cmd = "cat $velvet_OutDir/assemble_$filename/contigs.fa >> $outFile";
            system("$cmd");
            $contigFiles{$seqID} = "$outFile";
        }
    }
    return \%contigFiles;
}

sub get_UnmappedReadsCmd{
    my %mapFile = %{$_[0]};
    my %frhitFile = %{$_[1]};
    my $suffix = $_[2];
    my @cmds;
    foreach my $sampleId ( keys %totalReads){
        if(exists $mapFile{$sampleId}){
            my $file = $mapFile{$sampleId};
            my $dirname = dirname($file);
            my $INFILE = $frhitFile{$sampleId};
            my $MAPFILE = $file;
            my $OUTFILE = "$dirname/${sampleId}.$suffix";
            my $unmapped_cmd = "perl $Bin/FR-Hit_get_unmapped.pl $INFILE $MAPFILE $sampleId > $OUTFILE";
            push(@cmds, $unmapped_cmd);
            $splitFile{$sampleId} = [$OUTFILE];
        }
    }
    return \@cmds;
}

sub get_combineContigsCmd{
    my %seqFile = %{$_[0]};
    my %cdhitFile = %{$_[1]};
    my %contigFile = %{$_[2]};
    my $suffix = $_[3];
    my @cmds;
    foreach my $sampleId (keys %seqFile){
        if(exists $seqFile{$sampleId}){
            my $INFILE1 = $cdhitFile{$sampleId};
            my $INFILE2 = $contigFile{$sampleId};
            my $dirname = dirname($seqFile{$sampleId});
            my $OUTFILE = "$dirname/$sampleId.$suffix";
            my $tmp_cmd = "cat $INFILE1 $INFILE2 > $OUTFILE";
            push(@cmds, $tmp_cmd);
            $contigFiles{$sampleId} = $OUTFILE;
        }
    }
    return \@cmds;
}

sub get_BlastAdaptCmd{
    my %contigFile = %{$_[0]};
    my $suffix = $_[1];
    my @cmds;

    foreach my $sampleId (keys %inputSeqFile){
        if(exists $contigFile{$sampleId}){
            my $INFILE = $contigFile{$sampleId};
            my $ADAPTFILE = "$inputDir/$sampleId/Adapters_$sampleId.fna";
            my $dirname = dirname($contigFile{$sampleId});
            my $OUTFILE = "$dirname/$sampleId.$suffix";
            #my $tmp_cmd = "makeblastdb -in $ADAPTFILE -dbtype nucl; /srv/cgs/local/ncbi-blast/latest/bin/blastn -db $ADAPTFILE";
            my $tmp_cmd = "makeblastdb -in $ADAPTFILE -dbtype nucl; blastn -db $ADAPTFILE";
            $tmp_cmd = $tmp_cmd . " -query $INFILE -outfmt 6 -dust no -evalue 1e-4 | $Bin/CleanAdapterContigs.pl $INFILE > $OUTFILE";
            push(@cmds, $tmp_cmd);
            $contigFiles{$sampleId} = $OUTFILE;
        }
    }
    return \@cmds;
}

sub remove_temporary_files{
    File::Temp::cleanup();
}

sub get_Logfile_Name {
    return "$inputDir/parfums.log";
}

sub parse_command_line {
    my $help;
    my $result = GetOptions (
        "bcfile=s"          => \$bc_file,
        "dir=s"             => \$inputDir,
	"dblink=s"	    => \$dbLink,
        "paramVelveth1:s"   => \$paramVelveth1,
        "paramVelvetg1:s"   => \$paramVelvetg1,
        "help"              => \$help,
        "unchecked=s"         => \$sackedIDs
    ) or $logger -> logdie("Error in command line arguments: $!");

    &usage()if ($help);
}

sub usage {

print<<EOF;
    PaRFuMs velet_assembly, by Manish Boolchandani (manish\@wustl.edu),

    This program reads _NoVector fasta file and splits it into several
    smaller files, based on seqIds in the barcode file and then run velvet assembly on
    each file.

    usage: $0 --bcfile FILE --dir OUTPUT_DIR [--help]

    Arguments:

    --bcfile FILE       - Barcodes file name. (Barcodes and Identifiers separate by tab)
    --dir OUTPUT_DIR    - Directory in which all output files will be saved.
    --help              - This helpful help screen.

EOF
exit 1;
}


