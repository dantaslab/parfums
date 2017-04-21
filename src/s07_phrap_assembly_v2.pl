#!/usr/bin/perl -w

####################################################################################################
# Created by:           Manish Boolchandani
# Creation Date:        2015/03/09
# Version:              2.0
# Description:          Phrap assembly
#
# Modify By:		Manish Boolchandani
# Modification Date:	2016/03/16
# Remarks:		Changed logger commands, so that it prints only high level information and
#			detailed logs will be printed, when logger is SET to DEBUG mode.
# Modification Date: 	2016/05/12
# Remarks:		1) Added input parameter "resannotate" for resfam annotation mapping file
#                       2) Changed annotation function, such that program first concatenates contigs 
#			from all samples and run Molly's annotation script (Resfams) 
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

my(%inputSeqFile, %phrapFile, %contigFile, %splitFile, %totalReads, %annotationFile , %blastCOGFile, %cdhitFile, %frhitFile, %mapFile, %graphFile, %cogFile);
my(@params, @uncheckedIDs);
my($bc_file, $inputDir, $paramPhrap1, $paramPhrap2, $jobFile, $sackedIDs, $resAnnotate,$dbLink);

#------------------------------------------------------------------------#
# START OF MAIN PROGRAM                                                  #
#------------------------------------------------------------------------#
&parse_command_line();

#setup log file configuration path
my $maindir = dirname($Bin);
my $log_conf = "$maindir/config_files/log4perl.conf";
Log::Log4perl::init_once($log_conf);
my $logger = Log::Log4perl -> get_logger("Phrap_Assembly");

$logger -> info("Step 5: PHRAP ASSEMBLY STARTED");

my $MAXSIZE = 80;
my $counter = 0;
my $template = "PhrapRun1_XXXX";
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
$logger ->info("UncheckedIDs: @uncheckedIDs\tSackedIDs: $sackedIDs");
$logger -> info("Reading noVector file");
%inputSeqFile = &read_inputSequences("noVector.fasta");

$logger -> info("Reading ForPhrap file");
%phrapFile = &read_inputSequences("ForPhrap1.fasta");

#print %inputSeqFile;
#print "\n------------------\n";
#print %phrapFile;
#print "\n------------------\n";

#-----------------------------------------------------------------------#
# Phrap Round 1								#
#-----------------------------------------------------------------------#
$paramPhrap1 = "-minmatch 25 -maxmatch 40 -bandwidth 1 -minscore 30";

my $paramPhrap = $paramPhrap1;
$logger -> info("Writing script file for Phrap");
$jobFile = &generate_script("phrap", \%phrapFile, "phrapContigs.fasta"); # generate script file to run phrap assembly on cluster
$logger -> info("Phrap Assembly Round - 1 Started");

PARFuMs_Subs::submit_jobArray($jobFile, "PhrapRun1", $tempDir, $logger); # Submit array of jobs on SGE cluster using qsub
$logger -> info("Phrap Round-1 Completed");

#### Map reads to the contigs and clean for potential new chimeras
$jobFile = &generate_script("cd-hit-est", \%contigFile, "phrap.cdhit1");
PARFuMs_Subs::submit_jobArray($jobFile, "CD-hit1", $tempDir, $logger);
$logger -> info("CD-hit-est Run complete");

### Map reads again and link between contigs, phrap again.
$jobFile = &generate_script("fr-hit", \%cdhitFile, "phrap_Map1.frhit");
PARFuMs_Subs::submit_jobArray($jobFile, "FR-hit1", $tempDir, $logger);
$logger -> info("FR-hit Run complete");

$jobFile = &generate_script("link-contigs", \%frhitFile, "ForPhrap2.fasta");
PARFuMs_Subs::submit_jobArray($jobFile, "Link-contigs", $tempDir, $logger);
$logger -> info("Link-contigs Run complete");

#-----------------------------------------------------------------------#
# Phrap Round 2                                                         #
#-----------------------------------------------------------------------#
$logger -> info("Phrap Assembly Round-2 Started");

$template = "PhrapRun2_XXXX";
$tempDir = tempdir( $template, DIR => "$inputDir/tmp", CLEANUP => 0 );

#Change the permissions of the temporary directory in order to make it accessible by everyone
if( -d $tempDir ) {
        $logger -> debug("$tempDir exists");
        chmod(0755, $tempDir) or $logger -> die("Unable to change the permissions of dir $tempDir");
} else {
        $logger -> debug("$tempDir does not exist");
}

$jobFile = &generate_script("phrap",\%phrapFile, "phrapContigs2.fasta"); # generate script file to run phrap assembly on cluster
PARFuMs_Subs::submit_jobArray($jobFile, "PhrapRun2", $tempDir, $logger); # Submit array of jobs on SGE cluster using qsub
$logger -> info("Phrap Round-2 Completed");

$jobFile = &generate_script("cd-hit-est", \%contigFile, "phrap.cdhit2");
PARFuMs_Subs::submit_jobArray($jobFile, "CD-hit2", $tempDir, $logger);
$logger -> info("CD-hit2 Run complete");

&split_fastafile("phrap.cdhit2", $MAXSIZE, 1);

$jobFile = &generate_script("blastCOG", \%splitFile, "blastCOG");
$logger -> info("blastCOG Run started");
PARFuMs_Subs::submit_jobArray($jobFile, "blastCOG", $tempDir, $logger);
$logger -> info("blastCOG Run complete");

$logger -> info("Merging Blast output");
my $blastCogRef = &merge_BlastOutput(\%splitFile, "blastCOG");
$logger -> info("Merging Finished");
%blastCOGFile = %$blastCogRef;

$jobFile = &generate_script("parse-blast", \%blastCOGFile, "annotation.txt");
PARFuMs_Subs::submit_jobArray($jobFile, "ParseBlastCOG", $tempDir, $logger);
$logger -> info("parse-blast run completed");

$MAXSIZE = 10000;
&split_fastafile("noVector.fasta", $MAXSIZE, 0);
$logger -> info("Split Fasta Files completed. MaxSize = $MAXSIZE");

$jobFile = &generate_script("make-Blastdb",\%splitFile);
PARFuMs_Subs::submit_jobArray($jobFile,"makeBlastdb",$tempDir,$logger);
$logger -> info("Created Blast database");

$jobFile = &generate_script("blastForStitch", \%annotationFile, "BlastForStitch");
PARFuMs_Subs::submit_jobArray($jobFile, "BlastForStitch", $tempDir, $logger);
$logger -> info("blastforStitch Run completed");

my $blastRef = &merge_BlastOutput(\%splitFile, "BlastForStitch");
$logger -> info("Merging Blast output");
my %blastStitch = %$blastRef;

$jobFile = &generate_script("make-graph", \%blastStitch, "graph.dot");
PARFuMs_Subs::submit_jobArray($jobFile, "MakeGraph", $tempDir, $logger);
$logger -> info("makeGraph run completed");

$jobFile = &generate_script("stitch-contigs", \%graphFile, "Stitched.fasta");
PARFuMs_Subs::submit_jobArray($jobFile, "StitchContigs", $tempDir, $logger);
$logger -> info("Stitch-contigs run completed");

$jobFile = &generate_script("cd-hit-est", \%contigFile, "cdhit-last.fasta");
PARFuMs_Subs::submit_jobArray($jobFile, "CD-hit2", $tempDir, $logger);
$logger -> info("CD-hit2 run completed");

$jobFile = &generate_script("fr-hit", \%cdhitFile, "Map5.frhit");
PARFuMs_Subs::submit_jobArray($jobFile, "FR-hit2", $tempDir, $logger);
$logger -> info("FR-hit2 run completed");

##### Rename Final contigs adding meadian coverage and length
$jobFile = &generate_script("rename-contigs", \%frhitFile, "last-contigs.fasta");
PARFuMs_Subs::submit_jobArray($jobFile, "RenameContigs", $tempDir, $logger);
$logger -> info("RenameContigs run completed");

#-----------------------------------------------------------------#
# Annotation of Contigs                                           #
#-----------------------------------------------------------------#
#$logger -> info("Step 6: Started PIPELINE FOR ANNOTATION OF CONTIGS");
#
#$template = "Annotate_XXXX";
#$tempDir = tempdir( $template, DIR => "$inputDir/tmp", CLEANUP => 0 );
#
#$logger -> info("Running Blast against COG database");
#
###### Blast against COG to do functional annotation
#$MAXSIZE = 5;

# subroutine to read input sequences from barcode file
sub read_inputSequences{
    my $suffix = $_[0];
    my %seqFile;
    open (my $IN,"<$bc_file") or $logger -> logdie("$!");
    chomp(my @BC = <$IN>);
    foreach my $line (@BC){
        my @seqID = split(/\s+/,$line);
        if (not grep {$_ eq $seqID[1]} @uncheckedIDs) {
            my $filename = "$seqID[1].$suffix";
            my $filepath = "$inputDir/$seqID[1]/$filename";
            if( -e $filepath){
                $seqFile{$seqID[1]} = $filepath;
            }
            else {
                $logger -> logdie("File doesnt exists: $filepath");
            }
        }
    }
    foreach my $k (keys %seqFile) {
        print "$k\t$seqFile{$k}\n";
	$logger -> debug("Writing seqFile hash");
	$logger -> debug("$k\t$seqFile{$k}");
    }
    return %seqFile;
}

sub split_fastafile {

    my $suffix = shift;
    my $MAXSIZE = shift;
    my $multiFlag = shift;
    my %processFiles;
    my $splitdir = "$tempDir/splitFastaFiles";
    if (not -d $splitdir) {
        mkdir $splitdir;
    }
    foreach my $sampleId (keys %inputSeqFile){
        my $filename = "$sampleId\.$suffix";
        my $fastaFile = "$inputDir/$sampleId/$filename";
        if(-e $fastaFile){
            $mapFile{$sampleId} = $fastaFile;
            chomp (my $readCount = `grep -c '^>' $fastaFile`);
            $totalReads{$sampleId} = $readCount;
            if ($readCount > $MAXSIZE){
                $processFiles{$sampleId} = $fastaFile;
            }
            else{
                system("cp $fastaFile $splitdir");
                $splitFile{$sampleId} = [ "$splitdir/$filename" ];
            }
        }
        else{
            $logger -> logdie("File doesn't exists: $fastaFile");
        }
    }
    $logger -> debug("SplitFasta run for $suffix: MAXSIZE: $MAXSIZE");
    my %tmpFasta = PARFuMs_Subs::split_FastaFiles(\%processFiles, $splitdir, $MAXSIZE, $logger, $multiFlag);
    @splitFile{keys %tmpFasta} = values %tmpFasta;
}

sub generate_script {
    my $runType = $_[0];
    my %seqFile = %{$_[1]};
    #my $suffix = $_[2];
    my $suffix="";
    my $jobfile = "${inputDir}/tmp/qsub_scripts/${runType}_script_p${counter}.sh";
    my $jobfilename = basename($jobfile);
    my @cmds;
    my $param = $paramPhrap;
    if(scalar(@_) == 3){
        $suffix =$_[2];
    }

    open (my $OUT, "> $jobfile") or $logger -> logdie($!);
    $logger -> info("Generating a script file to run $runType: $jobfilename");

    if($runType =~ "phrap"){
        $logger -> info("Writing Phrap script file");
        my $cmdRef = get_PhrapCmd($param, \%seqFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType =~ "cd-hit-est") {
        $logger -> info("Writing CD-HIT script file");
        my $cmdRef = get_CDhitCmd(\%seqFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType =~ "fr-hit" ){
        $logger -> info("Writing FR-Hit script file");
        my $cmdRef = get_FRhitCmd(\%inputSeqFile, \%cdhitFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType =~ "link-contigs"){
        $logger -> info("Writing script to link phrap contigs file");
        my $cmdRef = get_LinkContigsCmd(\%cdhitFile, \%frhitFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType =~ "blastCOG"){
        $logger -> info("Writing script file to blast against COG");
        my $cmdRef = get_BlastCOGCmd(\%seqFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType eq "parse-blast"){
        $logger -> info("Writing script to parse blast COG output");
        my $cmdRef = get_parseBlastCOGCmd(\%seqFile, \%cdhitFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType =~ "make-Blastdb"){
    	$logger -> info("Writing script file to make blast db");
	    my $cmdRef = get_makeBlastdb(\%seqFile, \%cdhitFile);
        @cmds = @$cmdRef;
    }
    elsif($runType =~ "blastForStitch"){
        $logger -> info("Writing script file to blast cdhit file");
        my $cmdRef = get_BlastForStitchCmd(\%seqFile, \%splitFile, \%cdhitFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType =~ "make-graph"){
        $logger -> info("Writing script to make graph");
        my $cmdRef = get_MakeGraphCmd(\%seqFile, \%cdhitFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType =~ "stitch-contigs"){
        $logger -> info("Writing script to stitch contigs");
        my $cmdRef = get_StitchContigsCmd(\%seqFile, \%cdhitFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType =~ "rename-contigs"){
        $logger -> info("Writing script to rename contigs");
        my $cmdRef = get_RenameContigsCmd(\%seqFile, \%cdhitFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType =~ "annotate"){
        $logger -> info("Writing script file to annotate contigs");
        my $cmdRef = get_AnnotateContigsCmd(\%seqFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType eq 'parse-blast-final') {
    	$logger -> info("Writing script to parse blast COG Final output");
	    my $cmdRef = get_parseBlastCOGCmd(\%seqFile, \%contigFile, $suffix);
	    @cmds = @$cmdRef;
    }
    else {
        $logger -> logdie("Invalid runType: $runType");
    }

    if ( @cmds > 0 ){
        foreach my $cmd (@cmds){
            print $OUT "$cmd\n";
        }
    }
    else {
        $logger -> logdie("Nothing in the command list to submit on cluster!");
    }
    $counter++;
    close($OUT);
    return $jobfile;
}

sub get_PhrapCmd {
    my $parameter = $_[0];
    my %seqFile = %{$_[1]};
    my $suffix = $_[2];
    my @cmds;
    
    $logger -> debug("inside get phrap cmd function");
    foreach my $sampleId ( keys %seqFile){
        my $file = $seqFile{$sampleId};
        if( -e $file ){
            my $filename = basename($file);
            my $dirname = dirname($file);
            my $INFILE = $file;
            my $OUTFILE = "$dirname/$sampleId\.$suffix";
            my $tmp_cmd = "phrap $parameter $INFILE && cat $file\.contigs $file\.singlets $file\.problems > $OUTFILE";
            $logger -> debug("tmp_cmd=$tmp_cmd");
	    $contigFile{$sampleId} = $OUTFILE;
            push(@cmds, $tmp_cmd);
        }
    }
    return \@cmds;
}

sub get_FRhitCmd {
    my %inputSeqFile = %{$_[0]};
    my %cdhitFile = %{$_[1]};
    my $suffix = $_[2];
    my @cmds;
    foreach my $sampleId ( keys %cdhitFile){
        my $file = $cdhitFile{$sampleId};
        print "$sampleId\t$file\n";
        if(-e $file){
            my $dirname = dirname($file);
            my $INFILE = $file;
            my $ALIGNFILE = $inputSeqFile{$sampleId};
            my $OUTFILE = "$dirname/${sampleId}\.$suffix";
            my $tmp_cmd = "fr-hit -d $INFILE -a $ALIGNFILE -m 30 -o $OUTFILE";
               $frhitFile{$sampleId} = $OUTFILE;
            push(@cmds, $tmp_cmd);
        }
    }
    return \@cmds;
}

sub get_CDhitCmd {
    my %seqFile = %{$_[0]};
    my $suffix = $_[1];
    my @cmds;
    foreach my $sampleId ( keys %seqFile){
        my $file = $seqFile{$sampleId};
        if( -e $file){
            my $filename = basename($file);
            my $dirname = dirname($file);
            my $INFILE = $file;
            my $OUTFILE = "$dirname/${sampleId}\.$suffix";
            my $tmp_cmd = "cd-hit-est -i $INFILE -o $OUTFILE -g 1 -r 1 -c 0.9";
            $cdhitFile{$sampleId} = $OUTFILE;
            push(@cmds, $tmp_cmd);
        }
    }
    return \@cmds;
}

sub get_LinkContigsCmd {
    my %cdhitFile = %{$_[0]};
    my %frhitFile = %{$_[1]};
    my $suffix = $_[2];
    my @cmds;
    foreach my $sampleId ( keys %cdhitFile){
        my $file = $cdhitFile{$sampleId};
        if(-e $file){
            my $dirname = dirname($file);
            my $MAPFILE = $frhitFile{$sampleId};
            my $CDHITFILE = $file;
            my $OUTFILE = "$dirname/$sampleId\.$suffix";
            my $tmp_cmd = "perl $Bin/FR-Hit_link_light_v2.pl $MAPFILE $CDHITFILE $sampleId $OUTFILE";
            $phrapFile{$sampleId} = $OUTFILE;
            push(@cmds, $tmp_cmd);
        }
    }
    return \@cmds;
}

sub get_BlastCOGCmd {
    my %splitFile = %{$_[0]};
    my $suffix = $_[1];
    my @cmds;
    foreach my $sampleId ( keys %splitFile){
        if(exists $splitFile{$sampleId}){
            my @files = @{$splitFile{$sampleId}};
            foreach my $file (@files){
                my $filename = basename($file, ".fasta");
                my $dirname = dirname($file);
                my $INFILE = $file;
                my $OUTFILE = "$dirname/$filename\.$suffix";
                my $tmp_cmd = "blastx -outfmt 6 -evalue 1E-4 -seg no -db $dbLink/COG_with_len.fa -query $INFILE > $OUTFILE";
		#$splitFile{$sampleId} = $OUTFILE;
                push(@cmds, $tmp_cmd);
            }
        }
    }
    return \@cmds;
}

sub merge_BlastOutput {
    my %splitFile = %{$_[0]};
    my $suffix = $_[1];
    my %blastCOGFile;
    foreach my $sampleId (keys %splitFile){
	    $logger -> debug("sampleId: $sampleId");
            my @allfiles = @{$splitFile{$sampleId}};
            my $combinedFile = "$inputDir/$sampleId/$sampleId\.$suffix";
            if(-e $combinedFile){ 
		$logger -> debug("$combinedFile exists!! will be deleted"); 
		unlink $combinedFile or $logger -> logdie("Unable to delete file $combinedFile : $!");
	    }
            foreach my $file (@allfiles) {
		my $filename = basename($file, ".fasta");
		my $dirname = dirname($file);
		my $tmpFile = "$dirname/$filename\.$suffix";
		#$logger->info("fileName inside merge_blastoutput function: $tmpFile");
		#$logger->info("combined_file_name:$combinedFile");
                my $cmd = "cat $tmpFile >> $combinedFile";
		my $out = `$cmd`;
		#print "$out\n";
		$logger -> debug("Command to Merge: $cmd");
                $blastCOGFile{$sampleId} = "$combinedFile";
            }
    }
    return \%blastCOGFile;
}

sub get_parseBlastCOGCmd{
    my %seqFile = %{$_[0]};
    my %cdhitFile = %{$_[1]};
    my $suffix = $_[2];
    my @cmds;
    foreach my $sampleId ( keys %seqFile){
        my $file = $seqFile{$sampleId};
        if( -e $file){
            my $filename = basename($file);
            my $dirname = dirname($file);
            my $BLASTCOG_FILE = $file;
            my $CDHIT_FILE = $cdhitFile{$sampleId};
            my $PREPARSE_FILE = "$dirname/$sampleId\.preparse_COG";
            my $OUTPUT_FILE = "$dirname/$sampleId\.$suffix";
            my $tmp_cmd = "$Bin/Preparse_Cog_Blast1_Gautam.pl $dbLink/COG.mappings.v8.2.txt $CDHIT_FILE $BLASTCOG_FILE > $PREPARSE_FILE ;";
               $tmp_cmd .= " $Bin/Posparse_Cog_Blast1_Gautam.pl $PREPARSE_FILE > $OUTPUT_FILE";
            $annotationFile{$sampleId} = $OUTPUT_FILE;
            push(@cmds, $tmp_cmd);
        }
    }
    return \@cmds;
}

sub get_makeBlastdb{
    my $cdhitFile = %{$_[0]};
    my $splitFile = %{$_[1]};
    my @cmds;
    foreach my $sampleId ( keys %splitFile){
        if(exists $splitFile{$sampleId}){
            my @allfiles = @{$splitFile{$sampleId}};
            #foreach my $file (@allfiles){
             	my $CDHIT_FILE = $cdhitFile{$sampleId};
		        $logger -> debug("make blast database: input file => $CDHIT_FILE");
                my $tmp_cmd = "makeblastdb -in $CDHIT_FILE -dbtype nucl;";
                push(@cmds, $tmp_cmd);
            #}
        }
    }
    return \@cmds;
}


sub get_BlastForStitchCmd{
    my %annotationFile = %{$_[0]};
    my %splitFile = %{$_[1]};
    my $cdhitFile = %{$_[2]};
    my $suffix = $_[3];
    my @cmds;
    foreach my $sampleId ( keys %splitFile){
        if(exists $splitFile{$sampleId}){
            my @allfiles = @{$splitFile{$sampleId}};
            foreach my $file (@allfiles){
                my $filename = basename($file, ".fasta");
                my $dirname = dirname($file);
                my $INFILE = $file;
                my $CDHIT_FILE = $cdhitFile{$sampleId};
                my $OUTFILE = "$dirname/$filename\.$suffix";
		$logger -> debug("Blast_to_Stitch: inFile:outFile\t$INFILE\t$OUTFILE");
                #my $tmp_cmd = "makeblastdb -in $CDHIT_FILE -dbtype nucl;";
                my $tmp_cmd = " blastn -query $INFILE -db $CDHIT_FILE -outfmt 6 -evalue 1e-7 -dust no -out $OUTFILE";
                push(@cmds, $tmp_cmd);
            }
        }
    }
    return \@cmds;
}

sub get_MakeGraphCmd {
    my %seqFile = %{$_[0]};
    my %cdhitFile = %{$_[1]};
    my $suffix = $_[2];
    my @cmds;
    foreach my $sampleId ( keys %seqFile){
        my $file = $seqFile{$sampleId};
        if( -e $file){
            my $filename = basename($file);
            my $dirname = dirname($file);
            my $INFILE = $file;
            my $CDHITFILE = $cdhitFile{$sampleId};
            my $OUTFILE = "$dirname/${sampleId}_$suffix";
            my $tmp_cmd = "$Bin/MakeGraph_Gautam.pl $CDHITFILE $INFILE > $OUTFILE";
            $graphFile{$sampleId} = $OUTFILE;
            push(@cmds, $tmp_cmd);
        }
    }
    return \@cmds;
}

sub get_StitchContigsCmd {
    my %seqFile = %{$_[0]};
    my %cdhitFile = %{$_[1]};
    my $suffix = $_[2];
    my @cmds;
    foreach my $sampleId ( keys %seqFile){
        my $file = $seqFile{$sampleId};
        if( -e $file){
            my $filename = basename($file);
            my $dirname = dirname($file);
            my $INFILE = $file;
            my $CDHIT_FILE = $cdhitFile{$sampleId};
            my $ANNOTATION_FILE = $annotationFile{$sampleId};
            my $OUTFILE = "$dirname/${sampleId}\.$suffix";
            my $tmp_cmd = "$Bin/Stitch_Contigs_Gautam.pl $CDHIT_FILE $INFILE $ANNOTATION_FILE > $OUTFILE";
            $contigFile{$sampleId} = $OUTFILE;
            push(@cmds, $tmp_cmd);
        }
    }
    return \@cmds;
}

sub get_RenameContigsCmd {
    my %seqFile = %{$_[0]};
    my %cdhitFile = %{$_[1]};
    my $suffix = $_[2];
    my @cmds;
    foreach my $sampleId ( keys %seqFile){
        my $file = $seqFile{$sampleId};
        my $filename = basename($file);
        my $dirname = dirname($file);
        my $INFILE = $file;
        my $CDHITFILE = $cdhitFile{$sampleId};
        my $OUTFILE = "$dirname/${sampleId}\.$suffix";
        my $tmp_cmd = "$Bin/FR-Hit_RenameFinalContigs.pl $INFILE $CDHITFILE > $OUTFILE";
        $contigFile{$sampleId} = $OUTFILE;
        push(@cmds, $tmp_cmd);
    }
    return \@cmds;
}

sub get_AnnotateContigsCmd{
    my %seqFile = %{$_[0]};
    my $suffix = $_[1];
    my @cmds;
    foreach my $sampleId ( keys %seqFile){
        my $file = $seqFile{$sampleId};
        if(-e $file){
            my $filename = basename($file,".fasta");
            my $dirname = dirname($file);
            my $INFILE = $file;
	    my $OUTFILE = "$dirname/resfam_annotation";
            #my $OUTFILE = "$dirname/$filename\.$suffix";
	    my $tmp_cmd;
            if ($resAnnotate) {
	    	$tmp_cmd = "annotate_functional_selections.py -contigs $INFILE -m $resAnnotate --resfams --prefix $sampleId -o $OUTFILE -f";
            }
	    else {
		$tmp_cmd = "annotate_functional_selections.py -contigs $INFILE  --resfams --prefix $sampleId -o $OUTFILE -f";
	    }
	    push(@cmds, $tmp_cmd);
        }
    }
    return \@cmds;
}

sub get_Logfile_Name {
    return "$inputDir/parfums.log";
}

sub parse_command_line {
    my $help;

    if (!@ARGV) {
        print "$0: Argument required.\n";
        $logger -> logdie("Error in commmand: $!");
        exit 1;
    }

    my $result = GetOptions (
        "bcfile=s"          => \$bc_file,
        "dir=s"             => \$inputDir,
	"dblink=s"	    => \$dbLink,
        "paramPhrap1:s"     => \$paramPhrap1,
        "paramPhrap2:s"     => \$paramPhrap2,
        "help"              => \$help,
        "unchecked:s"       => \$sackedIDs,
	"resannotate:s"	    => \$resAnnotate
    ) or $logger -> logdie("Error in command line argument: $!");

    &usage()if ($help);

}

sub usage {

print<<EOF;
    PaRFuMs phrap_assembly, by Manish Boolchandani (manish\@wustl.edu),

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


