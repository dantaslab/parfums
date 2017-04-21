#!/usr/bin/perl -w

####################################################################################################
# Created by:           Manish Boolchandani
# Creation Date:        2016/05/12
# Version:              1.0
# Description:          annotation of contigs
#
# Modify By:		Manish Boolchandani
# Modification Date:	2016/05/12
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
use File::Basename;
use Log::Log4perl;
use Cwd;
use PARFuMs_Subs qw(split_FastaFiles submit_jobArray get_FRhit_cmd get_CDhit_cmd);

my(%inputSeqFile, %phrapFile, %contigFile, %splitFile, %totalReads, %annotationFile , %blastCOGFile, %cdhitFile, %frhitFile, %mapFile, %graphFile, %cogFile);
my(@params, @uncheckedIDs);
my($bc_file, $inputDir, $runResfam, $jobFile, $sackedIDs, $resAnnotate,$dbLink);

#------------------------------------------------------------------------#
# START OF MAIN PROGRAM                                                  #
#------------------------------------------------------------------------#
&parse_command_line();

#setup log file configuration path
my $maindir = dirname($Bin);
my $log_conf = "$maindir/config_files/log4perl.conf";
my $resfam_conf = "$maindir/config_files/fmg_config.txt";

Log::Log4perl::init_once($log_conf);
my $logger = Log::Log4perl -> get_logger("Annotation");

$logger -> info("Step 6: Started PIPELINE FOR ANNOTATION OF CONTIGS");

my $MAXSIZE = 5;
my $counter = 0;

if ($sackedIDs) {
    @uncheckedIDs = split(/\,/, $sackedIDs);
}

if (!$runResfam){
   $runResfam = 1;
}

$logger -> info("Reading last-contig files");
%inputSeqFile = &read_inputSequences("last-contigs.fasta");

#-----------------------------------------------------------------#
# Annotation of Contigs                                           #
#-----------------------------------------------------------------#

my $template = "Annotate_XXXX";
my $tempDir = tempdir( $template, DIR => "$inputDir/tmp", CLEANUP => 0 );

#Change the permissions of the temporary directory in order to make it accessible by everyone
if( -d $tempDir ) {
        $logger -> debug("$tempDir exists");
        chmod(0755, $tempDir) or $logger -> die("Unable to change the permissions of dir $tempDir");
} else {
        $logger -> debug("$tempDir does not exist");
}

$logger -> info("Running Blast against COG database");

##### Blast against COG to do functional annotation
&split_fastafile("last-contigs.fasta", $MAXSIZE, 0);
$jobFile = &generate_script("blastCOG", \%splitFile, "blastCOGFinal");
PARFuMs_Subs::submit_jobArray($jobFile, "blastCOGFinal2", $tempDir, $logger);

$logger -> info("Merging BlastCOG output");

my $cogRef = &merge_BlastOutput(\%splitFile, "blastCOGFinal");
$logger -> info("blastCOG final and Merging output completed");
%cogFile = %$cogRef;

#New set of code added for annotation-final
$jobFile = &generate_script("parse-blast-final", \%cogFile, "annotation-final.txt");
PARFuMs_Subs::submit_jobArray($jobFile, "ParseBlastCOG", $tempDir, $logger);
$logger -> info("parse-blast-final run completed");

$logger -> info("Using Molly Gibson's script + HMMER3 + HMM database");
##### Molly Gibson's Python-based GeneFinder + HMM-annotater, supplementing COG-based annotation #####
$jobFile = &generate_script("annotate", \%inputSeqFile, "annotation.txt");
PARFuMs_Subs::submit_jobArray($jobFile, "Annotate", $tempDir, $logger);
$logger -> info("Annotation completed!!");

#### Now do some cleaning. In particular big Blast files
#system("cd $maindir/; rm -r *.nhr *.nin *.nsq *.clstr");


#------------------------------------------------------------------------#
# SUBROUTINES                                                            #
#------------------------------------------------------------------------#

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
	
    $logger -> debug("Writing seqFile hash");
    foreach my $k (keys %seqFile) {
        #print "$k\t$seqFile{$k}\n";
	$logger -> debug("$k => $seqFile{$k}");
    }
    
    return %seqFile;
}

sub split_fastafile{

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
    my $jobfile = "${inputDir}/tmp/qsub_scripts/${runType}_script_a${counter}.sh";
    my $jobfilename = basename($jobfile);
    my @cmds;
    my $param = $runResfam;
    if(scalar(@_) == 3){
        $suffix =$_[2];
    }

    open (my $OUT, "> $jobfile") or $logger -> logdie($!);
    $logger -> info("Generating a script file to run $runType: $jobfilename");

    if($runType =~ "blastCOG"){
        $logger -> info("Writing script file to blast against COG");
        my $cmdRef = get_BlastCOGCmd(\%seqFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType eq 'parse-blast'){
        $logger -> info("Writing script to parse blast COG output");
        my $cmdRef = get_parseBlastCOGCmd(\%seqFile, \%cdhitFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType =~ "annotate"){
        $logger -> info("Writing script file to annotate contigs");
        my $cmdRef = get_AnnotateContigsCmd(\%seqFile, $suffix);
        @cmds = @$cmdRef;
    }
    elsif($runType eq 'parse-blast-final') {
    	$logger -> info("Writing script to parse blast COG Final output");
	my $cmdRef = get_parseBlastCOGCmd(\%seqFile, \%inputSeqFile, $suffix);
	@cmds = @$cmdRef;
    }
    else {
        $logger -> logdie("Invalid runType: $runType");
    }

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
	$logger -> debug("sampleid => $sampleId");
        my $file = $seqFile{$sampleId};
        if( -e $file){
            my $filename = basename($file);
            my $dirname = dirname($file);
            my $BLASTCOG_FILE = $file;
            my $CDHIT_FILE = $cdhitFile{$sampleId};
            my $PREPARSE_FILE = "$dirname/$sampleId\.preparse_COG";
            my $OUTPUT_FILE = "$dirname/$sampleId\.$suffix";
	    $logger -> debug("filename => $filename cdhit_file => $CDHIT_FILE preparse_file => $PREPARSE_FILE");
            my $tmp_cmd = "$Bin/Preparse_Cog_Blast1_Gautam.pl $dbLink/COG.mappings.v8.2.txt $CDHIT_FILE $BLASTCOG_FILE > $PREPARSE_FILE ;";
               $tmp_cmd .= " $Bin/Posparse_Cog_Blast1_Gautam.pl $PREPARSE_FILE > $OUTPUT_FILE";
            $annotationFile{$sampleId} = $OUTPUT_FILE;
            push(@cmds, $tmp_cmd);
        }
    }
    return \@cmds;
}

sub get_AnnotateContigsCmd{
    my %seqFile = %{$_[0]}; 	#seqFile{sampleID} = Last contig filename
    my $suffix = $_[1];			
    my @cmds;

    my $OUTDIR = "$inputDir/resfam_annotation";
    if(! -d $OUTDIR){
	system("mkdir $OUTDIR");
    }
    my $OUTFILE = "$OUTDIR/merged_contigs.fasta";
    if ( -e $OUTFILE ) {
        unlink($OUTFILE) or $logger -> logdie("$OUTFILE: $!");
    }
    open (my $OUT, "> $OUTFILE") or $logger -> logdie("$!");
    foreach my $sampleId ( keys %seqFile){
        my $file = $seqFile{$sampleId};
        if(-e $file){
	    open (my $fastaIN, "< $file") or $logger -> logdie("$!");
	    while(my $line = <$fastaIN>){
		chomp($line);
	    	if ($line =~/^>/){
			my @headerArr = split(/>/, $line);
			my $headerline = ">". $sampleId . "_" . $headerArr[1] ;
			print $OUT "$headerline\n";   
		}
		else {
			print $OUT "$line\n";
		}
	    }
	    close($fastaIN);
            #my $filename = basename($file,".fasta");
            #my $dirname = dirname($file);
            #my $INFILE = $file;
            #my $OUTFILE = "$dirname/$filename\.$suffix";
            
        }
    }
    close($OUT);
    
    my $tmp_cmd;
    if ($resAnnotate) {
       $tmp_cmd = "$Bin/annotate_parfum_contigs.py -contigs $OUTFILE -config $resfam_conf -m $resAnnotate --resfams -o $OUTDIR --prefix merged_contigs -f";
       #srun /scratch/gdlab/manish/resistome_collab/src/annotate_functional_selections.py -contigs /scratch/gdlab/manish/parfums/demo/06_testrun_MB/resfam_annotation/merged_contigs.fasta -m /scratch/gdlab/manish/parfums/demo/02_annotation_mappingfile.txt --resfams  -o /scratch/gdlab/manish/parfums/demo/06_testrun_MB/resfam_annotation --prefix mergedContigs -f
       $logger -> debug("annotation command => $tmp_cmd")
    }
    else {
       $tmp_cmd = "$Bin/annotate_parfum_contigs.py -contigs $OUTFILE -config $resfam_conf --resfams -o $OUTDIR --prefix merged_contigs -f";
    }
    
    push(@cmds, $tmp_cmd);
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
        "runResfam:s"       => \$runResfam,
        "help"              => \$help,
        "unchecked:s"       => \$sackedIDs,
	"resannotate:s"	    => \$resAnnotate
    ) or $logger -> logdie("Error in command line argument: $!");

    &usage()if ($help);

}

sub usage {

print<<EOF;
    PaRFuMs annotate_v1.pl, by Manish Boolchandani (manish\@wustl.edu),

    This program reads *.last-contigs.fasta file and splits it into several
    smaller files, based on seqIds in the barcode file and then blast against COG
    database and run ResFam annotation
    each file.

    usage: $0 --bcfile FILE --dir OUTPUT_DIR [--help]

    Arguments:

    --bcfile FILE       - Barcodes file name. (Barcodes and Identifiers separated by tab)
    --resannotate FILE	- Annotation mapping file (id, library, and anitbiotic separated by tab )
    --dir OUTPUT_DIR    - Directory in which all output files will be saved.
    --help              - This helpful help screen.

EOF
exit 1;
}


