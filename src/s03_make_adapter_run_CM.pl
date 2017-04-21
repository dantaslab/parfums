#!/usr/bin/perl -w
################################
#   Author: Manish Boolchandani
#   Date: 2015/02/18
#   Purpose: Make adapter file based on barcodes
###############################

use strict;
use FindBin qw($Bin);
use File::Basename;
use File::Temp qw(tempdir tempfile);
use File::Spec::Functions qw(catdir);

use lib catdir(dirname($Bin), 'lib');
use Getopt::Long;
use Log::Log4perl;
use Cwd;
use PARFuMs_Subs qw(split_FastaFiles submit_jobArray);

#Global variables
my $MAXSIZE = 200000;
my ($adapA, $adapB, $endVec1, $endVec1R, $endVec2, $endVec2R);
my ($bc_file, $inputDir, $sackedIDs);
my (%spFasta, %totalReads, %cm_outfile);
my @uncheckedIDs;


#START OF THE PROGRAM
&parse_command_line();

#setup log file configuration path
my $dirname = dirname($Bin);
my $log_conf = "$dirname/config_files/log4perl.conf";
Log::Log4perl::init_once($log_conf);
my $logger = Log::Log4perl->get_logger("Adapter_CM_Run");
$logger ->info("STEP 2: REMOVING ADAPTER SEQUENCES");
$logger -> debug("Current directory:$dirname");
print("\nStep 2: Making adapter files and CrossMatch run\n");

if ($sackedIDs) {
	@uncheckedIDs = split(/\,/, $sackedIDs);
}
#Defined Variables
my $template = "adap_CMRun_XXXX";
my $tempDir = tempdir($template, DIR => "$inputDir/tmp/", CLEANUP => 0);
if( -d $tempDir ) {
	$logger -> debug("$tempDir exists");
	chmod(0755, $tempDir) or $logger -> die("Unable to change the permissions of dir $tempDir");
} else {
	$logger -> debug("$tempDir does not exist");
}

$logger -> debug("Writing adapter file");
&update_adapters();

$logger -> debug("Split Fasta Files for CrossMatch run");
&read_FastaFiles();
$logger -> debug("Making jobscript file");
my $cm_jobfile = &make_CrossMatch_Jobscript();
PARFuMs_Subs::submit_jobArray($cm_jobfile, "adap_CM_run", $tempDir, $logger);

$logger -> debug("Merging clean fasta output files");
foreach my $id (keys %cm_outfile) {
	my $clean_filename = "$inputDir/$id/$id.fasta.clean";
	if (-e $clean_filename) {
		unlink $clean_filename or $logger ->die("Unable to remove pre-existing $clean_filename");
	}
	foreach my $filename (@{$cm_outfile{$id}}) {
		system("cat $filename >> $clean_filename");
	}
	$logger ->info("$clean_filename is formed");
}
#PARFuMs_Subs::merge_CrossMatch_Output(\%cm_outfile, $inputDir);

exit 0;
#END OF THE PROGRAM

sub make_CrossMatch_Jobscript {
    my $cm_jobfile = "$tempDir/CM_run_jobscript.sh";
    $logger -> debug("Generating a script file to run cross_match:CM_run_jobscript.sh");
    open (my $OUT, "> $cm_jobfile") or $logger -> logdie("$!");
    my $codepath = "$dirname/src/CleanAdapter-Illumina_PE_mod.pl";
    my $CM_cmd = PARFuMs_Subs::getCM_cmd();
    foreach my $id (sort {$totalReads{$a} <=> $totalReads{$b}} keys %totalReads) {
        if (not grep {$_ eq $id} @uncheckedIDs) { 
		if (exists $spFasta{$id}) {
            		my @files = @{$spFasta{$id}};
            		foreach my $filename (@files) {
                		my $name = fileparse($filename, ".fasta");
                		my $tmp_cmd = $CM_cmd;
                		my $adapFile = "$inputDir/$id/Adapters_$id.fna";
                		my $clean_outFile = "$tempDir/splitFastaFiles/$name.fasta.clean";
                		$tmp_cmd =~ s/FASTA_FILE_1/$filename/g;
                		$tmp_cmd =~ s/FASTA_FILE_2/$adapFile/g;
                		$tmp_cmd =~ s/CODE_TO_REMOVE/$codepath/g;
                		#my $cm_cmd = "cross_match $filename $inputDir/$id/Adapters_$id.fna" . $cm_param . "> $filename.cmoutput";
                		print $OUT "$tmp_cmd\n";
                		if (exists $cm_outfile{$id}) {
                		    push(@{$cm_outfile{$id}}, $clean_outFile);
                		}
                		else {
                		    $cm_outfile{$id} = [ $clean_outFile ];
                		}
			}
            	}
       		 else {
            		$logger->info("$id does not exist in splitFasta");
        	}
	}
	else {
		$logger ->info("$id is not included in analysis and is sacked");
    	}
    }
    close($OUT);
    return $cm_jobfile;
}

sub read_FastaFiles {
    my %processFiles;
    my $splitDir = "$tempDir/splitFastaFiles";
    if (not -d $splitDir) {
        mkdir $splitDir;
    }
    $logger -> debug("Value of MaxSize inside make adapter run: $MAXSIZE");
    open (my $in, "< $bc_file") or die "Could not find barcode file at $bc_file: $!";
    chomp(my @BC = <$in>);
    foreach my $line (@BC) {
        my @seqID = split(/\s+/, $line);
        my $fastaFile = "$inputDir/$seqID[1]/$seqID[1].fasta";
        if (-e $fastaFile) {
            chomp(my $readCount = `grep '^>' $fastaFile | wc -l`);
            $totalReads{$seqID[1]} = $readCount;
            
	    $logger -> debug("$seqID[1] --> $readCount");
            
	    if ($readCount > $MAXSIZE) {
                $processFiles{$seqID[1]} = $fastaFile;
            }
            else {
                system("cp $fastaFile $splitDir");
                $spFasta{$seqID[1]} = [ "$splitDir/$seqID[1].fasta" ];
            }
        }
    }
    my %tmpFasta = PARFuMs_Subs::split_FastaFiles(\%processFiles, $splitDir, $MAXSIZE, $logger);
    @spFasta{keys %tmpFasta} = values %tmpFasta;
}

sub update_adapters {
    open (my $bcIN, "< $bc_file") or die "$!";
    chomp(my @bcData = <$bcIN>);
    foreach my $bcline (@bcData) {
        my @bcInfo = split(/\s+/, $bcline);
        my $revBarcode = reverse($bcInfo[0]);
        $revBarcode =~ tr/ACGTacgt/TGCAtgca/;

        #Concatenating the reverse barcode in adapter sequences
        my $tmp_adapA = $revBarcode.$adapA;
        my $tmp_adapB = $revBarcode.$adapB;

        #Replacing X's and Y's from reverse and forward barcode respectively
        my $tmp_endVec1 = $endVec1; $tmp_endVec1 =~ s/XXXX/$revBarcode/;
        my $tmp_endVec1R = $endVec1R; $tmp_endVec1R =~ s/YYYY/$bcInfo[0]/ ;
        my $tmp_endVec2 = $endVec2; $tmp_endVec2 =~ s/XXXX/$revBarcode/ ;
        my $tmp_endVec2R = $endVec2R; $tmp_endVec2R =~ s/YYYY/$bcInfo[0]/ ;

        #Writing these barcodes in file
        $logger -> debug("Adapters_$bcInfo[1].fna");

        open (my $out, "> $inputDir/$bcInfo[1]/Adapters_$bcInfo[1].fna") or $logger -> logdie("$!");
        print $out ">Vector1\n$tmp_endVec1\n>Vector1R\n$tmp_endVec1R\n>Vector2\n$tmp_endVec2\n>Vector2R\n$tmp_endVec2R\n";
        print $out ">Adapter1\n$tmp_adapA\n>Adapter2\n$tmp_adapB\n";
        close($out);
    }
}

sub get_Logfile_Name {
    return "$inputDir/parfums.log";
}

sub parse_command_line {
    my $help;
    my $result = GetOptions ("bcfile=s" => \$bc_file,
                    "dir=s" => \$inputDir,
                    "adapA=s" => \$adapA,
                    "adapB=s" => \$adapB,
                    "endVec1=s" => \$endVec1,
                    "endVec1R=s" => \$endVec1R,
                    "endVec2=s" => \$endVec2,
                    "endVec2R=s" => \$endVec2R,
		    "unchecked=s" => \$sackedIDs,
                    "help" => \$help
            );
}
