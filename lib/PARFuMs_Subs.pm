package PARFuMs_Subs;

#######################################################
#   Package: PARFuMs_Subs
#   Author: Manish Boolchandani
#   Date: 2015/02/28
#   Purpose: Harbors general functions used multiple
#   times like splitFasta, runCrossMatch
#
######################################################

use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename;
use Capture::Tiny ':all';
use lib $Bin;
use Log::Log4perl;
use Getopt::Long;
use POSIX;

our @EXPORT = qw();
our @EXPORT_OK = qw(split_FastaFiles merge_CrossMatch_Output getCM_cmd submit_jobArray frhit_clean);
our %EXPORT_TAGS = ( ALL => [@EXPORT_OK]);

##setup log file configuration path
#my $dirname = dirname($Bin);
#my $log_conf = "/opt/apps/labs/gdlab/software/parfums/1.0/config_files/log4perl.conf";
#Log::Log4perl::init_once($log_conf);
#my $log_subs = Log::Log4perl->get_logger("Perl_Module");
#$log_subs -> info("Inside Perl module");
#
#sub get_Logfile_Name {
#    return "/scratch/gdlab/manish/parfums/160127_testrun/01_testrun_MB/parfums.log";
#}

sub split_FastaFiles {
    my %filenames = %{$_[0]};
    my $outDir = $_[1];
    #print "tmpdir in splitFasta=> $tmpdir";
    my $MAXSIZE = $_[2];
    my $multi = $_[4];
    my $log_subs = $_[3];
    my @suffList = qw(.fasta .clean .fna);
    my %spFasta;
    $log_subs -> debug("splitFasta inside perl module");
    if ($multi) {
        %spFasta = &split_MultiLineFasta(\%filenames, $outDir, $MAXSIZE, $log_subs);
    }
    else {
        foreach my $id (keys %filenames) {
            my $file = $filenames{$id};
            my @splitFasta;
            my $iterator = 0; my $filepart = 1;
            open (my $fastaIN, "< $file") or die "$!";
            open (my $OUT, "> $outDir/$id-part$filepart.fasta") or die "$!";
            my ($readID, $seq);
            while(my $line = <$fastaIN>) {
                chomp($line);
                if ($line =~ /^>/) {
                    $readID = $line;
                    chomp($seq = <$fastaIN>);
                    if ($iterator < $MAXSIZE) {
                        print $OUT "$readID\n$seq\n";
                        $iterator++;
                    }
                    else {
                        close($OUT);
                        push(@splitFasta, "$outDir/$id-part$filepart.fasta");
                        #$log_subs -> debug("Created $id-part$filepart.fasta");
                        $filepart++;
                        $iterator = 0;
                        open ($OUT, "> $outDir/$id-part$filepart.fasta") or die "$!";
                        print $OUT "$readID\n$seq\n";
                        $iterator++;
                    }
                }
            }
            push(@splitFasta, "$outDir/$id-part$filepart.fasta");
            #$log_subs -> debug("Created $id-part$filepart.fasta");
            $spFasta{$id} = [ @splitFasta ];
        }
    }
    return %spFasta;
}

sub split_MultiLineFasta {
    my %filenames = %{$_[0]};
    my $outdir = $_[1];
    my $MAXSIZE = $_[2];
    my $log_subs = $_[3];
    my %spFasta;
    $log_subs -> debug("splitDebug: outdir:$outdir and MAX:$MAXSIZE");
    foreach my $id (keys %filenames) {
        my $file = $filenames{$id};
	$log_subs ->debug("Reading file in Multi-Line split: $file");
        open (my $fastaIN, "< $file") or die "$!";
        my $readID; my %fastaSeq;
        while(my $line = <$fastaIN>) {
            chomp($line);
            if ($line =~ /^>/) {
                $readID = $line;
                $readID =~ s/>//;
                $fastaSeq{$readID} = "";
                next;
            }
            $fastaSeq{$readID} = $fastaSeq{$readID}.$line;
        }
	close($fastaIN);
	my @allIDs = keys(%fastaSeq);
	$log_subs ->debug("IDs: @allIDs");
        my $iterator = 0; my $filepart = 0;
        open (my $OUT, "> $outdir/$id-part$filepart.fasta");
        my @splitFasta;
        foreach my $readID (keys %fastaSeq) {
            if ($iterator < $MAXSIZE) {
                print $OUT ">$readID\n$fastaSeq{$readID}\n";
                $iterator++;
            }
            else {
                close($OUT);
                push(@splitFasta, "$outdir/$id-part$filepart.fasta");
                #$log_subs -> debug("Created $id-part$filepart.fasta");
                $filepart++; $iterator = 0;
                open ($OUT, "> $outdir/$id-part$filepart.fasta");
                print $OUT ">$readID\n$fastaSeq{$readID}\n";
            }
        }
        push(@splitFasta, "$outdir/$id-part$filepart.fasta");
        #$log_subs -> debug("Created $id-part$filepart.fasta");
        $spFasta{$id} = [ @splitFasta ];
    }
    return %spFasta;
}

sub submit_jobArray_SGE {
    my $jobscript_file = $_[0];
    my $jobName = $_[1];
    my $workDir = $_[2];
    #$log_subs -> info("Running $jobName on cluster using qsub");
    if (not(-d $workDir)) {
        mkdir $workDir;
    }
    chdir($workDir);
    my $cmd = "nq $jobscript_file | qsub -P long -N $jobName -sync y";
    my $job_return = `$cmd`;
    #$log_subs -> info("$job_return");
    my $jobID = $1 if ($job_return =~ /^Your job-array (\d+)\./);
    #$log_subs -> info("job $jobID completed");
}

sub submit_jobArray {
    my $job_file = $_[0];
    my $job_name = $_[1];
    my $work_dir = $_[2];
    my $log_subs = $_[3];
    
    my $step_size = 1;
    #$log_subs -> debug("Running $job_name on cluster using slurm\n");
    
    if (not(-d $work_dir)) {
        mkdir $work_dir;
    }
    
    chdir($work_dir);
    my $main_dir = dirname($Bin);
    
    my $linecnt_cmd = "wc -l $job_file | awk '{print \$1}'";
    #$log_subs -> debug("cmd => $linecnt_cmd"); 
    
    my $linecnt = capture{system("$linecnt_cmd");};
    
    chomp($linecnt);
    #$log_subs -> debug("linecnt => $linecnt");
    
    if ($linecnt > 5000){	
    	$step_size = ceil(int($linecnt) / 5000);
    }
    $log_subs -> debug("step size = $step_size"); 

    my $cmd = "$main_dir/src/submit_job_array_slurm.sh $job_file $step_size | sbatch -J $job_name --mem=30000M";
    $log_subs -> debug("cmd => $cmd");
    
    my $job_return = capture{ system("$cmd");};
    $log_subs -> debug("job_return => $job_return");
    
    my $jobID = $1 if ($job_return =~ /^Submitted batch job (\d+)/);
    my $file = "$job_name\_$jobID\_check.txt"; 
    
    if (-e $file){
	unlink $file or die("could not delete $file");
    }
    
    my $job_cmd = "sbatch --dependency=singleton --job-name=$job_name $main_dir/src/complete.sh $work_dir\/$file | awk \'{print \$NF}\'";
    #$log_subs -> debug("job_cmd => $job_cmd");
    $job_return = capture{ system("$job_cmd");};
    #$log_subs -> debug("dependent job id = $job_return");
    &wait_job_to_finish($jobID, $file, $log_subs, $linecnt);
    $log_subs -> debug("job $jobID completed");
}

sub wait_job_to_finish{
	
	my $job_id = $_[0];
 	my $flag = 1;
	my $completed_file = $_[1];
	my $log_subs = $_[2];
	my $job_cnt = $_[3]; 
	
	my $status = 0;
	#$log_subs ->debug("$job_id inside function");
	#$log_subs -> debug("dependent_file => $completed_file");
	while($flag){
		if(-s $completed_file){
			$log_subs -> debug("job $job_id finished..Will now check for errors");
			$flag = 0;
			my $success_cnt = 0;
			my $err_cnt = 0;
			my @files = glob("slurm-$job_id\_*");
			#$log_subs -> debug(@files);
			foreach my $file (@files){
				#$log_subs -> debug("slurm filename=>$file");
				my $cmd = "grep -h -c \'^Well done!! Job finished\' $file";
				my $err_cmd = "grep -h -c \'error\' $file";		
				my ($stdout, $stderr, $exit) = capture {system("$cmd");};
				#$log_subs -> debug("filename=> $file stdout=> $stdout stderr=> $stderr exit=> $exit");
				if (!$stderr){
					chomp($stdout);
					#$log_subs -> debug("filename success count => $stdout");
					$success_cnt += $stdout;
				}

				($stdout, $stderr, $exit) = capture{ system("$err_cmd");};
				if (!$stderr){
					chomp($stdout);
					#$log_subs -> debug("filename error count => $stdout");
					$err_cnt += $stdout;
				}
				
			}
			
			 $log_subs -> debug("success_cnt => $success_cnt");
			 $log_subs -> debug("err_cnt => $err_cnt");
			
			if($success_cnt == $job_cnt && $err_cnt == 0){
				$log_subs -> debug("job $job_id completed successfully");
				
			}
			else{
				$log_subs -> debug("Error occurred while running the job $job_id");
				exit 1;
			}
		}

		sleep(5);

	}
	return ;
}
	
sub wait_job_to_finish_old {	
	my $jobid = shift;
	my $log_subs = shift;
	$log_subs -> info("Waiting for job: $jobid to finish!");
	#print "jobid $jobid inside function\n";
	my $flag = 1; 
	my $counter = 1;
	my $prev_status = "";
	while ($flag){
		sleep(5);
		#Check job status
		#my $status = `squeue -j $jobid`;
		my $st_check_cmd = "squeue -j $jobid -h -o %t";
		my ($stdout, $stderr, $exit) = capture { system("$st_check_cmd"); };
		if ($stderr =~ /timed out/) {
			sleep(5);
			($stdout, $stderr, $exit) = capture { system("$st_check_cmd"); };
		}
		chomp($stdout);
		#$log_subs -> debug("current_status : $stdout,$stderr,$exit");
		
		if ($stdout =~ /R/) {
			if($counter || $stdout ne $prev_status){
				$log_subs -> debug("Job $jobid is running");
				$counter = 0;
				$prev_status = $stdout;
			}
		}
		elsif ($stdout =~ /CA|F|NF/) {
			$log_subs -> debug("Job $jobid failed. Exiting!!");
		}
		elsif ($stdout =~ /PD/) {
			if ($counter || $stdout ne $prev_status){
				$log_subs -> debug("Job $jobid is pending");
				$counter = 0;
				$prev_status = $stdout;
			}
		} 
		elsif ($stdout eq "CD" || $stdout eq "") {
			sleep(5);	
			$log_subs -> info("Checking whether all jobs completed successfully");
			my $job_detail = `sacct -P -n -j $jobid --format jobid,jobname,elapsed,ReqMem,alloccpus,state,exitcode`;
			#$log_subs -> debug("job_detail = $job_detail");
			my @lines = split /\n/, $job_detail;
			#print $lines[0];
			foreach my $line(@lines) {
				my @data = split /\|/, $line;
				my $msg = $data[5];
				#get exit codes for the job in an array
				my $exit_code = $data[6];
				chomp($exit_code);
				if($exit_code ne "0:0" || $msg ne "COMPLETED"){
					$log_subs -> info("Job: $jobid terminated with an error");
					#print "Check parfums.log file for details\n";
					exit 1;	
				}
			}
			$flag = 0;
			#exit 1;	
		}
		else {
			$log_subs -> debug("Job $jobid state: $stdout")	
		} 		
	} 
}

sub merge_CrossMatch_Output {
    my %CMfiles = %{$_[0]};
    my $inputDir = $_[1];
    my $fileCount = 0;
    foreach my $seqID (keys %CMfiles) {
        my @allfiles = @{$CMfiles{$seqID}};
        if (@allfiles > 1) {
            my $name = basename($allfiles[0], ".cmoutput");
            my $combFile = "$inputDir/$seqID/$name.cmoutput";
            $combFile =~ s/-part\d+//g;
            my @cmout_files = @allfiles;
            #my @cmout_files = map {$_ = $_ . ".cmoutput" } @allfiles;
            my @flagFiles = map { -e $_ ? 1 : 0 } @cmout_files;
            if (grep {$_ != 1} @flagFiles) {
                die "Cross Match not complete: Please check error log\n";
                #$log_subs -> logdie("CrossMatch not complete: Please check error log");
            }
            else {
		if (-e $combFile) { unlink $combFile or die "Cannot delete pre-existing $combFile"; }
                foreach my $cmfile (@cmout_files) {
                    my $cmd = "cat $cmfile >> $combFile";
                    system("$cmd");
                }
                #$log_subs -> info("$seqID merged file formed");
                $fileCount++;
            }
        }
        else {
            system ("cp", "$CMfiles{$seqID}[0]", "$inputDir/$seqID/");
            $fileCount++;
        }
    }
    #$log_subs ->info("Total $fileCount files processed");
}

sub getCM_cmd {
    my $CM_cmd = "cross_match FASTA_FILE_1 FASTA_FILE_2 -gap1_only -minmatch 6 -minscore 10 -gap_init -3 | perl CODE_TO_REMOVE FASTA_FILE_1";
    return $CM_cmd;
}

sub get_frhit_cmd {
    my $frhit_cmd = "fr-hit -a FASTA_FILE_1 -d FASTA_FILE_2 -o OUTPUT_FILE -c 70 -m 40 -r 0";
    return $frhit_cmd;
}

sub frhit_clean {
    my %frhit_out = %{$_[0]};
    my $log_subs = $_[1];
    #my $fastafile = shift;
    my %matchReads;
    my $id_count = scalar(keys %frhit_out);
    foreach my $id (keys %frhit_out) {
        my $frhit_file = $frhit_out{$id};
        if (!-z $frhit_file) {
            $log_subs -> debug("$frhit_file found. Processing Information");
            open (my $in, "< $frhit_file") or die "$!";
            while(my $line = <$in>) {
                my @vals = split(/\s+/, $line);
                my $readID = $vals[0];
                $readID =~ s/\_\d$//;
                if (exists $matchReads{$id}{$readID}) {
                    next;
                }
                else {
                    $matchReads{$id}{$readID} = 1;
                }
            }
        }
        else {
            #$log_subs->info("$frhit_file is empty");
            $matchReads{$id}{"No_Match_Reads"} = 1;
        }
    }
    my $matchIDs = scalar(keys %matchReads);
    $log_subs -> debug("$matchIDs in matchReads");
    return \%matchReads;
}
