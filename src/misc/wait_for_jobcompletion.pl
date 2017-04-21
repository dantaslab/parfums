#!/usr/bin/perl -w

use strict;
use lib "/scratch/gdlab/manish/parfums/parfum_testing/1.0/lib";
use FindBin qw($Bin);
use Capture::Tiny ':all';
use File::Basename;
use File::Find;
print("This is a test\n");
my $jobscript_file = "blastCOG_script_p6.sh";
#my $jobscript_file = "cd-hit-est_script_p1.sh";
my $jobName = "test_job";
my $workDir = "/scratch/gdlab/manish/parfums/temp";
&submit_jobArray($jobscript_file, $jobName, $workDir);

sub submit_jobArray {
    my $job_file = $_[0];
    my $job_name = $_[1];
    my $work_dir = $_[2];
    #my $log_subs = $_[3];
    print("Running $job_name on cluster using slurm\n");
    if (not(-d $work_dir)) {
        mkdir $work_dir;
    }
    chdir($work_dir);
    #if ( $jobName eq "Annotate"){
    #    system("cp $Bin/lib/.gm_key $workDir");
    #}
    #my $script = dirname($Bin);
    my $cmd = "$Bin/src/submit_job_array_slurm.sh $job_file 1 | sbatch -J $job_name --mem=30000";
    print("cmd => $cmd\n");
    my $job_return = `$cmd`;
    print("job_return => $job_return\n");
    my $jobID = $1 if ($job_return =~ /^Submitted batch job (\d+)/);
    my $file = "test.txt"; 
    if (-e $file){
	unlink $file or die("could not delete $file");
    }
    #$job_return = `sbatch --dependency=afterok:$jobID complete.sh | awk \'{print \$NF}\'`;
    $job_return = `sbatch --dependency=singleton --job-name=$job_name complete.sh $jobID | awk \'{print \$NF}\'`;
    print("dependent job id = $job_return");
    &wait_job_to_finish($jobID);
    print("job $jobID completed\n");
}

sub wait_job_to_finish{
	my $job_id = $_[0];
 	my $flag = 1;
	my $completed_file = "test.txt";
	#print @files;
	my $status = 0;
	my $total_cnt = 0;
	print("$job_id inside function\n");
	while($flag){
		if(-s $completed_file){
			print("job $job_id finished\nWill now check for errors\n");
			$flag = 0;
			$total_cnt = 0;
			my @files = glob("slurm-$job_id\_*");
			foreach my $file (@files){
				print("slurm filename=>$file\n");
				my $cmd = "grep -h -c \'Finished\' $file";		
				my ($stdout, $stderr, $exit) = capture {system("$cmd");};
				print("filename=> $file stdout=> $stdout stderr=> $stderr exit=> $exit \n");
				if (!$stderr){
					chomp($stdout);
					print ("filename error count => $stdout\n");
					$total_cnt += $stdout;
				}
			}
			
			print "total_cnt => $total_cnt\n";
			
			if($total_cnt == scalar(@files)){
				print ("job $job_id completed successfully\n");
				exit 0;
			}
			else{
				print("Error occurred while running the job $job_id\n");
				exit 1;
			}
		}

		sleep(5);

	}
	
 exit 0;
}
