##############################################################
#                                                            #
# How To Run PARFuMS (v1.1.0, newer version) on HTCF Cluster   #
#                                                            #
##############################################################

Note: PARFuMS pipleline was modified and recoded in order to make it work on the HTCF cluster. All the parameters are kept same as the previous version. 

==================================
1. Requirements: 
==================================
	
	1.1 Load parfum module
	------------------------

	In order to run PARFuMs on HTCF cluster, you need to first load parfums module. Currently, the default version is 1.0. 

	Note: As we improve the pipeline, the default version may change. Please make sure, you have loaded the correct parfums version.

		<code>module load parfums</code>

	The above command will load parfums (default v1.0) and all other dependencies required to run parfums. You can list all modules loaded using the following command:

 		<code> module list</code>

		Output:
		Currently Loaded Modules:
  		1) cd-hit/4.6.4   3) ncbi-blast/2.2.31+   5) phrap/1.090518   7) hmmer/3.1b1         9) pandas/0.15.2  11) scikit-learn/0.16.1  13) numpy/1.9.1
  		2) fr-hit/0.7.1   4) velvet/1.2.10        6) parfums/1.0      8) metagenemark/2.10  10) scipy/0.15.1   12) biopython/1.65 

	
	1.2. Input files required to run the pipeline:
	--------------------------------------------

	A) Barcode Mapping file: Create a tab delimited '.txt' file listing each of your samples and the associated barcode. The file should consists of a tab-delimited columns: Barcode and Sample Name. Here is the location of an example file here: "./examples/01_testrun/01_barcode_mappingfile.txt"
	
	B) Annotation Mapping File]]: Create tab-delimited file with columns: id, library, antibiotic. This file is required to run ResFams. Here is the location of an example file here: "./examples/01_testrun/02_annotation_mappingfile.txt"
	
	C) Config File : Config file is used to define multiple parameters used in PARFuMS run. Each line starts with the variable name that is used "as-is" in the code, followed by values stored in these variables. Any changes to the variable names might lead to error. Variables defined by this method include: Forward input file (FW_file), Reverse Input File (RC_file), Barcode file (BC_file), Output Directory (DIR), Barcode Position information (bol/eol), Adapter & Vector Sequence Information etc. Here is the location of an example file: "./examples/01_testrun/03_config_file.txt"
 
	Note: Most of this information remains unchanged during different runs. You might only need to change these variables FW_file, RC_file, BC_file, resAnnotate path in most cases.

===================================
2. Usage: 
===================================

	perl parfum_wrapper.pl --config [CONFIG_FILE] --step [1:6|0|]  

 	PARFUMS version 1.0, by Manish Boolchandani (manish@wustl.edu),
 	This program is a wrapper to run all steps of the PARFUMS pipeline.

 	Usage: perl parfum_wrapper.pl --config [CONFIG_FILE] --step [1:6|0|N1:N2];

 	Please find format of CONFIG_FILE in README.
 	Each step corresponds to a specific stage on PARFUMS pipeline.
 	Step information is entered to run a complete pipeline (enter 0), run a specific step (any number 1 to 6), or specific set of steps (2:5 to run Step 2,3,4,5).

 	Description of each step:
	-------------------------
 	Step 1: Preprocessing input FastQ files
			Reads FW- and RC- FASTQ files and splits it into smaller files based on barcode matching. It further converts each small FASTQ into FASTA file format.
 	
	Step 2: Adapter Cleaning
			 Creates adapter file for each specific ID and trim sequences based on adapter matching.
	
	Step 3: Clean PhiX and Rev Vector reads
			Identify reads associated with PhiX and RevVector and removes them from FASTA file. (Output filename: ID.noVector.fasta)
 	
	Step 4: Velvet Assembly
			Assemble reads using velvet assembler

 	Step 5: Phrap Assembly

 	Step 6: Annotation


====================================
3. Demo: 
====================================

	See "examples" folder for a demo run. The folder "01_testrun" within examples directory contains sequence files, configuration file, scripts to run parfums as well as  output directory. 

