#!/usr/bin/perl -w

# Written by Alejandro Reyes
# Receives a Fasta File and a list of names and retrieve the sequences that begin with the name
# or the sequences that dont begin with that name.
#
# Usage: GetSequences.pl "Names list" "Fasta File" > "output"
#
# Input:
# Output: 
# Note:
# Created: Aug 11 / 08.
# Last-updated:

use strict;

if (@ARGV != 2) {
	die "\nUsage: GetSequences.pl <Names list> <Fasta File> > <output>\n\n";
}


#Reads the names to get

my $file = shift @ARGV;
open (IN, "<$file") or die ("Couldn't open file: $file\n");

my %names = ();

while (my $lines=<IN>) {
	if ($lines =~ /\S+/){
	chomp $lines;
        #my @entry = split (/\-/, $lines);
	my @entry = split (/\s+/, $lines);
	$names{$entry[0]}=1;
	}	
}	
close IN;


#Reads the Fasta Sequences

$file = shift @ARGV;
open (IN, "<$file") or die ("Couldn't open file: $file\n");

my $counter = 0;
while (my $line = <IN>) {
	if ($line =~ /^>/){
		$counter = 0;
		chomp $line;
		$line =~ s/>//;			#Only if the name does not contain >, otherwise comment this line
		my @temp = split (/\s+/, $line);
#		my @temp = split (/\-/, $line);

		if (exists $names{$temp[0]}){
		#if (!(exists $names{$temp[0]})){ 	#Modification for don't exist, if the name is not in the list, print the sequence
			print ">$line\n";
			$counter = 1;
		}
	}
	
	elsif (($line !~ /^\@/) && ($counter == 1)){
		print $line;
	}
	
}
close IN;
exit;
