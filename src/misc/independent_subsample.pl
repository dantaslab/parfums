#!/usr/bin/perl -w
	use strict;

	
	&check_point()

sub check_point {
	my $dirpath = $_[0];
	my $suffix = $_[1];
	foreach my $seqID (keys %inputSeqFile) {
		if (exists $splitFiles{$seqID}) {
			my @allFiles = @{$splitFiles{$seqID}};
			foreach my $filename (@allFiles) {
				if (-z "$filename.$suffix") {
					
			}
		}
	}
	
