#!/usr/bin/perl
## written by Toshiya Ando @ NIBB, Japan
## 2019.03.04

use strict;
use List::Util qw(max);
use List::Util qw(min);
use Getopt::Long qw(:config auto_help);
use Pod::Usage;
use File::Basename 'basename', 'dirname';
use Try::Tiny;

# VERSION

=head1 SYNOPSIS

	filter_repeat.pl <-i file> [options]

	Example 1: Default settings (Batch)
		filter_repeat.pl -i [directory including input files (*.easyfig.out)]

	Example 2: Default settings (Single)
		filter_repeat.pl -s -i [an input file]

	Example 3: Custom settings (Single)
		filter_repeat.pl -s -i [an input file] -c [cutoff threshold, integer] -m [minimal repetitive hit size to be desplayed, integer]

=head1 DESCRIPTION

	Filter multiple hits between the two DNA sequences.

=cut

# OPTIONS

GetOptions(	'i=s' => \my $input,
			'c:i' => \my $cutoff,
			'm:i' => \my $min,
			'f:s' => \my $file,
			's' => \my $single,
			 '<>' => \&unknownArgHandler,
			 'help' => \my $help) or pod2usage(-verbose => 1) or pod2usage(q(-verbose) => 1);
pod2usage(-verbose => 2) if $help;
pod2usage(-verbose => 1) unless (defined($input));

#	set the smallest size of BLAST hits with multiple hits
#	Larger repeats will be displayed in the result file.

$cutoff = 1250 if($cutoff == 0);

#	set the threshold to filter too short BLAST hits
#	<250 may be better for interspecific comparison

$min = 300 if($min == 0);


sub unknownArgHandler{
	print STDERR "undefined Option! @_\n";
	exit;
}


# MAIN8

# TRY
try{

	my @file = ();

	if($single){
		my $dir = dirname $input or die;
		my $fname = basename $input or die;
		chdir("$dir") or die("$dir: $!");
		die ("$input not found : $!") unless(-f $fname);
		push @file, $input;
	}else{
		chdir("$input") or die("$input: $!");
		@file = glob("*.easyfig.out") or die ("$input/*.easyfig.out not found:$!");
	}

# Eliminate too short hits
# alternative:	system 'awk -F\'\t\' \'($4 > '.$min.'){print $0}\' '.$file.' > temp0.txt';

	foreach my $file (@file){
		open IN0, "$file" or die ("$file: $!");
		open OUT0, ">> temp0.txt";

		while(<IN0>){
			my @F = split(/\t/, $_);
			print OUT0 $_ if($F[3] >= $min);
		}
		close IN0; close OUT0;

# Sort hits by position information of the query seuqnece
# alternative:	system 'sort -n -k 7 temp0.txt > temp1.txt';

		open IN1, "temp0.txt";
		open OUT1, ">> temp1.txt";

		&sort_table(\*IN1, \*OUT1, 6);

		close IN1; close OUT1;

# Eliminate neighboring, overlapping hits against different regions in the target sequence

		open IN2, "temp1.txt";
		open OUT2, ">> temp2.txt";

		my $preend = "";
		my $prelen = "";
		my $prst = "";
		my $flag = 0;
		my $count = 0;
		my $pflag = 0;

		while(<IN2>){
		chomp $_;
		my @temp = split /\t/, $_;

		if($count == 0){
			$preend = max($temp[6],$temp[7]);
			$prelen = max($temp[6],$temp[7]) - min($temp[6],$temp[7]);
			$prst = $_;
			$count ++;
			next;
		}

		if(($temp[6] + $temp[7])/2 > $preend){
			if($flag == 0){
			print OUT2 $prst."\n";
			$pflag = 1;
			}elsif($flag == 1){
			$flag = 0;
			if($prelen > $cutoff){
				print OUT2 $prst."\n";
				$pflag = 1;
			}
			}else{
			print "error";
			}
		}else{
			$flag = 1;
			if($prelen > $cutoff){
			print OUT2 $prst."\n";
			$pflag = 1;
			}
		}

		$preend = max($temp[6],$temp[7]);
		$prst = $_;
		$prelen = max($temp[6],$temp[7]) - min($temp[6],$temp[7]);
	#		print $preend."\t".$temp[6]."\t".$temp[7]."\t".$flag."\t".$prelen."\t".$pflag."\t\t".$_."\n";
		$pflag = 0;
		}

		if($flag == 0){
		print OUT2 $prst."\n";
		}
		if($prelen > $cutoff){
		print OUT2 $prst."\n";
		$pflag = 1;
		}

		close IN2;
		close OUT2;

# Sort hits by position information of the target seuqnece
# alternative:	system 'sort -n -k 9 temp2.txt > temp3.txt';
#=pod

		open IN3, "temp2.txt";
		open OUT3, ">> temp3.txt";

		&sort_table(\*IN3, \*OUT3, 8);

		close IN3; close OUT3;
#=cut


# Eliminate neighboring, overlapping hits against different regions in the query sequence

		open IN4, "temp3.txt";
		open OUT4, ">> temp4.txt";


		my $preend = "";
		my $prelen = "";
		my $prst = "";
		my $flag = 0;
		my $count = 0;

		while(<IN4>){
		chomp $_;
		my @temp = split /\t/, $_;

		if($count == 0){
			$preend = max($temp[8],$temp[9]);
				$prelen = max($temp[8],$temp[9]) - min($temp[8],$temp[9]);
			$prst = $_;
			$count ++;
			next;
		}

		if(($temp[8] + $temp[9])/2 > $preend){
			if($flag == 0){
			print OUT4 $prst."\n";
			$pflag = 1;
			}elsif($flag == 1){
			$flag = 0;
					if($prelen > $cutoff){
						print OUT4 $prst."\n";
				$pflag = 1;
			}
			}else{
			print "error";
			}
		}else{
			$flag = 1;
			if($prelen > $cutoff){
			print OUT4 $prst."\n";
			$pflag = 1;
			}
		}

	#		print $preend."\t".$temp[8]."\t".$temp[9]."\t".$flag."\t".$prelen."\t".$pflag."\t\t".$_."\n";

		$preend = max($temp[8],$temp[9]);
		$prelen = max($temp[8],$temp[9]) - min($temp[8],$temp[9]);
		$prst = $_;
		$pflag = 0;
		}

		if($flag == 0){
		print OUT4 $prst."\n";
		}
		if($prelen > $cutoff){
		print OUT4 $prst."\n";
		$pflag = 1;
		}


		close IN4; close OUT4;


# Sort hits by position information of the query seuqnece again
# alternative:	system 'sort -n -k 7 temp4.txt > '.$file.'.filtered.txt';


		unlink $file.".filtered.txt" if(-f $file.".filtered.txt");

		open IN5, "temp4.txt";
		open OUT5, ">> ".$file.".filtered.txt";
		&sort_table(\*IN5, \*OUT5, 6);

		close IN5; close OUT5;

# remove temporal files
		unlink <temp[0-4].txt>;
	}


# Report
	print "#=== Filter Duplicates ===\n";
	print "# source:\n";
	foreach my $file (@file){
		print "#\t$file\n";
	}
	print "# filtered files:\n";
	foreach my $file (@file){
	  print "#\t".$file.".filtered.txt\n";
	}

	print "# script:\n#\tfilter_repeat_2.pl\n";
	print "# date:\n#\t".localtime."\n";
	print "# author:\n#\tToshiya Ando <toshiya.ando.ca\@gmail.com>\n";
	print "#\n";
	#end	
# CATCH
} catch {
	warn "caught error: $_";
};


sub sort_table{
	local *IN = $_[0];
	local *OUT = $_[1];
	my $col = $_[2];
	my @tmp = ();
	while(<IN>){
			my @F = split(/\t/, $_);
			my $hash = {
						"pos" => $F[$col],
						"row" => $_
			};
			push(@tmp, $hash);
		}
		my @tmp2 = sort {$a->{"pos"} <=> $b->{"pos"}} @tmp;
		
		foreach my $tmp (@tmp2){
				print OUT $tmp->{'row'};
	}
}
