#!/usr/bin/perl -w
use strict;

my ($infile)=@ARGV;

#Edit this one if needed
my $genomefastafile="/home/CAM/jcotney/GENOME/mm9/dna/mm9_nh.fa";

open(IN,$infile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Viewpoint/;
	my @array=split/\t/;
	
	unless (-e $array[0]) {
		mkdir $array[0];
		mkdir "$array[0]/primers";
		mkdir "$array[0]/bam";
	}
	
	my @conditions=split(";",$array[3]);
	my @replicates=get_replicates(@conditions);
	my @bams=split(";",$array[4]);
	my $conditionsR=convert_list_toR(@conditions);
	my $replicatesR=convert_list_toR(@replicates);
	my $bamsR=convert_list_toR(@bams);
	
	#print parameter file
	my $parameterfile=$array[0]."/".$array[0]."_parameters.txt";
	open(OUT,">$parameterfile") || die $!;
	print OUT <<EOF;

#parameters
#necessary files/folders
projectname= "$array[0]"
primerFile = "primers/$array[0]_primers.fas"
firstprimer = "$array[0]_forward"
referenceGenomeFile = "$genomefastafile"
bamFilePath = "bam"
bamfiles=$bamsR

conditions=$conditionsR
replicates=$replicatesR

reSequence1="$array[5]"
reSequence2="$array[6]"

resultfile="$array[7]"
scatterplotfile=sub(".txt","_scatterplot.pdf",resultfile)
fitplotfile=sub(".txt","_fit.pdf",resultfile)
peakplotfile=sub(".txt","_peaks.pdf",resultfile)

EOF
	close OUT;
	
	#print primer file
	open(OUT,">$array[0]/primers/$array[0]_primers.fas") || die $!;
	print OUT ">$array[0]_forward\n";
	print OUT $array[1],"\n";
	print OUT ">$array[0]_reverse\n";
	print OUT $array[2],"\n";
	close OUT;
}
close IN;

sub convert_list_toR {
	my @input=@_;
	return "c(\"".join("\",\"",@input)."\")";
}

sub get_replicates {
	my @input=@_;
	my %cons;
	my @reps;
	foreach my $c (@input) {
		$cons{$c}++;
		push @reps, $cons{$c};
	}
	return @reps;
}
