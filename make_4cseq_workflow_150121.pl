#!/usr/bin/perl -w
use strict;

my ($infile, $folder, $sh1, $sh2) = @ARGV;

#edit this one if needed
my $fourceseq="/home/CAM/jcotney/TOOLS/FourCSeq/run_fourcseq_workflow_150316.R";

my $folder=$folder;

#move bam files
#generate simplequeque script
#my $folder="/home/jy344/Scratch/Projects/4CSeq/RCData3";
#open(IN,"4cSeq_input_edi.txt") || die $!;
#open(OUT,">rcdata3_fourcseq.sh") || die $!;
#open(OUT2,">rcdata3_mvbams.sh") || die $!;

open(IN,$infile) || die $!;
open(OUT,">$sh1") || die $!;
open(OUT2,">$sh2") || die $!;

while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Viewpoint/;
	my @array=split/\t/;
	print OUT "source ~/.bashrc;cd $folder/$array[0];Rscript $fourceseq $array[0]_parameters.txt\n";
	print OUT2 "mv bam/*$array[0]*.bam $array[0]/bam/\n";
}
close OUT;
close OUT2;
