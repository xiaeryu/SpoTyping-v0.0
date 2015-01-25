#!usr/bin/perl

## Copyright (C) 2014 Xia Eryu (xiaeryu@nus.edu.sg).
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, see 
## http://www.opensource.org/licenses/gpl-3.0.html

## SpoTyping.pl
## --------------------------------
## Please report bugs to:
## xiaeryu@nus.edu.sg

=head1 NAME

 SpoTyping - in silico prediction of Mycobacterium tuberculosis spoligotypes

=head1 SYNOPSIS

 perl SpoTyping.pl [OPTIONS]... FASTQ_1 [FASTQ_2]

 -h             displays help message and exit
 --swift        swift mode, either [on] or [off]
                Defulat: on
 --min          minimum number of exact reads to support existence of loci
                Default: 5
 --rmin         minimum number of inexact reads, allowing for 1 mismatch to support existence of loci
                Default: 6
 --output       basename of all output files generated
                Default: SpoTyping
 --blast        specified path to bin directory containing blastn and makeblastdb executables
                Default: takes the system specified blastn and makeblastdb if available
 -d             enable debug mode, keeping all intermediate files for checking
                Default is off
 FASTQ_1        input FASTQ read1 file(mandatory)
 FASTQ_2        input FASTQ read 2 file (optional for paired end reads)

=head1 AUTHOR

 Xia Eryu, xiaeryu@nus.edu.sg

=head1 VERSION

 1.00

=cut


#use strict;
#use warnings;
use Getopt::Long;
use File::Basename qw<basename dirname>;
use Pod::Usage;

## Global variables.
my $dir=dirname $0; 		# script directory
my $setlength=50*5000000;	# base input cut-off for swift mode

## Option variables.
my $opt_help;                   # help message
my $opt_debug;                  # debug mode, default is off
my $swift='on';             	# swift mode, default is on
my $input1;			# input fastq file 1
my $input2;			# input fastq file 2
my $min=5;			# minimum number of exact reads to support existence of spacer sequence
my $min_relax=6;		# minimum number of approximate reads, allowing for 1 mismatch to support existence of spacer sequence
my $output="SpoTyping";		# basename of output files
my $blast;			# basename of bin directory containing the NCBI BLAST+ executables

## Getting options.
Getopt::Long::Configure('bundling');
if(
!GetOptions(
	"h"       => \$opt_help,
	"d"       => \$opt_debug,
	"swift:s" => \$swift,
        "output=s"=> \$output,
        "min=i"   => \$min,
        "rmin=i"  => \$min_relax,
	"blast=s" => \$blast,
) || $opt_help
){
	pod2usage(-verbose=>2) if ($opt_help);
	pod2usage(-verbose=>0 );
}

## Check options.
my $num_files=scalar(@ARGV);
if($num_files<1||$num_files>2){
	pod2usage(q/invalid num of input fastq files specified/);
}
unless(defined $ARGV[0] && -e $ARGV[0]){
	pod2usage(q/invalid FASTQ_1 file specified/);
}

$input1=$ARGV[0];

if($num_files==2){
	unless(defined $ARGV[1] && -e $ARGV[1]){
		pod2usage(q/invalid FASTQ_2 file specified/);
	}    
	$input2=$ARGV[1];
}

if(defined $blast){
    unless(-e "$blast/makeblastdb" && -e "$blast/blastn"){
        pod2usage(q/user specified blast directory does not contain blastn and makeblastdb executables/);
    }
}else{
    my $path_to_blastn=`which blastn`;
    my $path_to_makeblastdb=`which makeblastdb`;
    unless ($path_to_blastn && $path_to_makeblastdb){
        pod2usage(q/unable to find system blastn and|or makeblastdb executables, please specify blast directory using script parameters/);
    }
}

## Checking name of existing tmp files.
my $tmpfile=0;
while(-e "$output.SpoTyping.tmp.$tmpfile"){
	$tmpfile++;
}

##########################################################
## Create a fasta file with the reads concatenated.
##########################################################
my $count=0;
open(TMP,"> $output.SpoTyping.tmp.$tmpfile") or die "Cannot write tmp file to the current directory:$!\n";
print TMP "> Combine | $input1 $input2\n";

if($swift eq "on"){
	my $outlength;
	open(INPUT,$input1) or die "Cannot open input file:$!\n";
	while(<INPUT>){
        	chomp;
		if($outlength>$setlength){
			last;
		}elsif($count%4==1){
			s/\s+//g;
        	        print TMP "$_";
			$outlength+=length($_);
	        }
        	$count=($count+1)%4;
	}
	close INPUT;
	$count=0;

	if(defined($input2) && $outlength<$setlength){
        	open(INPUT,$input2) or die "Cannot open input file:$!\n";
	        while(<INPUT>){
        	        chomp;
			if($outlength>$setlength){
                		last;
	        	}elsif($count%4==1){
				s/\s+//g;
        	                print TMP "$_";
				$outlength+=length($_);
                	}
	        	$count=($count+1)%4;
		}
        	close INPUT;
	}
}else{
	open(INPUT,$input1) or die "Cannot open input file:$!\n";
        while(<INPUT>){
                chomp;
                if($count%4==1){
			s/\s+//g;
                	print TMP "$_";
                }
                $count=($count+1)%4;
        }
        close INPUT;
	$count=0;

        if(defined($input2)){
                open(INPUT,$input2) or die "Cannot open input file:$!\n";
                while(<INPUT>){
                        chomp;
                        if($count%4==1){
				s/\s+//g;
                                print TMP "$_";
                        }
                	$count=($count+1)%4;
		}
                close INPUT;
        }
}
print TMP "\n";
close TMP;

##########################################################
### Blast the spacers against the concatenated fasta file.
##########################################################
if($blast=~/\w+/){
	system("$blast/makeblastdb -in $output.SpoTyping.tmp.$tmpfile -out $output.SpoTyping.tmp.$tmpfile -dbtype nucl");
	system("$blast/blastn -query $dir/ref/spacer.fasta -db $output.SpoTyping.tmp.$tmpfile -task blastn -dust no -outfmt 7 -max_target_seqs 1000000 > $output.SpoTyping.tmp.$tmpfile.blast.out");
}else{
	system("makeblastdb -in $output.SpoTyping.tmp.$tmpfile -out $output.SpoTyping.tmp.$tmpfile -dbtype nucl");
	system("blastn -query $dir/ref/spacer.fasta -db $output.SpoTyping.tmp.$tmpfile -task blastn -dust no -outfmt 7 -max_target_seqs 1000000 > $output.SpoTyping.tmp.$tmpfile.blast.out");
}

my %record;
my %record_relax;

##########################################################
### Parsing blast output.
##########################################################
open(BLAST,"$output.SpoTyping.tmp.$tmpfile.blast.out") or die "Cannot open input file SpoTyping.tmp.$tmpfile.blast.out:$!\n";
while(<BLAST>){
        chomp;
        if(!/^#/){
                my @tmp=split(/\s+/);
                if($tmp[2]==100 && $tmp[3]==25){
                        $record{$tmp[0]}++;
                        $record_relax{$tmp[0]}++;
                }elsif(($tmp[2]==96 && $tmp[3]==25) || ($tmp[2]==100 && $tmp[3]==24)){
                        $record_relax{$tmp[0]}++;
                }
        }
}
close BLAST;

##########################################################
### Writing to the output file and the log file.
##########################################################
my $logname=$output.".log";
open(LOG, ">> $logname") or die "Cannot write log file to the output directory:$!\n";
my @storage;
print LOG "## $input1 $input2\n";
print LOG "Spacer\tError-free_number\t1-error-tolerant_number\tCode\n";
for(my $i=1;$i<=43;$i++){
        my $signal;
        if($record{"Spacer".$i}>=$min || $record_relax{"Spacer".$i}>=$min_relax){
                $signal=1;
        }
        push(@storage,$signal);
        printf LOG ("Spacer%d\t%d\t%d\t%d\n",$i,$record{"Spacer".$i},$record_relax{"Spacer".$i},$signal);
}

open(OUT,">> $output") or die "Cannot write to the output file:$!\n";
if(defined($input2)){
	print OUT "$input1&$input2\t";
}else{
	print OUT "$input1\t";
}
for(my $i=0;$i<=39;$i=$i+3){
        my $print=4*$storage[$i]+2*$storage[$i+1]+$storage[$i+2];
        print OUT "$print";
}
printf OUT ("%d",$storage[42]);
print OUT "\n";

close OUT;
close LOG;

unlink glob("$output.SpoTyping.tmp.$tmpfile*") unless (defined $opt_debug);
