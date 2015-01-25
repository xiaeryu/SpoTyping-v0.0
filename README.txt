SpoTyping is a software for predicting spoligotype from sequencing reads.

>>> Prerequisites:
1. Perl
2. Perl modules: Getopt::Long, File::Basename, Pod::Usage;
3. BLAST

>>> Input:
Fastq file or pair-end fastq files.

>>> Output:
In the output file specified:	predicted spoligotype in the format of octal code.
In the output log file:		count of hits from BLAST result for each spacer sequence. 

>>> Usage:
perl SpoTyping.pl [OPTIONS]... FASTQ_1 [FASTQ_2]

An Example call:
perl SpoTyping.pl read_1.fastq read_2.fastq –output spo.out 

>>> Options:
-h             displays help message and exit
--swift        swift mode, either [on] or [off]
               [Defulat] on
--min          minimum number of exact reads to support existence of loci
               [Default] 5
--rmin         minimum number of inexact reads, allowing for 1 mismatch to support existence of loci
               [Default] 6
--output       basename of all output files generated
               [Default] SpoTyping
--blast        specified path to bin directory containing blastn and makeblastdb executables
               [Default]takes the system specified blastn and makeblastdb if available
-d             enable debug mode, keeping all intermediate files for checking
               [Default] off

FASTQ_1        input FASTQ read1 file(mandatory)
FASTQ_2        input FASTQ read 2 file (optional for paired end reads)

>>> Suggestions:
It's highly suggested to use the swift mode (set as the default).
If you do wish to take in all reads, it's suggested to estimated the coverage first. 
The --min is suggested to be set to 1/10 of the estimation to increase accuracy and to eliminate false positive. 
At the same time, --rmin should also be adjusted to be a bit larger than min.
