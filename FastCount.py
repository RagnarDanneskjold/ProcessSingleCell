# Aparna Rajpurkar
# This is the central script for the FastCount program, all other .py modules
# except TimetTest.py are DEPENDENCIES.
# imports
import sys
import pysam
import Bam, GTFparse, InputOutput, Reads, Tree

# set threading constant to False--we are not threading in this version
THREADING = False

# check input files and add to a files dictionary
files = InputOutput.getInput(sys.argv, THREADING)

# if something was wrong with the input files, print usage and exit
if not files:
    print("usage: python3 FastCount.py <path/to/GTF> <path/to/outfile> <path/to/bamfile_1> <path/to/bamfile_2> <...>",file=sys.stderr)
    sys.exit(1)

# parse the GTF file into a data structure for later use
genes = GTFparse.parseGTFFile(files['gtffile'], len(files['bams']))

# iterate each input BAM file and process it
for i in range(len(files['bams'])):
    # Process BAM file
    Bam.processBam(files['bams'][i],i,genes)
    # output at each BAM file in case there is a crash / for temporary
    # checking out output
    InputOutput.printOutput(files['outfile']+"_TMP_bamlen_"+str(i),genes,files['bams'])

# print the final output
InputOutput.printOutput(files['outfile']+"_all",genes,files['bams'])
