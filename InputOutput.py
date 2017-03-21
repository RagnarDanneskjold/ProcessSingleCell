# Aparna Rajpurkar
# Input Output module for FastCount.py program
# imports
import Bam # my Bam.py module
import sys
import os.path

def getInputThread(args,inputdict):
    """If this is the Threaded version of the program, get input differently than\
            normal version"""
    # check number of arguments and validity of file paths        
    if len(args) < 4 or not os.path.isfile(args[1]) or not os.path.isfile(args[4]):
        return {}

    # set the thread number from input arguments
    inputdict['thread'] = args[3]

    # open the list of BAM files
    bamlist = open(args[4], "r")

    # append each bamfile name to the dictionary of BAM files
    for line in bamlist:
       filename = line.rstrip('\n\r') 
       inputdict['bams'].append(filename)

    # close the file
    bamlist.close()
 
def getInput(args,THREADING):
    """Default Input function for the FastCount.py program. Checks validity of\
            input, parses input into useful dictionary, and returns it"""

    # initialize the input dictionary
    inputdict = {
            'gtffile':'',
            'bams':list(),
            'outfile':''
            }

    # check if we're running the threaded version of the program
    # process input differently if so
    if THREADING:
        if getInputThread(args,inputdict) == {}:
            return {}
    elif len(args) < 3 or not os.path.isfile(args[1]) or not os.path.isfile(args[3]):
        # if file paths are not valid or argument number is incorrect, return
        # empty dictionary
        return {}
    else:
        # open list of BAM files
        bamlist = open(args[3], "r")

        # iterate over each BAM filename, and add to input dictionary 
        # list of BAM files
        for line in bamlist:
            filename = line.rstrip('\n\r') 
            inputdict['bams'].append(filename)

        # close the file
        bamlist.close()
 
    # set GTFFile and Output filename fields of the input dictionary
    inputdict['gtffile'] = args[1]
    inputdict['outfile'] = args[2]

    # check whether there were duplicates in the input using a SET.
    # if so, return empty dictionary
    if len(set(inputdict['bams'])) != len(inputdict['bams']):
        print("Duplicate bam files in input!",file=sys.stderr)
        return {}

    # check each bam file input for validity
    # Using the Bam.py module checkBam() function
    for bam in inputdict['bams']:
        if not os.path.isfile(bam) or not Bam.checkBam(bam):
            print("ERROR: bamfile[",bam,"] did not check out!", file=sys.stderr)
            return {}

    # return the input dictionary
    return inputdict


def printOutput(outfilename, gene_elements, bam_names):
    """output function that for FastCount.py that goes through every GeneTree\
            object and calls its output method\
            As well as prints the header to the output file"""
    # open the output file
    out_fp = open(outfilename, 'w')

    # initialize header string
    header = "chr\tstart\tend\tstrand\tgene_id"

    # construct header from bam file names
    for filename in bam_names:
        header += "\t" + filename

    header += "\n"
    # write to file
    out_fp.write(header)

    # iterate over every chromosome and strand, calling the output method
    # for each GeneTree
    # which handles outputing the read count for each gene
    for chrom in sorted(gene_elements):
        for strand in sorted(gene_elements[chrom]):
            gene_elements[chrom][strand].writeTree(out_fp)

    # close output file
    out_fp.close()
