# Aparna Rajpurkar
# Bam module for FastCount.py

# imports #
import pysam # an Object Oriented module for reading BAM files
import Reads # my Reads.py module

## CONSTANTS ##
LINES_TO_CHECK = 20
UNSORTED_PERC = 0.5
##

def checkBam(filename):
    """function that checks BAM file for validity and name-sorting"""
    # open input BAM file with pysam
    test_fp = pysam.AlignmentFile(filename, "rb")

    # initialize count variables used to check sorting
    pair_count = 0
    prev_name = ""
    # limit the number of lines we check so speed is not too affected
    unsorted_limit = int(UNSORTED_PERC * LINES_TO_CHECK)
    unsort_count = 0

    # iterate over specified number of lines in the bam file 
    for algn in test_fp.head(LINES_TO_CHECK):
        # get name from pysam Alignment Object
        curr_name = algn.query_name

        # Check name sorting by comparing this name to previous name
        if curr_name == prev_name:
            pair_count += 1
        elif prev_name != "" and pair_count != 2:
            # if we don't have pairs of 2, add to unsorted count
            unsort_count += 1
                
        # if we have too many unsorted lines, return False
        if unsort_count >= unsorted_limit:
            test_fp.close()
            return False

        prev_name = curr_name

    # close file
    test_fp.close()

    # return true if unsorted count is within tolerance
    return True

def processBam(bamfile, bam_num, genes):
    """Function that processes each line of an input BAM file and calls Tree functions\
            to search for gene interval matches"""

    # open BAM file
    bam_fp = pysam.AlignmentFile(bamfile, "rb")
    # initialize pair counting varibles 
    prev_reads = list()
    prev_read_name = ""

    # iterate over every read in the BAM file
    for read in bam_fp.fetch(until_eof = True):
        # check quality using a function from Reads.py
        # if quality does not check out, go to next read
        if not Reads.readQualityCheck(read):
            continue

        # check pairedness
        if read.query_name == prev_read_name:
            # if this read is a pair with the previous read, add it to the list
            prev_reads.append(read)
        else:
            # if we're on a new read, process previous pair & check interval 
            # overlaps with GTF file genes
            if len(prev_reads) == 2:
                # get chromosome ID
                chrom = bam_fp.get_reference_name(prev_reads[0].reference_id)
                # parse paired-end read into a fragment dictionary
                frag = Reads.parseFragment(prev_reads,chrom)

                if frag != {} and chrom in genes:
                    # if fragment is valid and chromosome exists in GTF datastructure
                    # then look for overlap with genes
                    genes[chrom][frag['strand']].overlapInterval(frag,bam_num)

            # reset pair counters
            prev_read_name = read.query_name
            prev_reads = [read]


    # be sure to check the LAST pair in the BAM file!
    # this will be skipped in the loop
    # same process as what we do in the loop
    if len(prev_reads) == 2:
        chrom = bam_fp.get_reference_name(prev_reads[0].reference_id)
        frag = Reads.parseFragment(prev_reads,chrom)

        if frag != {} and chrom in genes:
            genes[chrom][frag['strand']].overlapInterval(frag, bam_num)

    # close file
    bam_fp.close()
