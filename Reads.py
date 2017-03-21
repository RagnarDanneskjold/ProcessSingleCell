# Aparna Rajpurkar
# Reads module for FastCount.py program

def parseFragment(readPair,chrom):
    """function to take a raw BAM read pair and return a parsed fragment dictionary\
            with all the useful information"""

    # Paired end reads need reads to be pointing in different directions
    # if this is not the case, return empty dictionary
    if readPair[0].is_reverse == readPair[1].is_reverse:
        return {}
    
    # intitialize the fragment dictionary for this read.
    fragment = {
        'chrom':chrom,
        'strand':'',
        'start':0,
        'end':0
    }
 
    # initialize the reads as Read 1 or Read 2
    read1 = readPair[0]
    read2 = readPair[1]

    # switch them if the second read in the pair is actually Read 1
    if readPair[1].is_read1:
        read1 = readPair[1]
        read2 = readPair[0]

    # infer strandedness of the fragment from orientation
    # of Read 1 Read 2 pair, and set the fields of the fragment dictionary
    # accordingly
    if not read1.is_reverse:
        fragment['start'] = read1.reference_start
        fragment['end'] = read2.reference_end
        fragment['strand'] = '+'
    else:
        fragment['start'] = read2.reference_start
        fragment['end'] = read1.reference_end
        fragment['strand'] = '-'

    # return the parsed fragment
    return fragment

def readQualityCheck(read):
    """Function to check the quality of a raw read"""
    # As pysam's Aligment object contains all the quality information about a read,
    # check its fields for unacceptable values
    # if this read is unacceptable, return False, else return True
    if ( read.is_secondary or 
        read.is_qcfail or 
        read.is_unmapped or 
        read.mate_is_unmapped or  
        not read.is_paired or 
        not read.is_proper_pair 
        ):

        return False

    return True

