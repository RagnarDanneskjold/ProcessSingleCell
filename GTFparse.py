# Aparna Rajpurkar
# This is the GTFparse module of the FastCount.py program
# imports
import re
import operator
from Tree import GeneTree # my Object Oriented classes

def parseValidLine(line):
    """Function that parses a line of a GTF file and returns a useful \
            data structure of its fields"""
    # split line into array
    line_ar = line.rstrip('\n\r').split('\t')

    # if line is too short, return nothing
    if len(line_ar) < 9:
        return {}

    # Grab the gene ID from the line using regular expressions
    id_string = line_ar[8]
    gene_id = re.search(r'gene_id \"(.+?)\";', id_string).group(1)

    # construct the results dictionary for this line
    result = {
            'chrom': line_ar[0],
            'feature': line_ar[2],
            'start': int(line_ar[3]),
            'end': int(line_ar[4]),
            'strand': line_ar[6],
            'gene_id': gene_id
        }

    # We are only interested in gene and exon features, so return 
    # nothing if not gene or exon
    if result['feature'] != "gene" and result['feature'] != "exon":
        return {}

    # return the result dictionary
    return result

def parseGTFFile (gtffile, bam_num):
    """Function that handles parsing the GTF file and intializing the GeneTree\
            Objects for each chromosome and strand"""

    # open the GTF file and initialize local variables
    gtf_fp = open(gtffile,"r")
    parsedData = dict()
    curr_gene = 0

    # iterate over every line in the GTF file
    for line in gtf_fp:
        # skip if this is a header line
        if line.startswith('#'):
            continue

        # parse line into fields dictionary
        fields = parseValidLine(line)

        # skip if we could not parse, or feature is not of interest
        if not fields:
            continue
        
        # if we're on a new chromosome, initialize its GeneTree objects
        if fields['chrom'] not in parsedData:
            # set this chromosome's strand dictionary
            parsedData[fields['chrom']] = dict()
            # for each strand, intitialize a GeneTree object
            # which will store all entries for its genes
            parsedData[fields['chrom']]['+'] = GeneTree(fields['chrom'],'+')
            parsedData[fields['chrom']]['-'] = GeneTree(fields['chrom'],'-')

        # if this feature is a gene, add it to the GeneTree
        if fields['feature'] == 'gene':
            # call the addNode method of the GeneTree object on this gene
            curr_gene = parsedData[fields['chrom']][fields['strand']].addNode(fields,fields['gene_id'], bam_num)
        else: # exon
            # if this is an exon, add it to the current gene's Tree
            parsedData[fields['chrom']][fields['strand']].addExon(fields,curr_gene)

    # close the GTF file
    gtf_fp.close()

    # for each chromosome and strand, call the GeneTree object's balance method
    # to ensure optimally efficient find() operations later
    for chrom in parsedData:
        for strand in parsedData[chrom]:
            parsedData[chrom][strand].balanceAll()

    # return our data structure
    return parsedData

