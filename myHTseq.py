import sys
import pysam
import re
import operator
from collections import defaultdict

## CONSTANTS ##

LINES_TO_CHECK = 20
UNSORTED_PERC = 0.5

## FUNCTIONS ##

def parseGTFFile (filename):
        gtf_fp = open(filename, "r")
        parsedData = dict()

        for line in gtf_fp:
                if line.startswith('#'):
                        continue
        
                line_ar = line.rstrip('\n\r').split('\t')
                chrom = line_ar[0]
                feature = line_ar[2]
                start = int(line_ar[3])
                end = int(line_ar[4])
                strand = line_ar[6]
        
                id_string = line_ar[8]
                gene_id = re.search(r'gene_id \"(.+?)\";', id_string).group(1)
        
                if feature == 'gene':
                        if chrom not in parsedData:
                                parsedData[chrom] = {
                                                '+':list(),
                                                '-':list()
                                                }
                        parsedData[chrom][strand].append({
                                        'gene_id':gene_id,
                                        'start':start,
                                        'end':end,
                                        'exons':list(),
                                        'intronic':list(),
                                        'reads':list()
                                })

                elif feature == 'exon':
                        element_index = len(parsedData[chrom][strand]) - 1
                        parsedData[chrom][strand][element_index]['exons'].append({
                                        'start':start,
                                        'end':end
                                        })
                
        
        gtf_fp.close()

        # guarantee sorted GTF data
        for chrom in parsedData:
                for strand in parsedData[chrom]:
                        parsedData[chrom][strand].sort(key=operator.itemgetter('start', 'end'))
                        for i in range(len(parsedData[chrom][strand])):
#                                print(chrom,strand,parsedData[chrom][strand][i]['start'], parsedData[chrom][strand][i]['end'])
                                parsedData[chrom][strand][i]['exons'].sort(key=operator.itemgetter('start', 'end'))

        return parsedData

def overlap(block1, block2):
        if ( block1['start'] < block2['end'] and 
                        block1['end'] > block2['start']
                        ):
                return True

        return False

def overlapBases(block1, block2):
        if not overlap(block1, block2):
                return 0

        if block1['start'] < block2['start']:
                if block1['end'] < block2['end']:
                        return block1['end'] - block2['start'] - 1 
                else:
                        return block2['end'] - block2['start'] - 1
        else:
                if block1['end'] < block2['end']:
                        return block1['end'] - block1['start'] - 1
                else:
                        return block2['end'] - block1['start'] - 1

def rangeBsearch(queryStart, queryEnd, data):
        if len(data) == 0:
#            print("ERROR: len(data) == 0!")
            return []
#        print("data[0]Start:", data[0]['start'], "data[0]End", data[0]['end'])
#        print("start", queryStart, "end", queryEnd)
        upper_bound = rangeBsearchUpper(queryStart, queryEnd, data)
        lower_bound = rangeBsearchLower(queryStart, queryEnd, data)

#        print("data:", len(data), "lower:", lower_bound, "upper:", upper_bound)

        if (lower_bound >= len(data) or
                        not overlap(data[lower_bound], {
                                'start':queryStart,
                                'end':queryEnd
                                })
                        ):
                return []

        return [lower_bound, upper_bound]


def rangeBsearchUpper(queryStart, queryEnd, data):
        upper = len(data) - 1
        lower = 0

        while lower <= upper:
                mid = (lower + upper) // 2

                if data[mid]['start'] > queryEnd:
                        upper = mid - 1
                else:
                        lower = mid + 1
        
        return lower

def rangeBsearchLower(queryStart, queryEnd, data):
        upper = len(data) - 1
        lower = 0

        while lower <= upper:
                mid = (lower + upper) // 2

                if data[mid]['end'] < queryStart:
                        lower = mid + 1
                else:
                        upper = mid - 1
        
        return lower

def checkBam(filename):
        test_fp = pysam.AlignmentFile(filename, "rb")

        pair_count = 0
        prev_name = ""
        unsorted_limit = int(UNSORTED_PERC * LINES_TO_CHECK)
        unsort_count = 0

        for algn in test_fp.head(LINES_TO_CHECK):
                curr_name = algn.query_name

                if curr_name == prev_name:
                        pair_count += 1
                elif prev_name != "" and pair_count != 2:
                        unsort_count += 1
                
                if unsort_count >= unsorted_limit:
                        test_fp.close()
                        return False

                prev_name = curr_name

        test_fp.close()

        return True

def parseFragment(readPair):
        fragment = {
                        'strand':'',
                        'start':0,
                        'end':0
                        }


#        print("name:", readPair[0].query_name)
#        print("readPair[0]:",readPair[0].reference_start, readPair[0].reference_end)
#        print("readPair[1]:",readPair[1].reference_start, readPair[1].reference_end)
 
        if readPair[0].is_reverse == readPair[1].is_reverse:
                print("ERROR READS ARE SAME DIRECTION",file=sys.stderr)
                return {}

        if readPair[0].is_read1:
                if not readPair[0].is_reverse:
                        fragment['start'] = readPair[0].reference_start
                        fragment['end'] = readPair[1].reference_end
                        fragment['strand'] = '+'
                else:
                        fragment['start'] = readPair[1].reference_start
                        fragment['end'] = readPair[0].reference_end
                        fragment['strand'] = '-'
        elif readPair[1].is_read1:
                if not readPair[1].is_reverse:
                        fragment['start'] = readPair[1].reference_start
                        fragment['end'] = readPair[0].reference_end
                        fragment['strand'] = '+'
                else:
                        fragment['start'] = readPair[0].reference_start
                        fragment['end'] = readPair[1].reference_end
                        fragment['strand'] = '-'
        
#       print(readPair[0].query_name, fragment)
        return fragment

def readQualityCheck(read):
        if ( read.is_secondary or 
                read.is_duplicate or 
                read.is_qcfail or 
                read.is_unmapped or 
                read.mate_is_unmapped or  
                not read.is_paired or 
                not read.is_proper_pair ):
                return False

        return True

def overlapGenes(gene_elements, prev_reads, chrom):
        if len(prev_reads) == 2 and chrom in gene_elements:
                # get start, end, strand of the FRAGMENT
                frag = parseFragment(prev_reads)
#                print("DEBUG:IN IF, chr:", chrom, "strand:", frag['strand'])
                if frag == {} or frag['end'] - frag['start'] < 50 or frag['end'] - frag['start'] > 600:
       #         if frag == {}:
                        return 

                overlap_indicies = rangeBsearch(frag['start'], frag['end'],

                                gene_elements[chrom][frag['strand']])

#                print("overlap_indicies = ", overlap_indicies)
                if overlap_indicies == []:
                        return 

                low_index = overlap_indicies[0]
                high_index = overlap_indicies[1]

#                print("DEBUG:read:",frag)
#                print("low_index:",low_index)
#                print("high_index:",high_index)

                if high_index != low_index :
                        # this read MAY overlap multiple elements
                        highest_overlap_index = -1
                        curr_highest = 0
                        found_equal = False
                        count = 0
                        high = high_index + 1

                        if high_index >= len(gene_elements[chrom][frag['strand']]):
                                high = high_index

                        for i in range(low_index, high):

                                score = overlapBases(gene_elements[chrom][frag['strand']][i],frag)
                                if score > 0:
                                    count += 1

                                if score > curr_highest:
                                        found_equal = False
                                        highest_overlap_index = i
                                        curr_highest = score
                                elif score == curr_highest:
                                        found_equal = True

                        if not found_equal and highest_overlap_index > -1:
#                        if count == 1:
                                gene_elements[chrom][frag['strand']][highest_overlap_index]['reads'].append(frag)

                else: # 1 element
                        gene_elements[chrom][frag['strand']][low_index]['reads'].append(frag)

# FIXME handle exons/introns later
        return 

## END FUNCTIONS ##

if (len(sys.argv) < 3):
        print("usage:", sys.argv[0], "<list with full paths of name-sorted bamfiles> <GENCODE gtf file> <outfile name>")
        sys.exit(1)

bamfiles_list_file = sys.argv[1]
gtffile = sys.argv[2]
outfilename = sys.argv[3]

print ("DEBUG: opening GTF file", file=sys.stderr)
gene_elements = parseGTFFile(gtffile)
print("DEBUG: done with GTF file", file=sys.stderr)

bamfiles = list()

bamlist_fp = open(bamfiles_list_file, "r")

for line in bamlist_fp:
        bamfiles.append(line.rstrip('\n\r'))

bamlist_fp.close()

for bam in bamfiles:
        if not checkBam(bam):
                print("ERROR: BAM file not usable.", file=sys.stderr)
                sys.exit(1)

        bam_fp = pysam.AlignmentFile(bam, "rb")

        prev_reads = list()
        prev_read_name = ""

        for read in bam_fp.fetch(until_eof = True):
                # quality check
                if not readQualityCheck(read):
                        continue
                
                if read.query_name == prev_read_name:
                        prev_reads.append(read)
                else:
                        if len(prev_reads) == 2:
                                #print("reference name:", bam_fp.get_reference_name(prev_reads[0].reference_id))
                                chrom = bam_fp.get_reference_name(prev_reads[0].reference_id)
                                overlapGenes(gene_elements, prev_reads, chrom)

                        prev_read_name = read.query_name
                        prev_reads = [read]
                

        if len(prev_reads) == 2:
                chrom = bam_fp.get_reference_name(prev_reads[0].reference_id)
                overlapGenes(gene_elements, prev_reads, chrom)
        bam_fp.close()

outfile = open(outfilename, 'w')
#outfile.write('chr\tstart\tend\tstrand\tgene_id\tcount\n')

#FIXME handle multiple counts for different bamfiles!
for chrom in gene_elements:
        for strand in gene_elements[chrom]:
                for gene in gene_elements[chrom][strand]:
                        outstr = chrom + '\t' + str(gene['start']) + '\t' + str(gene['end']) + '\t' + strand + '\t' + gene['gene_id'] + '\t' + str(len(gene['reads'])) + '\n'
                        outfile.write(outstr)

outfile.close()
                        
