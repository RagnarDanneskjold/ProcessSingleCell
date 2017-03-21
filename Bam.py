import pysam
import Reads

## CONSTANTS ##
LINES_TO_CHECK = 20
UNSORTED_PERC = 0.5
##

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

def processBam(bamfile, bam_num, genes):
    bam_fp = pysam.AlignmentFile(bamfile, "rb")
    prev_reads = list()
    prev_read_name = ""

    for read in bam_fp.fetch(until_eof = True):
#        print("Bam:",bam_num,"Read:",read.query_name)
        if not Reads.readQualityCheck(read):
            continue

        if read.query_name == prev_read_name:
            prev_reads.append(read)
        else:
            if len(prev_reads) == 2:
                chrom = bam_fp.get_reference_name(prev_reads[0].reference_id)
                frag = Reads.parseFragment(prev_reads,chrom)

                if frag != {} and chrom in genes:
                    genes[chrom][frag['strand']].overlapInterval(frag,bam_num)

            prev_read_name = read.query_name
            prev_reads = [read]


    if len(prev_reads) == 2:
        chrom = bam_fp.get_reference_name(prev_reads[0].reference_id)
        frag = Reads.parseFragment(prev_reads,chrom)

        if frag != {} and chrom in genes:
            genes[chrom][frag['strand']].overlapInterval(frag, bam_num)

    bam_fp.close()
