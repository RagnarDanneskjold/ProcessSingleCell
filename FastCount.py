import sys
import pysam
import Bam, FindGenes, GTFparse, InputOutput, Overlaps, Reads, Searches

files = InputOutput.getInput(sys.argv)

if not files:
    print("usage: python3 FastCount.py <path/to/bamfile> <path/to/GTF> <path/to/outfile>")
    sys.exit(1)

genes = GTFparse.parseGTFFile(files['gtffile'])

bam_fp = pysam.AlignmentFile(files['bamfile'], "rb")
prev_reads = list()
prev_read_name = ""

for read in bam_fp.fetch(until_eof = True):
    if not Reads.readQualityCheck(read):
        continue

    if read.query_name == prev_read_name:
        prev_reads.append(read)
    else:
        FindGenes.runOverlapGenes(prev_reads, bam_fp, genes)

        prev_read_name = read.query_name
        prev_reads = [read]

FindGenes.runOverlapGenes(prev_reads, bam_fp, genes)
bam_fp.close()

InputOutput.printOutput(files['outfile'],genes)
