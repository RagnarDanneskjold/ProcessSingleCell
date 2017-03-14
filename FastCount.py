import sys
import pysam
import Bam, FindGenes, GTFparse, InputOutput, Overlaps, Reads, Searches

files = InputOutput.getInput(sys.argv)

if not files:
    print("usage: python3 FastCount.py <path/to/bamfile> <path/to/GTF> <path/to/outfile>")
    sys.exit(1)

print("Start GTF file", file=sys.stderr)
gtf_fp = open(files['gtffile'],"r")
genes = GTFparse.parseGTFFile(gtf_fp)
gtf_fp.close()
print("End GTF file", file=sys.stderr)

print("Start BAM file", file=sys.stderr)
bam_fp = pysam.AlignmentFile(files['bamfile'], "rb")
prev_reads = list()
prev_read_name = ""

count = 0
for read in bam_fp.fetch(until_eof = True):
    if not Reads.readQualityCheck(read):
        continue

    if read.query_name == prev_read_name:
        prev_reads.append(read)
    else:
        count+=1
#        print("FindGenes read",count)
        FindGenes.runOverlapGenes(prev_reads, bam_fp, genes)

        prev_read_name = read.query_name
        prev_reads = [read]

FindGenes.runOverlapGenes(prev_reads, bam_fp, genes)
bam_fp.close()
print("End BAM file", file=sys.stderr)

print("Start outputting", file=sys.stderr)
InputOutput.printOutput(files['outfile'],genes)
print("End outputting", file=sys.stderr)
