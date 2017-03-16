import sys
import pysam
import Bam, FindGenes, GTFparse, InputOutput, Reads, Tree

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
badreads = dict()

for read in bam_fp.fetch(until_eof = True):

    print("on:", read.query_name)
#    if read.query_name == "HWI-ST999:184:C44V8ACXX:7:1101:1362:54252":
#        print("read of interest")

    if not Reads.readQualityCheck(read, badreads):
        print("skipping:", read.query_name)
        prev_read_name = ""
        continue

    if read.query_name == prev_read_name:
        print("appending:", read.query_name)
        print ("read query:", read.query_name, "prev_read", prev_read_name)
        prev_reads.append(read)
    else:
        #print("on find_genes:" prev_read_name)
        FindGenes.runOverlapGenes(prev_reads, bam_fp, genes)

        print("appending:",read.query_name)
        prev_read_name = read.query_name
        prev_reads = [read]

FindGenes.runOverlapGenes(prev_reads, bam_fp, genes)
bam_fp.close()
print("End BAM file", file=sys.stderr)

print("Start outputting", file=sys.stderr)
InputOutput.printOutput(files['outfile'],genes)
print("End outputting", file=sys.stderr)
