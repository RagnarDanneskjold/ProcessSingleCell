import sys
import pysam
import Bam, GTFparse, InputOutput, Reads, Tree

THREADING = False

files = InputOutput.getInput(sys.argv, THREADING)

if not files:
    print("usage: python3 FastCount.py <path/to/GTF> <path/to/outfile> <path/to/bamfile_1> <path/to/bamfile_2> <...>",file=sys.stderr)
    sys.exit(1)

#print("Start GTF file", file=sys.stderr)
genes = GTFparse.parseGTFFile(files['gtffile'], len(files['bams']))
#print("End GTF file", file=sys.stderr)

for i in range(len(files['bams'])):
#    print("Start BAM file", file=sys.stderr)
    Bam.processBam(files['bams'][i],i,genes)
    InputOutput.printOutput(files['outfile']+"_bamlen_"+str(i),genes,files['bams'])
#    print("End BAM file", file=sys.stderr)

#print("Start outputting", file=sys.stderr)
InputOutput.printOutput(files['outfile']+"_all",genes,files['bams'])
#print("End outputting", file=sys.stderr)
