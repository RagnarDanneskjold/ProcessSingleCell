import sys
import pysam
import Bam, GTFparse, InputOutput, Reads, Tree

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

for read in bam_fp.fetch(until_eof = True):
    if not Reads.readQualityCheck(read):
        continue

    if read.query_name == prev_read_name:
        prev_reads.append(read)
    else:
        if len(prev_reads) == 2:
            chrom = bam_fp.get_reference_name(prev_reads[0].reference_id)
            frag = Reads.parseFragment(prev_reads,chrom)

            if frag != {} and chrom in genes:
                genes[chrom][frag['strand']].overlapInterval(frag)

        prev_read_name = read.query_name
        prev_reads = [read]


if len(prev_reads) == 2:
    chrom = bam_fp.get_reference_name(prev_reads[0].reference_id)
    frag = Reads.parseFragment(prev_reads,chrom)

    if frag != {} and chrom in genes:
        genes[chrom][frag['strand']].overlapInterval(frag)

bam_fp.close()
print("End BAM file", file=sys.stderr)

print("Start outputting", file=sys.stderr)
InputOutput.printOutput(files['outfile'],genes)
print("End outputting", file=sys.stderr)
