import sys
import pysam
import re
from collections import defaultdict

if (len(sys.argv) < 2):
	print("usage:", sys.argv[0], "<bamfile> <gtf file>")
	sys.exit(1)

bamfile = sys.argv[1]
gtffile = sys.argv[2]

# process gtf file first
gtf_fp = open(gtffile, "r")

gene_elements = defaultdict()

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
		gene_elements[gene_id] = {
				'chrom':chrom,
				'start':start,
				'end':end,
				'strand':strand,
				'exons':list(),
				'intronic':list(),
				'reads':list()
				}
	elif feature == 'exon':
		gene_elements[gene_id]['exons'].append({
				'start':start,
				'end':end
				})
	

gtf_fp.close()

bam_fp = pysam.AlignmentFile(bamfile, "rb")

if not bam_fp.check_index():
	print("ERROR: no index!", file=sys.stderr)
	sys.exit(1)

for gene in gene_elements:
	overlap_reads = list()

	for read in bam_fp.fetch(gene_elements[gene]['chrom'], gene_elements[gene]['start'], gene_elements[gene]['stop']):
		if ( read.is_secondary or 
			read.is_duplicate or 
			read.is_qcfail or 
			read.is_unmapped or 
			read.mate_is_unmapped or  
			not read.is_paired or 
			not read.is_proper_pair ):

			print("ERROR: read fails. Continuing!", file=sys.stderr)
			continue

		overlap_reads.append(read)
	
	for overlap in overlap_reads:
		print("DEBUG:", gene, " my read is:\n[", overlap, "]", file=sys.stderr)




#		print(line)
#	line_list = line.split('\t')
#	read_name = line_list[0]
#	read_flag = line_list[1]
#	read_start_pos = line_list[3]
#	read_rnext = line_list[6]
#	read_pnext = line_list[7]

#	if (
#		(not read_flag & 0x2) or    # unmapped
#		(read_flag & 0x100) or      # multimapped
#		(read_flag & 0x200) or      # low quality
#		(read_flag & 400) or        # PCR duplicate
#		(not read_flag & 0x1)       # unpaired
#	): 
#		continue
