import sys
import pysam
import re
from collections import defaultdict

if (len(sys.argv) < 2):
	print("usage:", sys.argv[0], "<name-sorted bamfile> <GENCODE gtf file>")
	sys.exit(1)

bamfile = sys.argv[1]
gtffile = sys.argv[2]

# process gtf file first
gtf_fp = open(gtffile, "r")

gene_elements = defaultdict()

print ("DEBUG: opening GTF file", file=sys.stderr)

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

print("DEBUG: done with GTF file", file=sys.stderr)

bam_fp = pysam.AlignmentFile(bamfile, "rb")

print("DEBUG: checking BAM index", file=sys.stderr)
if not bam_fp.check_index():
	print("ERROR: no index!", file=sys.stderr)
	sys.exit(1)

print("DEBUG: starting to iterate over genes", file=sys.stderr)
for gene in gene_elements:
	print("DEBUG: gene is", gene, file=sys.stderr)
	#overlap_reads = list()
	overlap_reads = dict()

	for read in bam_fp.fetch(gene_elements[gene]['chrom'], gene_elements[gene]['start'], gene_elements[gene]['end']):
		if ( read.is_secondary or 
			read.is_duplicate or 
			read.is_qcfail or 
			read.is_unmapped or 
			read.mate_is_unmapped or  
			not read.is_paired or 
			not read.is_proper_pair ):

#			print("ERROR: read fails. Continuing!", file=sys.stderr)
			continue

		
#		print("DEBUG: appending read", file=sys.stderr)
#		overlap_reads.append(read)
		if read.query_name not in overlap_reads:
			overlap_reads[read.query_name] = list()

		overlap_reads[read.query_name].append(read)
	
	print("DEBUG: overlap_reads num:", len(overlap_reads), file=sys.stderr)
	flag = False
	for qname in overlap_reads:
#		print(qname)
#		for read in overlap_reads[qname]:
#			print("\tChecking a read:")
#			print("\t\tread1?:", read.is_read1)
#			print("\t\tread2?:", read.is_read2)
		if len(overlap_reads[qname]) != 2:
			flag = True
			print("Gene:", gene)
			print("Overlap_reads = ", len(overlap_reads[qname]))
			for read in overlap_reads[qname]:
				print("\tChecking a read:")
				print("\t\tread1?:", read.is_read1)
				print("\t\tread2?:", read.is_read2)
				mate_read = bam_fp.mate(read)
				print("\t\tmate:", mate_read)

				if ( mate_read.is_secondary or 
					mate_read.is_duplicate or 
					mate_read.is_qcfail or 
					mate_read.is_unmapped or 
					mate_read.mate_is_unmapped or  
					not mate_read.is_paired or 
					not mate_read.is_proper_pair ):
					print("\t\t\tmate fails QC")
				else:
					print("\t\t\tmate passes QC")
					print("\t\t\taligned to:", mate_read.get_reference_positions())
					print("\t\t\tgene positions: chr", gene_elements[gene]['chrom'], "start", gene_elements[gene]['start'], "end", gene_elements[gene]['end'])
					


	print("DEBUG: done with gene", gene, file=sys.stderr)

	if flag:
		bam_fp.close()
		sys.exit(1)
bam_fp.close()

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
