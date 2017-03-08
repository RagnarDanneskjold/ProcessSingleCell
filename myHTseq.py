import sys
import pysam
import re
import operator
from collections import defaultdict

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
				parsedData[chrom][strand][i]['exons'].sort(key=operator.itemgetter('start', 'end'))

	return parsedData

def overlap(block1, block2):
	if ( block1['start'] < block2['end'] and 
			block1['end'] > block2['start']
			):
		return True

	return False

def rangeBsearch(queryStart, queryEnd, data):
	upper_bound = rangeBsearchUpper(queryStart, queryEnd, data)
	lower_bound = rangeBsearchLower(queryStart, queryEnd, data)

	if (lower_bound >= len(data) or
			not overlap(data[lower_bound], {
				'start':queryStart,
				'end':queryEnd
				})
			):
		return []

	return data[lower_bound:upper_bound]


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



## END FUNCTIONS ##

if (len(sys.argv) < 2):
	print("usage:", sys.argv[0], "<name-sorted bamfile> <GENCODE gtf file>")
	sys.exit(1)

bamfile = sys.argv[1]
gtffile = sys.argv[2]

print ("DEBUG: opening GTF file", file=sys.stderr)
gene_elements = parseGTFFile(gtffile)
print("DEBUG: done with GTF file", file=sys.stderr)



"""
bam_fp = pysam.AlignmentFile(bamfile, "rb")

print("DEBUG: checking BAM index", file=sys.stderr)
if not bam_fp.check_index():
	print("ERROR: no index!", file=sys.stderr)
	sys.exit(1)


print("DEBUG: starting to iterate over genes", file=sys.stderr)
for gene in gene_elements:
	print("DEBUG: gene is", gene, file=sys.stderr)
	overlap_reads = dict()

	for read in bam_fp.fetch(gene_elements[gene]['chrom'], gene_elements[gene]['start'], gene_elements[gene]['end']):
		if ( read.is_secondary or 
			read.is_duplicate or 
			read.is_qcfail or 
			read.is_unmapped or 
			read.mate_is_unmapped or  
			not read.is_paired or 
			not read.is_proper_pair ):
			continue

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
					print("\t\t\tmate fails QC:")

					if mate_read.is_secondary:
						print("\t\t\t\tis_secondary")
					elif mate_read.is_duplicate:
						print("\t\t\t\tis_duplicate")
					elif mate_read.is_qcfail:
						print("\t\t\t\tis_qcfail")
					elif mate_read.is_unmapped:
						print("\t\t\t\tis_unmapped")
					elif mate_read.mate_is_unmapped:
						print("\t\t\t\tmate_is_unmapped")
					elif not mate_read.is_paired:
						print("\t\t\t\tNOT is_paired")
					elif not mate_read.is_proper_paired:
						print("\t\t\t\tNOT is_proper_paired")
					else:
						print("\t\t\t\tSomething else?")

				else:
					print("\t\t\tmate passes QC")
					print("\t\t\tread aligned to:", read.get_reference_positions())
					print("\t\t\tmate aligned to:", mate_read.get_reference_positions())
					print("\t\t\tgene positions: chr", gene_elements[gene]['chrom'], "start", gene_elements[gene]['start'], "end", gene_elements[gene]['end'])
					print("\t\t\tread overlap is:", read.get_overlap(gene_elements[gene]['start'], gene_elements[gene]['end']))
					print("\t\t\tmate overlap is:", mate_read.get_overlap(gene_elements[gene]['start'], gene_elements[gene]['end']))
					


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
	
		if feature == 'gene':
			parsedData[gene_id] = {
					'chrom':chrom,
					'start':start,
					'end':end,
					'strand':strand,
					'exons':list(),
					'intronic':list(),
					'reads':list()
					}
		elif feature == 'exon':
			parsedData[gene_id]['exons'].append({
					'start':start,
					'end':end
					})
"""
