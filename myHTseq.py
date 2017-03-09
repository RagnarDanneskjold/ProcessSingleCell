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

def overlapBases(block1, block2):
	if not overlap(block1, block2):
		return 0

	if block1['start'] < block2['start']:
		if block1['end'] < block2['end']:
			return block1['end'] - block2['start']
		else:
			return block2['end'] - block2['start']
	else:
		if block1['end'] < block2['end']:
			return block1['end'] - block1['start']
		else:
			return block2['end'] - block1['start']

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

	return [lower_bound, upper_bound]
#	return data[lower_bound:upper_bound]


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
	return fragment


## END FUNCTIONS ##

if (len(sys.argv) < 2):
	print("usage:", sys.argv[0], "<name-sorted bamfile> <GENCODE gtf file>")
	sys.exit(1)

bamfiles_list_file = sys.argv[1]
gtffile = sys.argv[2]

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
		if ( read.is_secondary or 
			read.is_duplicate or 
			read.is_qcfail or 
			read.is_unmapped or 
			read.mate_is_unmapped or  
			not read.is_paired or 
			not read.is_proper_pair ):
			continue
		
		if read.query_name == prev_read_name:
			prev_reads.append(read)
		else:
			if (len(prev_reads) > 0 and
					"chr" + str(prev_reads[0].reference_id) in gene_elements
					):
				# handle last pair
				# FIXME find a way to handle the last pair in the entire bam file
				# debug:
				if len(prev_reads) != 2:
					print ("ERROR! len(prev_reads) is",len(prev_reads), file=sys.stderr)
					for err_reads in prev_reads:
						print(err_reads)
						bam_fp.close()
						sys.exit(1)

				# get start, end, strand of the FRAGMENT

				frag = parseFragment(prev_reads)

				if frag == {}:
					# replace with this read
					prev_read_name = read.query_name
					prev_reads = [read]
					continue
			#parsedData[chrom][strand].append({
			#		'gene_id':gene_id,
			#		'start':start,
			#		'end':end,
			#		'exons':list(),
			#		'intronic':list(),
			#		'reads':list()
			#	})
				curr_chrom = prev_reads[0].reference_id
				curr_chrom = "chr" + str(curr_chrom)
				#print(curr_chrom, file=sys.stderr)

				overlap_indicies = rangeBsearch(frag['start'], frag['end'], 
						gene_elements[curr_chrom][frag['strand']])
				if overlap_indicies == []:
					prev_read_name = read.query_name
					prev_reads = [read]
					continue


				if overlap_indicies[0] - overlap_indicies[1] > 1:
					# this read overlaps multiple elements
					highest_overlap_index = []
					curr_highest = 0

					for i in range(overlap_indicies[0], overlap_indicies[1]):
						score = overlapBases(gene_elements[curr_chrom][frag['strand']][i],frag)
						if score >= curr_highest:
							highest_overlap_index.append(i)
							curr_highest = score

					if len(highest_overlap_index) > 1:
						# multiple with equal overlap
						# cannot assign?
						pass # FIXME temporary
					elif len(highest_overlap_index) == 1:
						gene_elements[curr_chrom][frag['strand']][highest_overlap_index[0]]['reads'].append(frag)
					# else overlap for all was 0, continue
				else: # 1 element
					gene_elements[curr_chrom][frag['strand']][overlap_indicies[0]]['reads'].append(frag)

# FIXME handle exons/introns later
			# replace with this read
			prev_read_name = read.query_name
			prev_reads = [read]




	bam_fp.close()


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
