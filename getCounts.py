import HTSeq
import sys

gtffile = sys.argv[1]
listOfBams = sys.argv[2:]

#print(listOfBams)

exons = HTSeq.GenomicArrayOfSets("auto", stranded=True)
gtf = HTSeq.GFF_Reader(gtffile, end_included=True)

for feature in gtf:
	if feature.type == "exon":
		exons[ feature.iv ] += feature.name

for bamfile in listOfBams:
	bamObj = HTSeq.BAM_Reader(bamfile)
	for alignment in bamObj:
		if alignment.aligned:
			iset = None
			for interval2, step_set in exons[ alignment.iv ].steps():
				if iset is None:
					iset = step_set.copy()
				else:
					iset.intersection_update( step_set )

				if len(iset) == 1:
					counts[ list(iset)[0] ] += 1

for name in sorted(counts.keys()):
	print (name, counts[name])
