import re
import operator
from Tree import GeneTree

def parseValidLine(line):
    line_ar = line.rstrip('\n\r').split('\t')

    if len(line_ar) < 9:
        return {}

    id_string = line_ar[8]
    gene_id = re.search(r'gene_id \"(.+?)\";', id_string).group(1)

    result = {
            'chrom': line_ar[0],
            'feature': line_ar[2],
            'start': int(line_ar[3]),
            'end': int(line_ar[4]),
            'strand': line_ar[6],
            'gene_id': gene_id
        }

    if result['feature'] != "gene" and result['feature'] != "exon":
        return {}

    return result

def parseGTFFile (gtffile, bam_num):
    gtf_fp = open(gtffile,"r")
    parsedData = dict()
    curr_gene = 0

    for line in gtf_fp:
        if line.startswith('#'):
            continue

        fields = parseValidLine(line)

        if not fields:
            continue
        
        if fields['chrom'] not in parsedData:
            parsedData[fields['chrom']] = dict()
            parsedData[fields['chrom']]['+'] = GeneTree(fields['chrom'],'+')
            parsedData[fields['chrom']]['-'] = GeneTree(fields['chrom'],'-')

        if fields['feature'] == 'gene':
            curr_gene = parsedData[fields['chrom']][fields['strand']].addNode(fields,fields['gene_id'], bam_num)
        else: # exon
            parsedData[fields['chrom']][fields['strand']].addExon(fields,curr_gene)

    gtf_fp.close()

    for chrom in parsedData:
        for strand in parsedData[chrom]:
            parsedData[chrom][strand].balanceAll()

    return parsedData

