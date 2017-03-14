import re
import operator
import sys # for sys.stderr
from Tree import MyTree

def parseValidLine(line):
    line_ar = line.rstrip('\n\r').split('\t')

    gene_of_interest = False

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

def parseGTFFile (gtf_fp):
    parsedData = dict()

    count = 0
    for line in gtf_fp:
        count+=1
        if line.startswith('#'):
            continue

        fields = parseValidLine(line)

        if not fields:
            continue
        
        if fields['chrom'] not in parsedData:
            parsedData[fields['chrom']] = {
                '+': {
                    'tree' : MyTree(),
                    'genes' : list()
                    },
                '-': {
                    'tree' : MyTree(),
                    'genes' : list()
                    }
                }

        if fields['feature'] == 'gene':
            curr_index = len(parsedData[fields['chrom']][fields['strand']]['genes'])

            parsedData[fields['chrom']][fields['strand']]['tree'].addNode(fields, curr_index)

            parsedData[fields['chrom']][fields['strand']]['genes'].append({
                'gene_id':fields['gene_id'],
                'start':fields['start'],
                'end':fields['end'],
                'exons': MyTree(),
                'reads': 0
                })

        else: # exon
            curr_index = len(parsedData[fields['chrom']][fields['strand']]['genes']) - 1
            
            # the index doesn't matter
            parsedData[fields['chrom']][fields['strand']]['genes'][curr_index]['exons'].addNode(fields, curr_index)
   
    print("balancing", file=sys.stderr)

    for chrom in parsedData:
        for strand in parsedData[chrom]:
            parsedData[chrom][strand]['tree'].balance()
            for i in range(len(parsedData[chrom][strand]['genes'])):
                parsedData[chrom][strand]['genes'][i]['exons'].balance()

    print("done balancing", file=sys.stderr)

    return parsedData
