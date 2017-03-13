import Bam
import sys

def getInput(args):
    if len(args) < 3:
        return {}
    
    inputdict = {
        'bamfile':args[1],
        'gtffile':args[2],
        'outfile':args[3]
        }

    if not Bam.checkBam(inputdict['bamfile']):
        print("ERROR: bamfile did not check out!", file=sys.stderr)
        return {}

    # FIXME check if gtf exists

    return inputdict


def printOutput(outfilename, gene_elements):
    outfile = open(outfilename, 'w')
    #outfile.write('chr\tstart\tend\tstrand\tgene_id\tcount\n')

    #FIXME handle multiple counts for different bamfiles!
    for chrom in gene_elements:
        for strand in gene_elements[chrom]:
            for gene in gene_elements[chrom][strand]:
                outstr = chrom + '\t' + str(gene['start']) + '\t' +\
                        str(gene['end']) + '\t' + strand + '\t' +\
                        gene['gene_id'] + '\t' + str(len(gene['reads'])) +\
                        '\n'

                outfile.write(outstr)

    outfile.close()
