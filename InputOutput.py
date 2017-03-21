import Bam
import sys
import os.path

def getInputThread(args,inputdict):
    if len(args) < 4 or not os.path.isfile(args[1]) or not os.path.isfile(args[4]):
        return {}

    inputdict['thread'] = args[3]

    bamlist = open(args[4], "r")

    for line in bamlist:
       filename = line.rstrip('\n\r') 
       inputdict['bams'].append(filename)

    bamlist.close()
 
def getInput(args,THREADING):
    inputdict = {
            'gtffile':'',
            'bams':list(),
            'outfile':''
            }

    if THREADING:
        if getInputThread(args,inputdict) == {}:
            return {}
    elif len(args) < 3 or not os.path.isfile(args[1]) or not os.path.isfile(args[3]):
        return {}
    else:
        bamlist = open(args[3], "r")

        for line in bamlist:
            filename = line.rstrip('\n\r') 
            inputdict['bams'].append(filename)

        bamlist.close()
 
    inputdict['gtffile'] = args[1]
    inputdict['outfile'] = args[2]

    if len(set(inputdict['bams'])) != len(inputdict['bams']):
        print("Duplicate bam files in input!",file=sys.stderr)
        return {}

    for bam in inputdict['bams']:
        if not os.path.isfile(bam) or not Bam.checkBam(bam):
            print("ERROR: bamfile[",bam,"] did not check out!", file=sys.stderr)
            return {}

    return inputdict


def printOutput(outfilename, gene_elements, bam_names):
    out_fp = open(outfilename, 'w')

    header = "chr\tstart\tend\tstrand\tgene_id"

    for filename in bam_names:
        header += "\t" + filename

    header += "\n"
    out_fp.write(header)

    for chrom in sorted(gene_elements):
        for strand in sorted(gene_elements[chrom]):
            gene_elements[chrom][strand].writeTree(out_fp)

    out_fp.close()
