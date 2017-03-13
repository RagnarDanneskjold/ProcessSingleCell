import Searches
import Overlaps
from Reads import parseFragment

def runOverlapGenes(prev_reads, bam_fp, genes):
    if len(prev_reads) != 2:
        return False

    frag = parseFragment(prev_reads)
    chrom = bam_fp.get_reference_name(prev_reads[0].reference_id)

    if frag != {} and chrom in genes:
        FindGenes.overlapGenes(genes, frag, chrom)
        return True

    return False

def overlapGenes(gene_elements, frag, chrom):
#    if len(prev_reads) != 2 or chrom not in gene_elements or frag == {}:
#        return
#    frag = parseFragment(prev_reads)

    overlap_indicies = rangeBsearch(frag['start'], frag['end'],
        gene_elements[chrom][frag['strand']])

    if overlap_indicies == []:
        return False

    low_index = overlap_indicies[0]
    high_index = overlap_indicies[1]
    overlap_gene_index = low_index

    if (high_index - low_index) > 1:
        curr_highest = 0
        found_equal = False

        for i in range(low_index, high_index):
            score = overlapBases(gene_elements[chrom][frag['strand']][i],frag)
            if score > curr_highest:
                found_equal = False
                overlap_gene_index = i
                curr_highest = score
            elif score == curr_highest:
                found_equal = True

        if found_equal:
            return False

    gene_elements[chrom][frag['strand']][overlap_gene_index]['reads'].append(frag)

# FIXME handle exons/introns later
    return True

