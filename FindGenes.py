import Searches
import Overlaps
from Reads import parseFragment

def runOverlapGenes(prev_reads, bam_fp, genes):
    if len(prev_reads) != 2:
        return False

    frag = parseFragment(prev_reads)
    chrom = bam_fp.get_reference_name(prev_reads[0].reference_id)

    if frag != {} and chrom in genes:
        result = overlapGenes(genes, frag, chrom)
        return True

    return False

def overlapGenes(gene_elements, frag, chrom):
    overlap_indicies = gene_elements[chrom][frag['strand']]['tree'].findNode(frag)

    if overlap_indicies == [] or len(overlap_indicies) > 1:
        return False

    test = gene_elements[chrom][frag['strand']]['genes'][overlap_indicies[0]]['exons'].findNodeBool(frag)
    if test:
        gene_elements[chrom][frag['strand']]['genes'][overlap_indicies[0]]['reads'] += 1
        return True

    return False
