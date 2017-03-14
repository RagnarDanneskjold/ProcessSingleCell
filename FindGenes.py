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
#        print("Found overlap:",result)
        return True

    return False

def overlapGenes(gene_elements, frag, chrom):
#    print(chrom,frag)
#    overlap_indicies = Searches.rangeBsearch(frag['start'], frag['end'],
#        gene_elements[chrom][frag['strand']])
#    print("overlapGenes()")

    overlap_indicies = gene_elements[chrom][frag['strand']]['tree'].findNode(frag)

#    print(overlap_indicies)

#    print("overlap_indicies:",overlap_indicies)
    if overlap_indicies == [] or len(overlap_indicies) > 1:
        return False

#    low_index = overlap_indicies[0]
#    high_index = overlap_indicies[1]
#    overlap_gene_index = low_index

#    if (high_index - low_index) > 1:
#        return False
#        curr_highest = 0
#        found_equal = False

#        for i in range(low_index, high_index):
#            score = Overlaps.overlapBases(gene_elements[chrom][frag['strand']][i],frag)
#            if score > curr_highest:
#                found_equal = False
#                overlap_gene_index = i
#                curr_highest = score
#            elif score == curr_highest:
#                found_equal = True
#
#        if found_equal:
#            return False

#    if overlapExons(gene_elements[chrom][frag['strand']]['genes'][overlap_indicies[0]]['exons'], frag):
        if gene_elements[chrom][frag['strand']]['genes'][overlap_indicies[0]]['exons'].findNodeBool(frag):
            gene_elements[chrom][frag['strand']]['genes'][overlap_indicies[0]]['reads'] += 1
            return True

    return False

#def overlapExons(gene, frag):
#    for exon in gene:
#        if Overlaps.overlap(exon,frag):
#            return True
#
#    return False
