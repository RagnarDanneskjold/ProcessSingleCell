from Reads import parseFragment

def runOverlapGenes(prev_reads, bam_fp, genes):
    if len(prev_reads) != 2:
        print("not_pair_reads!")
        return False

    frag = parseFragment(prev_reads)


#    if frag['end'] - frag['start'] > 600 or frag['end'] - frag['start'] < 50:
#        print("frag_length_error!")
#        return False

    chrom = bam_fp.get_reference_name(prev_reads[0].reference_id)

    if frag != {} and chrom in genes:
#        result = overlapGenes(genes, frag, chrom)
#        return True
        return genes[chrom][frag['strand']].overlapInterval(frag,prev_reads)       

    print("frag_or_chrom_read_error!")
    return False

#def overlapGenes(gene_elements, frag, chrom):
#    overlap_indicies = gene_elements[chrom][frag['strand']]['tree'].findNode(frag)
#
#    if overlap_indicies == [] or len(overlap_indicies) > 1:
#    if overlap_indicies == []:
#        print("not_overlap!")
#        return False
#    elif len(overlap_indicies) > 1:
#        print("multiple_overlaps!")
#        return False
#
#    test = gene_elements[chrom][frag['strand']]['genes'][overlap_indicies[0]]['exons'].findNodeBool(frag)
#    if test:
#        print("SUCCESS_read_added!")
#        gene_elements[chrom][frag['strand']]['genes'][overlap_indicies[0]]['reads'] += 1
#        return True
#
#    print("no_exons_overlap!")
#    return False
