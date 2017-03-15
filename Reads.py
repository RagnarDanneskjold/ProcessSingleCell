def parseFragment(readPair):
    fragment = {
        'strand':'',
        'start':0,
        'end':0
    }
 
    if readPair[0].is_reverse == readPair[1].is_reverse:
        print("ERROR READS ARE SAME DIRECTION",file=sys.stderr)
        return {}

    read1 = readPair[0]
    read2 = readPair[1]

    if readPair[1].is_read1:
        read1 = readPair[1]
        read2 = readPair[0]

    if not read1.is_reverse:
        fragment['start'] = read1.reference_start
        fragment['end'] = read2.reference_end
        fragment['strand'] = '+'
    else:
        fragment['start'] = read2.reference_start
        fragment['end'] = read1.reference_end
        fragment['strand'] = '-'

    return fragment

def readQualityCheck(read,badreads):
    if ( read.is_secondary or 
#        read.is_duplicate or 
        read.is_qcfail or 
        read.is_unmapped or 
        read.mate_is_unmapped or  
        not read.is_paired or 
        not read.is_proper_pair or
        read.query_name in badreads
        ):
        badreads[read.query_name] = True
        return False

    return True

