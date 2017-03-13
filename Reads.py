def parseFragment(readPair):
    fragment = {
        'strand':'',
        'start':0,
        'end':0
    }
 
    if readPair[0].is_reverse == readPair[1].is_reverse:
        print("ERROR READS ARE SAME DIRECTION",file=sys.stderr)
        return {}

    if readPair[0].is_read1:
        if not readPair[0].is_reverse:
            fragment['start'] = readPair[0].reference_start
            fragment['end'] = readPair[1].reference_end
            fragment['strand'] = '+'
        else:
            fragment['start'] = readPair[1].reference_start
            fragment['end'] = readPair[0].reference_end
            fragment['strand'] = '-'
    elif readPair[1].is_read1:
        if not readPair[1].is_reverse:
            fragment['start'] = readPair[1].reference_start
            fragment['end'] = readPair[0].reference_end
            fragment['strand'] = '+'
        else:
            fragment['start'] = readPair[0].reference_start
            fragment['end'] = readPair[1].reference_end
            fragment['strand'] = '-'
        
    return fragment

def readQualityCheck(read):
    if ( read.is_secondary or 
        read.is_duplicate or 
        read.is_qcfail or 
        read.is_unmapped or 
        read.mate_is_unmapped or  
        not read.is_paired or 
        not read.is_proper_pair ):
        return False

    return True

