def overlap(block1, block2):
    if ( block1['start'] < block2['end'] and 
        block1['end'] > block2['start']
        ):
        return True

    return False

def overlapBases(block1, block2):
    if not overlap(block1, block2):
        return 0

    if block1['start'] < block2['start']:
        if block1['end'] < block2['end']:
            return block1['end'] - block2['start'] - 1 
        else:
            return block2['end'] - block2['start'] - 1
    else:
        if block1['end'] < block2['end']:
            return block1['end'] - block1['start'] - 1
        else:
            return block2['end'] - block1['start'] - 1

def rangeBsearch(queryStart, queryEnd, data):
    if len(data) == 0:
        return []

    upper_bound = rangeBsearchUpper(queryStart, queryEnd, data)
    lower_bound = rangeBsearchLower(queryStart, queryEnd, data)

    if (lower_bound >= len(data) or
        not overlap(data[lower_bound], {
            'start':queryStart,
            'end':queryEnd
            })
        ):
        return []

    return [lower_bound, upper_bound]


def rangeBsearchUpper(queryStart, queryEnd, data):
    upper = len(data) - 1
    lower = 0

    while lower <= upper:
        mid = (lower + upper) // 2

        if data[mid]['start'] > queryEnd:
            upper = mid - 1
        else:
            lower = mid + 1
        
    return lower

def rangeBsearchLower(queryStart, queryEnd, data):
    upper = len(data) - 1
    lower = 0

    while lower <= upper:
        mid = (lower + upper) // 2

        if data[mid]['end'] < queryStart:
            lower = mid + 1
        else:
            upper = mid - 1
        
    return lower

