from Overlaps import overlap

def rangeBsearch(queryStart, queryEnd, data):
    if len(data) == 0:
        return []

    upper_bound = rangeBsearchUpper(queryStart, queryEnd, data)
    lower_bound = rangeBsearchLower(queryStart, queryEnd, data)

#    print("upper:",upper_bound, "lower:",lower_bound, "data len:",len(data))

    #print("queryStart", queryStart, "queryEnd", queryEnd)
    #print("lower bound:",data[lower_bound])

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

        if data[mid]['start'] >= queryEnd:
            upper = mid - 1
        else:
            lower = mid + 1
        
    return lower

def rangeBsearchLower(queryStart, queryEnd, data):
    upper = len(data) - 1
    lower = 0

    while lower <= upper:
        mid = (lower + upper) // 2
#        print(mid, "mid[start]:",data[mid]['start'],"mid[end]:", data[mid]['end'], "queryStart:",queryStart, "queryEnd:",queryEnd)
        if data[mid]['end'] <= queryStart:
            lower = mid + 1
        else:
            upper = mid - 1
 
# this won't work as long as we have overlaps
#    while( lower > 0 && overlap(data[lower-1], {'start':queryStart, 'end':queryEnd}) ):
#        lower -= 1

    return lower

