from Overlaps import overlap

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

        if data[mid]['end'] <= queryStart:
            lower = mid + 1
        else:
            upper = mid - 1
        
    return lower

