import pysam

## CONSTANTS ##
LINES_TO_CHECK = 20
UNSORTED_PERC = 0.5
##

def checkBam(filename):
    test_fp = pysam.AlignmentFile(filename, "rb")

    pair_count = 0
    prev_name = ""
    unsorted_limit = int(UNSORTED_PERC * LINES_TO_CHECK)
    unsort_count = 0

    for algn in test_fp.head(LINES_TO_CHECK):
        curr_name = algn.query_name

        if curr_name == prev_name:
            pair_count += 1
        elif prev_name != "" and pair_count != 2:
            unsort_count += 1
                
        if unsort_count >= unsorted_limit:
            test_fp.close()
            return False

            prev_name = curr_name

    test_fp.close()

    return True

