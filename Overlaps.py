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
            return block1['end'] - block2['start']  
        else:
            return block2['end'] - block2['start'] 
    else:
        if block1['end'] < block2['end']:
            return block1['end'] - block1['start'] 
        else:
            return block2['end'] - block1['start'] 
