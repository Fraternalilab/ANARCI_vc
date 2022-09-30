from .schemes import *
# Alphabet used for insertion (last (-1th) is a blank space for no insertion)
alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "AA", "BB", "CC", "DD", "EE", "FF", "GG", "HH", "II", "JJ", "KK", "LL", "MM", "NN", "OO", "PP", "QQ", "RR", "SS", "TT", "UU", "VV", "WW", "XX", "YY", "ZZ", " "]

def _number_regions_c(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions):
    """
    General function to number a sequence and divide it into different regions

    @param sequence: The sequence string
    @param state_vector: The list of states from the aligned hmm
    @param state_string: A string of states for the scheme relative to IMGT (this is X for a direct equivalence, I if needs to be treated as insertion)
    @param region_string: A string of characters that indicate which hmm states are in each regions for this scheme (i.e. how should the sequence be divided up)
    @param region_index_dict: A dictionary converting the characters in region string to an index of the regions.
    @param rels: The difference of the numbering integer at the *start* of each region
    @param n_regions: The number of regions
    @param exclude_deletions: A list of region indices for which deletion states should not be included. Typically the CDRs.
                              These will be reannotated in the scheme function. Also allows the reset of insertions.

    @return: A list of lists where each region has been numbered according to the scheme. Some regions will need renumbering. This should be taken care of after the function called.

    """

    #state_vector = smooth_insertions( state_vector ) # The smoothing function is currently not needed for the C region.

    _regions = [ [] for _ in range(n_regions) ]

    # Initialise the insertion index (-1 is a blank space) and the previous state.
    insertion = -1
    previous_state_id = 1
    previous_state_type = 'd'
    start_index, end_index  = None, None

    region = None

    # Iterate over the aligned state vector
    for (state_id, state_type ), si in state_vector:

        # Retrieve the region index
        region = region_index_dict[region_string[state_id-1]]

        if region==0:
            # treat region 0 as all insersion of number 1
            _regions[region].append(((1,alphabet[8-state_id]),sequence[si]))
        else:
            # Check the state_types
            if state_type == "m": # It is a match

                # Check whether this position is in the scheme as an independent state
                if state_string[state_id-1]=="I": # No, it should be treated as an insertion
                    if previous_state_type != 'd': # Unless there was a deletion beforehand in which case this should be a real pos.
                        insertion +=1 # Increment the insertion annotation index

                    rels[region] -= 1 # Update the relative numbering from the imgt states
                else: # Yes
                    insertion = -1 # Reset the insertions


                # Add the numbering annotation to the appropriate region list
                _regions[region].append( ( (state_id + rels[region], alphabet[insertion] ), sequence[si]  ) )

                previous_state_id = state_id # Record the previous state ID
                if start_index is None:
                    start_index = si
                end_index = si

                previous_state_type = state_type

            elif state_type == "i": # It is an insertion
                insertion +=1 # Increment the insertion annotation index

                # Add the numbering annotation to the appropriate region list
                _regions[region].append( ( (previous_state_id + rels[region], alphabet[insertion]), sequence[si]  ) )
                if start_index is None:
                    start_index = si
                end_index = si

                previous_state_type = state_type

            else: # It is a deletion
                previous_state_type = state_type

                # Check whether this position is in the scheme as an independent state
                if state_string[state_id-1]=="I": # No, therefore irrelevant to the scheme.
                    rels[region] -= 1 # Update the relative numbering from the imgt states
                    continue

                insertion = -1 # Reset the insertions
                previous_state_id = state_id # Record the previous state ID, should not be needed (no delete to insert state transition)


        # Reset the inssertion index if necessary and allowed. (Means the insertion code is meaningless and will be reannotated)
        if insertion >= 25 and region in exclude_deletions:
            insertion = 0

        assert insertion < 25, "Too many insertions for numbering scheme to handle" # We ran out of letters.

    return _regions, start_index, end_index

def gap_missing_c( numbering ):
    '''
    Place gaps when a number is missing. All except wolfguy are continuously numbered
    '''
    # Gaps placed where a number is not present
    num = [ ((0,' '),'-') ]
    for p, a in sum( numbering, [] ):
        if p[0] > num[-1][0][0]+1:
            for _i in range( num[-1][0][0]+1, p[0] ):
                if _i not in [32,33]+list(range(46,77)):
                    num.append( ((_i, ' '), '-' ) )
        num.append( (p,a) )
    return num[1:]

def number_imgt_c(state_vector, sequence):
    """
    Apply the IMGT numbering scheme for heavy or light chains

    Rules should be implemented using two strings - the state string and the region string.

    There are 128 states in the HMMs. Treat X as a direct match in IMGT scheme, I is an insertion. (All X's for IMGT)
    XXXXXXXX XXXXXXXXXXXXXXX III XXXXXXXXXXX XXXXXXXXXX XXXXXXX IIIIIII XXXXXXXX IIIIIIIIIIIIII XXXXXXXXXXXX II XXXXXXXX XXXXXXXXXXXXX XXXXXXXXXX
    00000000 111111111111111 222 33333333333 4444444444 5555555 6666666 77777777 88888888888888 999999999999 aa bbbbbbbb ccccccccccccc dddddddddd

    Regions:
    0: before A strand ( should be numbered as 1.x)
    1: A strand
    2: AB turn
    3: B strand
    4: BC loop
    5: C strand
    6: CD turn
    7: D strand
    8: DE turn
    9: E strand
    a: EF turn
    b: F strand
    c: FG loop
    d: G strand

    """

    # Set up the numbering

    # State string - 'X' means the imgt position exists in the scheme. 'I' means that it should be treated as an insertion of the previous number
    state_string =  'XXXXXXXXXXXXXXXXXXXXXXXIIIXXXXXXXXXXXXXXXXXXXXXXXXXXXXIIIIIIIXXXXXXXXIIIIIIIIIIIIIIXXXXXXXXXXXXIIXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    # Or could consider to change the first character into X, and the numbering starts at 0, and then do some modification after that.
    # length:128

    # Region string - regions that should be treated separately in putting the numbering together
    region_string = '00000000111111111111111222333333333334444444444555555566666667777777788888888888888999999999999aabbbbbbbbcccccccccccccddddddddddd'

    region_index_dict = {
        "0":0,
        "1":1,
        "2":2,
        "3":3,
        "4":4,
        "5":5,
        "6":6,
        "7":7,
        "8":8,
        "9":9,
        "a":10,
        "b":11,
        "c":12,
        "d":13
                         }

    # Define how the scheme's numbering differs from IMGT at the start of each region.
    # This is updated in the loop below
    rels              =  {0:0,
                          1:-8,
                          2:-8,
                          3:-11,
                          4:-11,
                          5:-9,
                          6:-9,
                          7:15,
                          8:15,
                          9:1,
                          10:1,
                          11:-1,
                          12:-1,
                          13:-1
                          }

    moving=0

    n_regions = 14 # one extra region is the fragment before the A strand (numbered as 1.x)

    #exclude_deletions = [1,3,5]
    exclude_deletions=[]

    # Needs to write another function for the C region: _number_regions_c
    # Since the current function won't allow the region to start with I? -Check this!
    # Or consider to re-write the region string: needs to check the "smooth" function to see if we can do this
    _regions, startindex, endindex = _number_regions_c(sequence, state_vector, state_string , region_string,  region_index_dict, rels, n_regions, exclude_deletions)

    #return _regions
    ###############
    # Renumbering #
    ###############

    _numbering = [ _regions[0], # before A strand
                   _regions[1], # A strand
                   _regions[2], # AB turn
                   _regions[3], # B strand
                   [],          # BC
                   _regions[5], # C strand
                   _regions[6], # CD
                   _regions[7], # D strand
                   [],          # DE
                   _regions[9], # E strand
                   _regions[10], # EF
                   _regions[11],# F strand
                   [],          # FG
                   _regions[13] # G strand

                 ]

    # BC loop: re-numbering: change the number 32& 33 into 34&35, and all the follwing numbering so that BC loop ends at 38
    for (state_id, state_type ), si in _regions[4]:
        n_state_id=state_id
        if state_id >=32:
            n_state_id+=2
        _numbering[4].append( ((n_state_id, state_type),si) )


    # DE loop: symmetric insersion
    deseq    = "".join([ x[1] for x in _regions[8] if x[1] != "-" ])
    delength = len(deseq)
    if delength > 14: return [], startindex, endindex # Too many insertions. Do not apply numbering.
    si = 0
    previous_state_id = 84
    for ann in get_imgt_cdr(delength+2, 2, 84, 86)[1:-1]:
        if ann is None:
            _numbering[8].append( ((previous_state_id+1, " "), "-"   ) )
            previous_state_id+=1
        else:
            _numbering[8].append( (ann, deseq[si] ) )
            previous_state_id = ann[0]
            si+=1


    # FG loop: similar case as CDR3 loop
    # _regions[12]
    fg_seq    = "".join([ x[1] for x in _regions[12] if x[1] != "-" ])
    fg_length = len(fg_seq)
    if fg_length > 30: return [], startindex, endindex
    # Too many insertions. Do not apply numbering. Because the longest FG loop is 25 in length, according to IMGT

    si = 0
    previous_state_id = 104
    for ann in get_imgt_cdr(fg_length, 13, 105, 118):
        if ann is None:
            _numbering[12].append( ((previous_state_id+1, " "), "-"   ) )
            previous_state_id+=1
        else:
            _numbering[12].append( (ann, fg_seq[si] ) )
            previous_state_id = ann[0]
            si+=1



    #return _numbering
    return gap_missing_c( _numbering ), startindex, endindex
