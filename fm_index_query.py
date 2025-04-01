# Burrows-Wheeler algorythm implemented for string matching through FM indexing

def BWT(T, DoReturnOffsets=False):
    # Initialize by adding an end character if it isn't already there
    if T[-1] != "$":
        T += "$"
    n = len(T)
    rotations = sorted(range(n), key=lambda i: T[i:]+T[:i]) # Store the starting index of each rotation, used later to find the offset of a query match
    L = "".join(T[i-1 % n] for i in rotations) # All values of i-1 are less than the total length n, so the modulo returns i-1, except for -1 which returns n-1, effectively looping around to the last character $ when i=0
    if DoReturnOffsets:
        return L, rotations
    else:
        return L

def FMIndexQuery(T, P):
    # INITIALIZATION
    L, offsets = BWT(T, True) # Start by converting T to BWT(T) for FM-indexing
    char_count = {char:L.count(char) for char in set(L)} # Number of times each character shows up in the seq
    char_Fpos = {char:sum(char_count[p] for p in set(L) if p < char) for char in set(L)} # For each character, sum the counts of all characters that precede it in the seq's alphabet, which results in the index of the character's first occurrence (0th rank) in the F column of the BWM
    # FM INDEX QUERY
    if P[-1] in T:
        prefix_indexes = [x for x in range(char_Fpos[P[-1]], char_Fpos[P[-1]]+char_count[P[-1]])] # Initial range is indexes of F that match the last character of P => All rotations of T that have P[-1] as a prefix
    else: 
        return None # If the last character is not even present in T, the search doesn't start, as the initial range would be empty
    query_char = -2
    while query_char >= (-len(P)):
        c = P[query_char]
        matching_prefix_indexes = [e for e in prefix_indexes if L[e] == c] # The range is restricted by querying for the previous character in the query string P
        prefix_indexes = matching_prefix_indexes
        if len(prefix_indexes) == 0: # Stop the search if no rotations have P's suffix as their prefix
            break
        bottom_rank = L[:prefix_indexes[0]].count(c)
        bottom_index = char_Fpos[c] + bottom_rank # FL mapping of the first matching character in L to F, all other characters follow
        top_index = bottom_index + len(prefix_indexes) # These two indexes define the range for the next query
        prefix_indexes = [x for x in range(bottom_index, top_index)]
        query_char -= 1
    return len(prefix_indexes), sorted([offsets[i] for i in prefix_indexes]) # Return number of hits + offset of each hit within T
