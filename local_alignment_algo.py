# Smith and Waterman algorithm for local sequence alignment
        
def assignScore(a, b, match, mismatch):
    if a == b:
        return match
    else:
        return mismatch
    
def traceback(A, B, i, j, M):    
    seqA = ""
    seqcomp = ""
    seqB = ""
    while M[i][j] != "STOP":
        if M[i][j] == "DIAG":
            seqA = A[i-1] + seqA
            seqB = B[j-1] + seqB
            if A[i-1] == B[j-1]:
                seqcomp = "*" + seqcomp
            else:
                seqcomp = "|" + seqcomp
            i -= 1
            j -= 1
        elif M[i][j] == "UP":
            seqA = A[i-1] + seqA
            seqB = "-" + seqB
            seqcomp = " " + seqcomp
            i -= 1
        elif M[i][j] == "LEFT":
            seqA = "-" + seqA
            seqB = B[j-1] + seqB
            seqcomp = " " + seqcomp
            j -= 1
    return [seqA, seqcomp, seqB]

def localAlignment(A, B, gap_pen, match, mismatch):
    m = len(A)                               # number of rows
    n = len(B)                               # number of columns
    sM = [[0]]                               # score matrix
    cM = [["STOP"]]                          # case matrix for traceback
    max_score = 0                            # stores the highest score in the matrix
    max_i = 0                                # the row (i) and column (j) of the highest score 
    max_j = 0
    # MATRIX INITIALIZATION
    for i in range(1, m+1):
        sM.append([0])
        cM.append(["STOP"])
    for j in range(1, n+1):
        sM[0].append(0)
        cM[0].append("STOP")
    # MATRIX FILLING
    for i in range(1, m+1):
        for j in range(1, n+1):
            v1 = sM[i-1][j-1] + assignScore(A[i-1], B[j-1], match, mismatch)
            v2 = sM[i-1][j] + gap_pen
            v3 = sM[i][j-1] + gap_pen
            best_score = max(0, v1, v2, v3)
            sM[i].append(best_score)
            if best_score > max_score:
                max_score = best_score
                max_i, max_j = i, j
            if best_score == 0:
                cM[i].append("STOP")
            elif best_score == v1:
                cM[i].append("DIAG")
            elif best_score == v2:
                cM[i].append("UP")
            elif best_score == v3:
                cM[i].append("LEFT")
    # TRACEBACK
    alignment = traceback(A, B, max_i, max_j, cM)
    return alignment, max_score


