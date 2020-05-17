"""
HELPER FUNCTION FOR INDEX()
"""

# pickling our string
import pickle

# handles SAM files/formatting
from simplesam import Reader, Writer, Sam

# object structure for our FM-Indexes
from fmIndex import fmIndex

# object structure for our suffixes
from suffix import Suffix

# object structure for alignment
from alignment import alignment

# object for ranks
from rank import rank

# dump fmIndexes into output, use this in align
def pickleDump(fmIndexes, output):
    f = open(output, "wb")
    pickle.dump(fmIndexes, f)
    f.close()

# load pickled data, use this in align
def pickleLoad(input):
    f = open(input, "rb")
    fmS = pickle.load(f)
    f.close()
    return fmS

# write to output in a very straight-forward string expression
def writeToOutput(fmIndexes, output):

    # open output file
    f = open(output, "w")


    # write it to output
    for fm in fmIndexes:

        # name, sequence, length, suffixArray, first/lastBWM, occTable
        name = fm.name
        seq = fm.seq
        seqLength = str(fm.length)
        suffixArray = fm.suffixArray
        firstBWM = fm.firstBWM
        lastBWM = fm.lastBWM
        occTable = fm.occTable

        # our entire serialized representation that we will write to output
        serialization = ""

        serialization += "Name: " + name + "\n"
        serialization += "Sequence: " + seq + "\n"
        serialization += "Length of sequence: " + seqLength + "\n"
        serialization += "Suffix Array: " + commaArray(suffixArray) + "\n"
        serialization += "First column of BWM: " + firstBWM + "\n"
        serialization += "Last column of BWM: " + lastBWM + "\n"
        serialization += "First Column Occ Table: " + "\n"

        # build Occ Table string
        # not alphabetized, should i?
        occToStringF = ""
        for key in occTable[0]:
            occToStringF += "  " + key
            occToStringF += ": "
            occToStringF += commaArray((occTable[0])[key])
            occToStringF += "\n"

        serialization += occToStringF

        serialization += "Last Column Occ Table: " + "\n"

        occToStringL = ""
        for key in occTable[1]:
            occToStringL += "  " + key
            occToStringL += ": "
            occToStringL += commaArray((occTable[1])[key])
            occToStringL += "\n"

        serialization += occToStringL + "\n\n"

        # write to file (NOT binary)
        f.write(serialization)

    # close because of file leaks
    f.close()

# given a array of ints, returns string representation
def commaArray(arrayOfInts):

    for i in range(len(arrayOfInts)):
        arrayOfInts[i] = str(arrayOfInts[i])

    commaArrayString = "[" + ', '.join(arrayOfInts) + "]"

    return commaArrayString


# given a sequence and its length, return the suffix array
def buildSuffixArray(seq, seqlength, isThisForBWM):

    # suffix array, we will return this at the end
    suffixArray = []

    # array of suffix objects
    allSuffixes = []

    # put all suffixes into data structure
    # should be 0 to seqLength - 1
    for offset in range(seqlength):

        if (isThisForBWM == False):

            # current suffix
            currentSuffix = seq[offset:(seqlength)]

            # create suffix object, this is useful for "duplicate suffixes"
            suffixObj = Suffix(currentSuffix, offset)

            # update allSuffixes
            allSuffixes.append(suffixObj)
        else:

            # current suffix (bwm)
            currentSuffix = seq[offset:(seqlength)] + seq[0:offset]

            # create suffix object, this is useful for "duplicate suffixes"
            suffixObj = Suffix(currentSuffix, offset)

            # update allSuffixes
            allSuffixes.append(suffixObj)
    
    # sort all suffix objects alphabetically into a list
    allSuffixes = sorted(allSuffixes)

    # are we using this method for BWM(), if so let's stop here
    if (isThisForBWM == True):

        # returns list that will be translated into BWM
        return allSuffixes # listAllSuffixes

    # otherwise let's build a suffix array
    else:

        # place into suffixArray
        for suffix in allSuffixes:

            # current key (suffix) and value (offset)
            key = suffix.getSuffix()
            offset = suffix.getOffset()

            # add value (offset)
            suffixArray.append(offset)

        return suffixArray

# return either the first or last col of a BWM
def BWM(seq, seqLength, col):

    # a matrix represented as a list of lists
    BWM = []

    # add all permutations of seq to BWM
    for rotation in buildSuffixArray(seq, seqLength, True):

        # get rotation suffix
        rotation = rotation.getSuffix()
        BWM.append([rotation])

    # col decides whether you will return the first of last colum of BWM
    if (col.lower() == "first"):
        first = ""
        for rotation in BWM:
            first += (rotation[0])[0:1]
        return first

    elif(col.lower() == "last"):
        last = ""
        for rotation in BWM:
            last += (rotation[0])[len(rotation[0]) - 1:len(rotation[0])]
        return last  
    else:
        print("BWM invalid col parameter!")

# returns our occTable
def occTable(seq, seqLength):

    # Key: character
    # Value: occ-table tally represented by an array
    occTable = dict()

    # Key: character
    # Value: occ-table tally represented by an array
    occTableF = dict()

    # Key: character
    # Value: Stores index of last time a tally was recorded
    tallyIndex = dict()

    # first and last columns of BWM
    F = BWM(seq, seqLength, "first")
    L = BWM(seq, seqLength, "last")

    # counter to keep track of what index we're on
    index = 0

    # build occTable for L
    for character in L:
        # if its the first time we're seeing the character
        if not character in occTable:

            # create tallyRank column for this new character
            # initializes tallyRank with all zeros
            occTable[character] = [0] * (seqLength)

            # store tallyRank, later we'll put it back in
            tallyRank = occTable[character]

            # first occurence
            for i in range(index, len(tallyRank)):
                tallyRank[i] = 1

            # update occTable
            occTable[character] = tallyRank

            # update tallyIndex
            tallyIndex[character] = index
               
        # if we've seen this character before, +1 occurence
        else:
            
            # get tallyRank table for current character
            tallyRank = occTable[character]

            # number to fill with
            fillNumber = tallyRank[index - 1] + 1

            # fill with occurence #
            for i in range(index, len(tallyRank)):
                tallyRank[i] = fillNumber

            # update occTable
            occTable[character] = tallyRank

            # update tallyIdex
            tallyIndex[character] = index

        index += 1

    # reset values
    tallyIndex = dict()
    index = 0

    # build occTable for F
    for character in F:

    # if its the first time we're seeing the character
        if not character in occTableF:

            # create tallyRank column for this new character
            # initializes tallyRank with all zeros
            occTableF[character] = [0] * (seqLength)

            # store tallyRank, later we'll put it back in
            tallyRank = occTableF[character]

            # first occurence
            for i in range(index, len(tallyRank)):
                tallyRank[i] = 1

            # update occTable
            occTableF[character] = tallyRank

            # update tallyIndex
            tallyIndex[character] = index
               
        # if we've seen this character before, +1 occurence
        else:
            
            # get tallyRank table for current character
            tallyRank = occTableF[character]

            # number to fill with
            fillNumber = tallyRank[index - 1] + 1

            # fill with occurence #
            for i in range(index, len(tallyRank)):
                tallyRank[i] = fillNumber

            # update occTable
            occTableF[character] = tallyRank

            # update tallyIdex
            tallyIndex[character] = index

        index += 1

    return (occTableF, occTable)




















"""
HELPER FUNCTIONS FOR ALIGN()    
"""

# convert F, L strings to rank objects
# helpful in getInterval()
def FLMap(F, L):

    # dict for keeping rank #s
    ranksMapF = dict()
    ranksMapL = dict()

    # F-Rank list
    rankF = []
    rankL = []

    for i in range(len(F)):

        characterF = F[i]
        characterL = L[i]

        # F
        if not characterF in ranksMapF:

            # add rank
            ranksMapF[characterF] = 0
            addRank = rank(characterF, ranksMapF[characterF])
            rankF.append(addRank)

            # update dict, increment rank
            value = ranksMapF[characterF]
            ranksMapF[characterF] = (value + 1)

        else:

            # add rank
            addRank = rank(characterF, ranksMapF[characterF])
            rankF.append(addRank)

            # update dict, increment rank
            value = ranksMapF[characterF]
            ranksMapF[characterF] = (value + 1)

        # L
        if not characterL in ranksMapL:

            # add rank
            ranksMapL[characterL] = 0
            addRank = rank(characterL, ranksMapL[characterL])
            rankL.append(addRank)

            # update dict, increment rank
            value = ranksMapL[characterL]
            ranksMapL[characterL] = (value + 1)

        else:

            # add rank
            addRank = rank(characterL, ranksMapL[characterL])
            rankL.append(addRank)

            # update dict, increment rank
            value = ranksMapL[characterL]
            ranksMapL[characterL] = (value + 1)

    F = rankF
    L = rankL

    return (F, L)

# getInterval for align()
# seq here is just from seedStart to seedEnd
# unbelieveably confusing and sloppy
def getInterval(fmIndex, seq):

    # given some fmIndex, lets extract relevant info for our backwards search!

    # our F and L colums in object rank() form
    F, L = FLMap(fmIndex.firstBWM, fmIndex.lastBWM)
    
    # our occTable for the F and L columns
    occF = fmIndex.occTable[0]
    occL = fmIndex.occTable[1]

    # a tuple of some interval [a, b)
    interval = (1000000, -1)

    # integer representing length of our match
    matchLength = 0

    # is this our first time thru loop?
    first = True

    # which column are we checking?
    check = ""

    # list of L-column ranks we need to find in F
    ranksToFind = []

    done = False

    # LF-Mapping
    # so let's loop backwards through our sequence, for our backward search
    # along the way let's count match length and keep track of interval
    for c in seq[::-1]:

        """
        print debugging purposes
        if first:
            print("Begin")
        print("curr char: " + c)
        print("Interval: [" + str(interval[0]) + ", " + str(interval[1]) + ")")
        for x in ranksToFind:
               print(x.character + ": " + str(x.rank))
        print("\n")
        """

        # if this is our first time thru then we need to find our interval
        if first:
            intervalFStart = 0
            intervalFEnd = 0
            for i in range(len(F)):
                curr_rank= F[i]
                if curr_rank.character == c:
                    intervalFStart = i

                    # set interval to [a, b)
                    for j in range(i, len(F)):
                        curr_rankJ = F[j]

                        if j == len(F) - 1 and curr_rankJ.character == c:
                            intervalFEnd = j + 1
                            done = True

                            # we have our interval
                            break

                        if not curr_rankJ.character == c:
                            intervalFEnd = j
                            done = True

                            # we have our interval
                            break
                if done:
                    break        

            interval = (intervalFStart, intervalFEnd)
            first = False

        # let's check L column
        else:

            # loop L-columns [a, b)
            for i in range(interval[0], interval[1]):
                curr_rank = L[i]
                if curr_rank.character == c:
                    ranksToFind.append(curr_rank)

            # we found our interval!
            if len(ranksToFind) == 0:
                break

            # now ranksToFind holds all ranks we need to find in F
            for i in range(len(F)):
                curr_rank = F[i]
                if curr_rank.__eq__(ranksToFind[0]):
                    interval = (i, i + len(ranksToFind))
                    ranksToFind = []
                    break      

        matchLength += 1


    return (interval, matchLength)

def refPos(fmIndex, interval, seedEnd, matchLength):

    # get suffixArray
    suffixArray = fmIndex.suffixArray

    # list of offsets for each match, so we know where in the 
    ref_pos = []
    
    # walk down our interval, and return the offset
    # found in the fmIndex's suffix array
    for i in range(interval[0], interval[1]):
        ref_pos.append(suffixArray[i])

    # return all ref_pos (offsets) in an array
    return ref_pos

# semi-global alignment
# seq: entire sequence
# fmIndex: single fm
# ref_pos: an offset in the ENTIRE sequence where we a match to our seed
# gap: padding for errors, idk?
def fittingAlignment(seq, fmIndex, ref_pos, gap):

    # ref sequence
    ref = fmIndex.seq

    # length of each sequence
    m = len(seq) + 1 # i
    n = fmIndex.length + 1 # j

    # build empty edit-distance matrix
    # OPT(n, m) = OPT(fm.length, len(seq))
    OPT = [[0 for x in range(n)] for y in range(m)]

    # match = 0, mismtach cost = 2, gap = 2 (maybe should be 5?)
    # for each OPT(i, j) compare subtring of seq and ref, find edit distance
    # and place edit distance into that (i, j)

    gap = 3
    mismatch = 1

    # let's fill in the the first column of OPT
    for i in range(m):
        OPT[i][0] = gap*i

    # let's fill in the first row of OPT
    for j in range(n):
        OPT[0][j] = gap*j

    # find editDistance and fill OPT
    for i in range(1, m- 1):
        for j in range(1, n- 1):
            if seq[i - 1] == ref[i - 1]:
                OPT[i][j] = OPT[i - 1][j - 1]
            else:
                OPT[i][j] = min(
                    OPT[i - 1][j - 1] + mismatch, 
                    OPT[i - 1][j] + gap, 
                    OPT[i][j - 1] + gap
                    )
    
    # construct our soltuion align
    maxLength = n + m
    i = m - 1
    j = n - 1

    seqPos = maxLength - 1
    refPos = maxLength - 1

    seqRes = [0] * (maxLength)
    refRes = [0] * (maxLength)

    while (not i == 0) or (not j == 0):

        #print("i: " + str(i - 1))
        #print("j: " + str(j - 1))

        if seq[i - 1] == ref[j - 1]:
            seqRes[seqPos] = seq[i - 1]
            refRes[refPos] = ref[j - 1]
            seqPos -= 1
            refPos -= 1
            i -= 1
            j -= 1
        elif ((OPT[i - 1][j - 1] + mismatch) == (OPT[i][j])):
            seqRes[seqPos] = seq[i - 1]
            refRes[refPos] = ref[j - 1]
            seqPos -= 1
            refPos -= 1
            i -= 1
            j -= 1
        elif ((OPT[i - 1][j] + gap) == (OPT[i][j])):
            seqRes[seqPos] = seq[i - 1]
            refRes[refPos] = '-'
            seqPos -= 1
            refPos -= 1
            i -= 1
        elif ((OPT[i][j - 1] + gap) == (OPT[i][j])):
            seqRes[seqPos] = '-'
            refRes[refPos] = ref[j - 1]
            seqPos -= 1
            refPos -= 1
            j -= 1
    
    while seqPos > 0:
        if i > 0:
            i -= 1
            seqRes[seqPos] = seq[i]
            seqPos -= 1
        else:
            seqRes[seqPos] = '-'
            seqPos -= 1

    while refPos > 0:
        if j > 0:
            j -= 1
            refRes[refPos] = ref[j]
            refPos -= 1
        else:
            refRes[refPos] = '-'
            refPos -= 1

    # removing '-' characters
    id = 1
    for i in range(maxLength - 1, 0, -1):
        if ((refRes[i] == '-') and (seqRes[i] == '-')):
            id = i + 1
            break
    
    # removing '-' characters
    align = ""
    alignPos = 0
    for i in range(1, len(seqRes)):
        if not seqRes[i] == '-':
            alignPos = i
            break
    
    for i in range(alignPos, len(seqRes)):
        align += seqRes[i]
    
    # return resulting alignment object
    # return alignment(score, align, pos)
    return alignment(OPT[m-1][n-1], align, ref_pos) # OPT, align, refRes

def writeToSam(index, output, out_sam):

    alignment = index.align

    # create sam
    QNAME = index.name
    FLAG = 0
    RNAME = 'MN988713.1'
    POS = alignment.score
    MAPQ = 255
    CIGAR = '100M'
    RNEXT = '*'
    PNEXT = 0
    TLEN = 0
    SEQ = index.seq
    QUAL = '*'

    sam = Sam(QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, [])

    # take index, write to sam file
    outFile = open(output, "a")
    out_sam = Writer(outFile, None)
    out_sam.write(sam)
    out_sam.close()



    return None   