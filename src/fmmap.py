#! /usr/bin/env python3

"""
TO DO:
    align()
"""

# handles SAM files/formatting
from simplesam import Reader, Writer

# using math lib
import math

# read gzip file
import gzip

# object structure for our FM-Indexes
from fmIndex import fmIndex

# read from command-line
import sys

# reads FASTA sequences
from readfq import readfq

# import all helper functions
from helperFunctions import *

# to read file extensions
import os

# object structure for alignment
from alignment import alignment

# our fm-index mapper tool
class fmmap:

    def __init__(self, seq, output):

        # our FASTA file containing our sequence
        self.seq = seq

        # our output file
        self.output = output
    
    """
    *index command*
    
        - Execute: ./fmmap.py index input output pickle
        - index: command
        - input: some sequence in .txt, .fa, .gz
        - output: some output file, preferrably .txt
        - pickle: boolean, if we're ultimately going to stuff this index into align then write true
    """
    def index(self, shouldWePickle):

        # our output file we will write to
        output = self.output

        # determine type of file extension
        fileType = os.path.splitext(self.seq)[1]

        # open based on whether its .txt, .fa, .gz
        if fileType == ".gz":

            # open our FASTA file sequence we will be reading (gzip)
            seq = gzip.open(self.seq, "r") # might need to be "rt"

        else:

            # open our FASTA file sequence we will be reading
            seq = open(self.seq, "r") # was "r"

        # array of seq names (i'd just like to keep track of them)
        seqNames = []

        # Array of fmIndex
            # each consisting of: [name, seq, seq-length, suffixArray, firstBWM, lastBWM, occTable]
        fmIndexes = []

        # using readfq to read sequence
        for name, sequenceRead, qual in readfq(seq):

            # replace any "N" reads with "A"
            # I just chose randomly
            # but maybe this is done in aligning?
            sequence = sequenceRead.replace("N", "A")

            # populate fmIndexes[]
            if (not name in seqNames):

                # add to list of seq names used
                seqNames.append(name)

                # build suffix array
                suffixArray = buildSuffixArray(sequence, len(sequence), False)

                # build firstBWM
                firstBWM = BWM(sequence, len(sequence), "first")

                # build lastBWM
                lastBWM = BWM(sequence, len(sequence), "last")

                # build occTable
                occT = occTable(sequence, len(sequence))

                # build our fmIndex
                currentFMIndex = fmIndex(name, sequence, len(sequence), suffixArray, firstBWM, lastBWM, occT, None)

                # add it to fmIndexes
                fmIndexes.append(currentFMIndex)



            # we should never reach here, but if we do we know simulated names can be duplicates
            else:

                # prints in red for fun
                print('\033[31m' + "ERROR, SEQUENCE NAME (" + name + ") WAS ALREADY USED")

                # reset to default color
                print('\033[39m') 

                # exit program
                sys.exit(1)

        if shouldWePickle.lower() == "false":

            # write to output file
            # this is if you'd like to see a cleaner serialized verison
            writeToOutput(fmIndexes, output)

            return fmIndexes
        else:

            # write to output file in binary, for efficiency
            pickleDump(fmIndexes, output)

            # i want to return the array of all fmIndexes in case some other function needs to use it
            return fmIndexes


    # align command
    # ./fmmap.py align indexBuilt readSequences output
    # indexBuilt: some .txt file containing the fm-indexes (pickled)
    # readSequences: some file containing a sequence
    # output: file we're writing to
    # using seed-and-extend paradigm
    def align(self, indexBuilt):

        # a list of all our indexes built, each being it's own fmIndex object (fmIndex.py)
        refIndex = pickleLoad(indexBuilt)

        # determine type of file extension
        fileType = os.path.splitext(self.seq)[1]

        # open based on whether its .txt, .fa, .gz
        if fileType == ".gz":

            # open our FASTA file sequence we will be reading (gzip)
            seq = gzip.open(self.seq, "r") # might need to be "rt"

        else:

            # open our FASTA file sequence we will be reading (fa or txt)
            seq = open(self.seq, "r") # was "r"

        
        ninf = float("-inf")
        seedSkip = lambda x: math.floor(x / 5.0)
        gap = 5
        for name, sequenceRead, qual in readfq(seq):

            # replace any "N" reads with "A", I just chose randomly
            # but maybe this is done in aligning?
            sequence = sequenceRead.replace("N", "A")

            # to help with aligning
            alignments = []
            seqLength = len(sequence)
            bestScore = ninf
            seedPos = 0
            skip = seedSkip(seqLength)

            for seedStart in range(0, seqLength, skip):
                seedEnd = min(seqLength, seedStart + skip)

                # for each index in our reference index?
                for bwt_index in refIndex:


	                # get_interval takes a string and performs backward search until 
	                # either (1) the entire string is matched or (2) the search interval
	                # becomes empty.  The second return value, matchLength, is the length of
	                # the query matched in backward search.  If the interval is non-empty
	                # then this is just equal to `skip` above.
                    interval, matchLength = getInterval(bwt_index, sequence[seedStart:seedEnd])


                    # given all the places where the seed matches, look for an alignment around it
	                # the ref_positions member of `bwt_index` will return positions on the reference 
	                # string corresponding to the *beginning of the read*, assuming there are no gaps
	                # in the alignment before the start of the seed (handling that is why we do fitting 
	                # alignment below).
                    for ref_pos in refPos(bwt_index, interval, seedEnd, matchLength):

                        # Perform a "fitting" alignment of the query (sequence) into the reference (bwt_index)
	                    # the returned alignment object contains the score, the alignment and the 
	                    # implied position where the query (sequence) begins under the alignment.

	                    # To perform the fitting alignment, you should "slice out" a region of the 
	                    # reference around the implied start position (ref_pos) with a bit of padding
	                    # (e.g. gap bases) before the first base where the read would start and after
	                    # the last base where the read would end.  This will ensure that the fitting_alignment
	                    # procedure can find the optimal position for the query within the reference
	                    # window that contains it, even if there are insertions or deletions in the read.
                        alignment = fittingAlignment(sequence, bwt_index, ref_pos, gap)

                        # set alignment in bwt_index
                        bwt_index.align = alignment

                        if alignment.score > bestScore:
                            bestScore = alignment.score
                            alignments = [bwt_index] #[alignment]
                        elif alignment.score == bestScore:
                            alignments.append(bwt_index) #append(alignment)
            
            # should I pass each alignment, our just the entire index
            # this empties the sam file
            f = open(self.output, "w")
            f.write("")
            f.close
            for bwt_index in alignments:
                writeToSam(bwt_index, self.output)
        

# our main driver
def main():
    
    # determine if a command was entered
    if (not len(sys.argv) > 1):
        print("please enter a valid command: index or align")
    else:

        # read from command-line
        command = sys.argv[1]

        # determine which command to use

        # index
        if (command == "index"):

            # path to fasta, output file, and boolean
            seq = sys.argv[2]
            output = sys.argv[3]
            shouldWePickle = sys.argv[4]

            # create fmmap with a sequence and output file
            FM = fmmap(seq, output)

            # call index command
            FM.index(shouldWePickle)

        # align
        elif (command == "align"):

            # .txt file containing pickle, path to fasta, and output file
            indexBuilt = sys.argv[2]
            seq = sys.argv[3]
            output = sys.argv[4]

            # create fmmap with a sequence and output file
            FM = fmmap(seq, output)
            FM.align(indexBuilt)

        # no valid commands
        else:
            print(command + " is not a valid command.")

# function-call on main()
main()