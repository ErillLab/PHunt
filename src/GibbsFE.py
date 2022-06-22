from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import motifs
from Bio.SeqUtils import GC as gc
import numpy as np 
import matplotlib.pyplot as plt
import json



# Define GibbsFE function
"""
GibbsFE: computes the Gibbs Free Energy of a DNA sequence
reads in:
- DNA sequence (BioPython Seq object)
- window length (to compute Free Energy over)
Assumes the DNA sequence only contains DNA characters (A, C, G, T)
"""
def GibbsFE(mySeq, wlen=None):
    mySeq = mySeq.upper()   #upper case sequence
    #Gibbs Free Sequence components defined in settings.json (manipulatable)
    with open('settings.json') as settings:
        json_file = json.load(settings)
    LUTfe = json_file["GibbsFE_Bases"]
    #convert sequence to vector
    myVec = []
    #return empty vector if sequence is less than 2 nucleotides
    mySeqlen = len(mySeq)
    if mySeqlen<2:
        return(myVec)
    #iterate over the sequence, grabbing two elements at time
    #with a +1 step size
    for p in range(mySeqlen-1):
        dinuc = mySeq[p]+mySeq[p+1]  #create dinucleotide
        myVec.append(LUTfe[dinuc])   #look up GE value in table
    
    #perform moving average, if requested (wlen!=None)
    if (wlen):
        myVeclen = len(myVec)
        #return original vector if window size over length
        if myVeclen<wlen:
            return(myVec)
        myMAvec = []
        for p in range(myVeclen-wlen+1):
            #compute average
            myMAvec.append(sum(myVec[p:p+wlen])/wlen)
        #return moving average vector
        return(myMAvec)
    else:
        return(myVec)  

"""
    
#create a sequence record/Test Code
mydna = SeqRecord(Seq('GTCAGGTGATTGACGAAGATGTCTATCCGATTCTGTCGCTGCAATCGTGCCTCGACAAGCGTGCGGCAAAAGGCGGCGTCTCACCGCAGCAGGTGGCGCAGGCGATTGCTTTTGCGCAGGCTCGGTTAGGGTAAGAACATTTATATGTATAAATTTGAGCCTGGCTTATCGCCGGGCTTTTTTATGGCAAAAAAAAGCGGATCCTGGAGATCCGCAAAAGTTCACGTTGGCTTTAGTTATTCGAGTTGAGAAACTCTCGAAACGGGCAGTGACTTCAAGGGTTAAAAGAGGTGCCGCTCCGTTTCTGTGAGCAATTATCAGTCAGAATGCTTGATAGGGAGCGCCGTTCATTGCTATTCTACCTATCGCCATGAACTATCGTGGCGATGGAGGATGGATAATGAATATTCGTGATCTTGAGTACCTGGTGGCATTGGCTGAACACCGCCATTTTCGGCGTGCGGCAGATTCCTGCCACGTTAGCCAGCCGACGCTTAGCGG'),\
                        id='NC_000913.3:4158090-4158590',\
                        description='PoxyR_Ecoli')    

#call GibbsFE function - Come Back
energies = GibbsFE(mydna,50)

#plot results
plt.plot(energies)
"""