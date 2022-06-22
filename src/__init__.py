from Bio import SeqIO
from Bio.SeqUtils import GC as gc
from Bio import motifs
from GibbsFE import GibbsFE
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
hist_pos = []
hist_neg = []
bin_edges = 0
with open('settings.json') as settings:
    json_file = json.load(settings)
def ener_sets():
    pos_ener = []
    neg_ener = []
    wsize = json_file["wsize"]
    entries_prom = Path(json_file["promoter_folder"] + '/')
    #look into altering code for multiple files in promoters folder
    for entries in entries_prom.iterdir():
        seq_objectp = SeqIO.parse(entries,'fasta')
    #open reference set to estimate energy distribution
    #process all sequences in reference set
    #list of average energies
    for aseq in seq_objectp:   
        es = GibbsFE(aseq, 0)   #get list of energies
        pos_ener.append(sum(es)/len(es))   #story average energy for the sequence
    entries_gen = Path(json_file["background_set"] + '/')
    for entries in entries_gen.iterdir():
        seq_objectg = SeqIO.read(entries,'fasta')
    for p in range(0,len(seq_objectg)-wsize+1,25):
        es = GibbsFE(seq_objectg[p:p+wsize], 0)
        neg_ener.append(sum(es)/len(es))   #store average energy for the window
    # Frequency distributions (on positive and negative set)
    # Defining the frequency distributions for the positive and the
    # negative sets requires binning.
    return(hist(pos_ener,neg_ener))

def hist(pos_ener,neg_ener):
    # Set number of bins
    n_bins = 30
    
    # Set pseudocount values
    pc = 1
    
    # Energy frequencies for the negative set
    hist_neg, bin_edges = np.histogram(neg_ener, range=(2.43 , 9.37), density=False, bins=n_bins)
    hist_neg = (hist_neg+pc)/(hist_neg+pc).sum()
    
    # Energy frequencies for the positive set
    hist_pos, bin_edges = np.histogram(pos_ener, range=(2.43 , 9.37), density=False, bins=n_bins)
    hist_pos = (hist_pos+pc)/(hist_pos+pc).sum()
    return hist_neg, hist_pos, bin_edges
    
    # Plot frequency distributions as histograms
    plt.bar(range(n_bins), hist_neg, alpha=0.5)
    plt.bar(range(n_bins), hist_pos, alpha=0.5)  

"""Predicts putative promoters in a given sequence, following the PromoterHunter approach.
   Inputs:
   - mySeq - DNA sequence [Seq, not SeqRecord object]
   - lmot - -35 motif object
   - lthrs - score threshold for -35 motif
   - rmot - -10 motif object
   - rthrs - score threshold for -10 motif
   - [minD, maxD] - range of spacer lengths
   - wsize - size of moving average window for Gibbs free energy
   - [lerg,rerg] - range surrounding -10 on which to compute energy score
   
   Returns:
   - List of High-scoring Motif Pairs (HSMP), with all the score and sequence information
     for the constituent motifs and energy
"""    
def pHunt(mySeq, lmot, lthrs, rmot, rthrs, minD, maxD, wsize, posdist, negdist, binedges):
    lerg=-wsize//2   #left window from -10 for energy
    rerg=+wsize//2    #right window from -10 for energy
    GC = gc(mySeq)/100
    #score sequence with left motif
    lscrs=lmot.pssm.calculate(mySeq).tolist()
    #score sequence with left motif
    rscrs=rmot.pssm.calculate(mySeq).tolist()
    #compute the Gibbs Free Energy vector
    fes = GibbsFE(mySeq, None)
    
    #list of spacer-adequate high-scoring motif pairs
    hsmp = []
    
    mySeqlen = len(mySeq)
    #go through sequence, first with left PSSM    
    for pl in range(mySeqlen-(lmot.pssm.length + minD + rmot.pssm.length)):
        #if score above threshold
        if lscrs[pl] > lthrs:
            remrange = pl+lmot.pssm.length+maxD+1 if (pl+lmot.pssm.length+maxD+1 < mySeqlen - rmot.pssm.length) else mySeqlen - rmot.pssm.length
            #go through sequence, now with right PSSM, up to spacer
            for pr in range(pl+lmot.pssm.length+minD,remrange):
                #if score above threshold
                if rscrs[pr] > rthrs:
                    element = {'lpos' : pl, 'lseq' : str(mySeq[pl:pl+lmot.pssm.length]),
                               'lscr' : lscrs[pl],
                               'rpos' : pr, 'rseq' : str(mySeq[pr:pr+rmot.pssm.length]),
                               'rscr' : rscrs[pr],
                               'spcr' : str(mySeq[pl+lmot.pssm.length:pr].lower()),
                              }
                    #compute overall (normalized) PSSM contribution
                    element['mscr'] = element['lscr'] + element['rscr']

                    #add Gibbs Free Energy component
                    #as average betweeen coordinates (plus lerg, rerg margins)
                    lrange = pr+lerg if pr+lerg>0 else 0
                    rrange = pr+rerg if pr+rerg<mySeqlen-rmot.pssm.length else mySeqlen-rmot.pssm.length
                    element['escr'] = sum(fes[lrange:rrange])/(rrange-lrange)
                    
                    p_pos = posdist[np.where(binedges>element['escr'])[0][0] - 1]
                    p_neg = negdist[np.where(binedges>element['escr'])[0][0] - 1]
                    
                    element['ellr'] = np.log2(p_pos / p_neg)
                    
                    #compute final score (global norm motif score + 2*norm energy score)
                    element['Fscr'] = element['mscr'] + element['ellr']
                    
                    #add high-scoring motif pair with all info to the list of HSMPs
                    hsmp.append(element)
                  
    #return list of spacer-adequate high-scoring motif pairs
    return(hsmp)
def motif_create(paths,i):
    '''
    Creates or reads a motif from a fasta, jaspar, or xml file.
    Parameters
    ----------
    file_name : String
        file name with motif or instances (Assumes the file directory is in the motif folder within working directory)
    i : integer
        If reading a jaspar or xml file, specifies the wanted motif. 
    Returns
    -------
    motif : Motif Object
        motif created or read from the file
    '''
    entries = Path(str(paths) + '/')
    for entry in entries.iterdir():
            if '.fas' in str(entry):         
                seqlist = []
                seqs = SeqIO.parse(str(entry),'fasta')
                for seq in seqs:
                    seqlist.append(seq.seq)
                motif = motifs.create(seqlist,'ACGT')
            elif '.jaspar' in str(entry):
                with open(str(entry)) as q:
                    motif = motifs.parse(q,'jaspar')[i]
            elif '.xml' in str(entry):
                with open(str(entry)) as w:
                    motif = motifs.parse(w,'meme')[i]
            else:
                print('Unrecognized file format')
    return(motif)
def go():
    listz = []
    listz = ener_sets()
    entries_ref = Path(json_file["reference_set"])
    for entries in entries_ref.iterdir():
        mySeq = list(SeqIO.parse(entries,'fasta'))
    entries_lmot = Path(json_file["left_motif_folder"])
    entries_rmot = Path(json_file["right_motif_folder"])
    lmot = motif_create(entries_lmot,0)
    lmot._pseudocounts = 0.05
    lthrs = lmot.pssm.distribution(precision=10**3).threshold_balanced()
    rmot = motif_create(entries_rmot,0)
    rmot._psuedocounts = 0.05
    rthrs = rmot.pssm.distribution(precision=10**3).threshold_patser()
    minD = json_file["min_spacer_value"]
    maxD = json_file["max_spacer_value"]
    wsize = json_file["wsize"]
    posdist = listz[1]
    negdist = listz[0]
    binedges = listz[2]
    mypromoters = pHunt(mySeq[0].seq,lmot,lthrs,rmot,rthrs,minD,maxD,wsize, posdist, negdist, binedges)
    mypromoters = sorted(mypromoters, key = lambda k: k['Fscr'], reverse=True)
    for promoter in mypromoters:
        print('Range:', promoter['lpos'],'-',promoter['rpos'], '|', promoter['rpos']-promoter['lpos']-6, \
              '- L:',round(promoter['lscr'],2),'R:',round(promoter['rscr'],2),\
              'FE:',round(promoter['escr'],2), 'FELLR:',round(promoter['ellr'],2), 'T:',round(promoter['Fscr'],2),\
             '- S:',promoter['lseq']+promoter['spcr']+promoter['rseq'])
'''
    entries_txt = Path(json_file["output_folder"] + "/")
    f = open(str(entries_txt) + "Output_Information.txt", "w")
    f.write(('Range:' + str(promoter['lpos']) +'-' + str(promoter['rpos']) + '|' + str(promoter['rpos']-promoter['lpos']-6) + 
          '- L:' + str(round(promoter['lscr'],2)) + 'R:' + str(round(promoter['rscr'],2)) + 'FE:' + str(round(promoter['escr'],2) + 'FELLR:' + str(round(promoter['ellr'],2)) + str('T:',round(promoter['Fscr'],2)) + '- S:' + str(promoter['lseq']+promoter['spcr']+promoter['rseq']))))
    f.close()
'''

if __name__ == "__main__":
    go()
