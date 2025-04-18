from parser import *
from global_alignment_algo import *
from local_alignment_algo import *
from fm_index_query import *

class MitochondrialDNA:
    
    def __init__(self, seq:str, ID:str, description:str = ""):
        self.seq = seq
        self.ID = ID
        self.description = description
    
    def get_subsequence(self, start:int, end:int):
        if start < 0 or end > len(self.seq):
            raise ValueError(f"Subsequence indices out of range: start={start}, end={end}, length={len(self.seq)}")
        else:
            return self.seq[start:end]
        
    def get_GC_content(self):
        gc_count = self.seq.count("G") + self.seq.count("C")
        return (gc_count/len(self.seq)) * 100
    
    def get_length(self):
        return len(self.seq)

class MotifFinder:
    
    def __init__(self, motif_seq:str):
        self.motif_seq = motif_seq
        
    def count_occurrences(self, target_seq:str):
        """Count occurrences of the motif in the sequence."""
        occurrences = FMIndexQuery(target_seq, self.motif_seq)[0]
        return occurrences
    
    def search_motif(self, target_seq:str):
        """Search for motifs in a sequence."""
        search = FMIndexQuery(target_seq, self.motif_seq)[1]
        return search

class SequenceAlignment:
        
    def __init__(self, seq1:str, seq2:str):
        self.seq1 = seq1
        self.seq2 = seq2
    
    def align_sequences(self, gap_pen=-2, match=1, mismatch=-1, algo:str="global"):
        """Align two mitochondrial DNA sequences"""
        if algo == "global":
            seq1_gapped, comparison, seq2_gapped = globalAlignment(self.seq1, self.seq2, gap_pen, match, mismatch)[0]
        elif algo == "local":
            seq1_gapped, comparison, seq2_gapped = localAlignment(self.seq1, self.seq2, gap_pen, match, mismatch)[0]
        return seq1_gapped, comparison, seq2_gapped
    
    def get_alignment_scores(self, gap_pen=-2, match=1, mismatch=-1, algo:str="global"):
        """Return the alignment scores."""
        if algo == "global":
            score = globalAlignment(self.seq1, self.seq2, gap_pen, match, mismatch)[1]
        elif algo == "local":
            score = localAlignment(self.seq1, self.seq2, gap_pen, match, mismatch)[1]
        return score

# Usage example (comment out later)
'''"""
data = parser('synthetic_mtDNA_dataset.fasta', 'fasta')
genomes = data[['id', 'seq', 'description']]

NC10_DNA = MitochondrialDNA(genomes.loc[9]['seq'], genomes.loc[9]['id'], genomes.loc[9]['description'])
subseq = NC10_DNA.get_subsequence(30,70)
gc = NC10_DNA.get_GC_content()
nt = NC10_DNA.get_length()
print(f" sequence between index 30 and 69: {subseq}\n the gc content is {gc}%\n this genome is {nt} nucleotides in length\n")

GATC = MotifFinder("GATC")
hits = GATC.search_motif(NC10_DNA.seq)
n = GATC.count_occurrences(NC10_DNA.seq)
print(f" GATC motif was found at these offsets in NC_10's sequence: {hits}\n total number of hits: {n}")

NC18_DNA = MitochondrialDNA(genomes.loc[17]['seq'], genomes.loc[17]['id'], genomes.loc[17]['description'])
comp_10_18 = SequenceAlignment(NC10_DNA.seq, NC18_DNA.seq)
gaps1, comparison, gaps2 = comp_10_18.align_sequences(gap_pen=-1)
score = comp_10_18.get_alignment_scores(gap_pen=-1)
print(f" pairwise global alignment between NC_10 and NC_18 genomes:\n {gaps1}\n {comparison}\n {gaps2}\n alignment score: {score}")
'''
