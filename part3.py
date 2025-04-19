# PART 3 of project

from models import SequenceAlignment, MitochondrialDNA, parser, MotifFinder
import numpy as np
import matplotlib.pyplot as plt

def load_genomes(filepath):
    """Load genomes from a FASTA file."""
    data = parser(filepath, 'fasta')
    genomes = [
        MitochondrialDNA(
            seq=data.loc[i]['seq'],
            ID=data.loc[i]['id'],
            description=data.loc[i]['description']
        )
        for i in range(len(data))
    ]
    return genomes

def align_genomes(genome1, genome2):
    """Align two genomes and return the alignment result."""
    aligner = SequenceAlignment(genome1.seq, genome2.seq)
    aligned_seq1, comparison_line, aligned_seq2 = aligner.align_sequences()
    return aligned_seq1, comparison_line, aligned_seq2

def visualize_differences_bar(seq1, seq2, label1, label2):
    """Visualize differences between two aligned sequences."""
    matches = mismatches = gaps = 0
    block_size = 60
    length = len(seq1)
    result = []

    for i in range(0, length, block_size):
        block1 = seq1[i:i+block_size]
        block2 = seq2[i:i+block_size]
        comp_line = ""
        
        for a, b in zip(block1, block2):
            if a == b and a != "-":
                comp_line += "*"
                matches += 1
            elif a == '-' or b == '-':
                comp_line += " "
                gaps += 1
            else:
                comp_line += "|"
                mismatches += 1
            
        result.append({
            "block1": block1,
            "block2": block2,
            "comp_line": comp_line
        })

    summary = {
        "matches": matches,
        "mismatches": mismatches,
        "gaps": gaps,
        "total": length
    }
    return result, summary

def find_motifs(genomes, motif):
    """
    Search for a motif in each genome and return the results.
    :param genomes: List of genome objects.
    :param motif: The motif to search for.
    :return: List of dictionaries containing motif search results for each genome.
    """
    results = []
    for genome in genomes:
        motif_finder = MotifFinder(motif)
        positions = motif_finder.search_motif(genome.seq)
        results.append({
            "id": genome.ID,
            "description": genome.description,
            "motif": motif,
            "count": motif_finder.count_occurrences(genome.seq),
            "positions": positions
        })
    return results

def compare_to_reference(reference, genomes):
    """Compare all genomes to a reference genome."""
    results = []
    for target in genomes:
        aligner = SequenceAlignment(reference.seq, target.seq)
        _, comp, _ = aligner.align_sequences()
        score = aligner.get_alignment_scores()
        matches = comp.count("*")
        total = len(comp)
        similarity = (matches / total) * 100 if total > 0 else 0
        results.append({
            "id": target.ID,
            "score": score,
            "similarity": similarity
        })
    return results

if __name__ == "__main__":
    print("This script is meant to be imported, not executed directly.")
