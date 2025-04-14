from flask import Flask, render_template, request, redirect, url_for
from Bio import SeqIO
import os
import matplotlib.pyplot as plt
from models import SequenceAlignment

app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
STATIC_FOLDER = 'static'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# Ensure necessary folders exist
if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)
if not os.path.exists(STATIC_FOLDER):
    os.makedirs(STATIC_FOLDER)

# Global variables to store data
records_cache = []  # Store parsed records temporarily
all_stats = []  # Store cumulative statistics for all uploaded files
file_labels = []  # Store file labels corresponding to sequences
all_motif_results = []  # Store all motif search results

# Utility functions
def parse_fasta(filepath):
    """
    Parse a FASTA file and return a dictionary of sequences.
    :param filepath: Path to the FASTA file.
    :return: Dictionary with sequence IDs as keys and descriptions/sequences as values.
    """
    records = SeqIO.parse(filepath, "fasta")
    parsed_data = {}
    for rec in records:
        parsed_data[rec.id] = {
            "description": rec.description,
            "sequence": str(rec.seq)
        }
    return parsed_data

def get_stats(parsed_data):
    """
    Calculate statistics for each sequence in the parsed FASTA data.
    :param parsed_data: Dictionary of parsed FASTA data.
    :return: List of dictionaries containing statistics for each sequence.
    """
    stats = []
    for seq_id, data in parsed_data.items():
        seq = data["sequence"]
        gc = 100 * (seq.count("G") + seq.count("C")) / len(seq)
        stats.append({
            "id": seq_id,
            "description": data["description"],
            "length": len(seq),
            "gc_content": round(gc, 2)
        })
    return stats

def find_motifs(parsed_data, motif):
    """
    Search for a motif in each sequence and return results.
    :param parsed_data: Dictionary of parsed FASTA data.
    :param motif: Motif to search for.
    :return: List of dictionaries containing motif search results for each sequence.
    """
    results = []
    for seq_id, data in parsed_data.items():
        seq = data["sequence"]
        positions = [i for i in range(len(seq)) if seq.startswith(motif, i)]
        results.append({
            "id": seq_id,
            "description": data["description"],
            "motif": motif,
            "count": len(positions),
            "positions": positions
        })
    return results

def align_sequences(parsed_data, id1, id2, gap_pen=-2, match=1, mismatch=-1, algo="global"):
    """
    Align two sequences by their IDs using the SequenceAlignment class.
    :param parsed_data: Dictionary of parsed FASTA data.
    :param id1: ID of the first sequence.
    :param id2: ID of the second sequence.
    :param gap_pen: Gap penalty for alignment.
    :param match: Match score for alignment.
    :param mismatch: Mismatch penalty for alignment.
    :param algo: Alignment algorithm ("global" or "local").
    :return: Formatted alignment result as a string.
    """
    seq1 = parsed_data[id1]["sequence"]
    seq2 = parsed_data[id2]["sequence"]
    alignment = SequenceAlignment(seq1, seq2)
    seq1_gapped, comparison, seq2_gapped = alignment.align_sequences(gap_pen, match, mismatch, algo)
    score = alignment.get_alignment_scores(gap_pen, match, mismatch, algo)
    formatted_alignment = (
        f"{seq1_gapped}\n"
        f"{comparison}\n"
        f"{seq2_gapped}\n"
        f"Alignment Score: {score}"
    )
    return formatted_alignment

def plot_gc(stats):
    """
    Generate a GC content histogram for all sequences with repeating colors after all colors are used.
    :param stats: List of dictionaries containing sequence statistics.
    :return: Path to the saved plot.
    """
    ids = [entry["id"] for entry in stats]
    gc = [entry["gc_content"] for entry in stats]

    # Generate a colormap with a fixed number of distinct colors (e.g., 20 colors)
    num_colors = 20  # Number of distinct colors in the colormap
    cmap = plt.cm.get_cmap('tab20', num_colors)  # Use 'tab20' colormap with 20 colors
    colors = [cmap(i % num_colors) for i in range(len(ids))]  # Cycle through colors using modular arithmetic

    plt.figure(figsize=(12, 6))
    plt.bar(ids, gc, color=colors, width=0.5)
    plt.xlabel("Sequence ID")
    plt.ylabel("GC Content (%)")
    plt.title("GC Content of Sequences")
    plt.ylim(0, 100)  # Set y-axis to percentage scale
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    plot_path = f'static/gc_histogram.png'
    plt.savefig(plot_path)
    plt.close()
    return plot_path

@app.route('/')
def home():
    """Home page with buttons for different functionalities."""
    return render_template('home.html')

@app.route('/statistics', methods=['GET', 'POST'])
def statistics():
    """Page for viewing FASTA sequence statistics."""
    global records_cache, all_stats
    stats = []
    error_message = None
    histogram_path = None

    if request.method == 'POST':
        file = request.files.get('fasta_file')
        if file:
            if not file.filename.endswith(('.fasta', '.fa')):
                error_message = "Invalid file format. Please upload a FASTA file."
            else:
                filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
                file.save(filepath)
                records_cache = parse_fasta(filepath)
                stats = get_stats(records_cache)
                all_stats = stats
                histogram_path = plot_gc(stats)

    return render_template('statistics.html', stats=all_stats, error_message=error_message, histogram_path=histogram_path)

@app.route('/alignment', methods=['GET', 'POST'])
def alignment():
    """Page for aligning two sequences by their IDs."""
    global records_cache
    alignment_result = None
    error_message = None

    if request.method == 'POST':
        file = request.files.get('fasta_file')
        id1 = request.form.get('id1')
        id2 = request.form.get('id2')
        gap_pen = int(request.form.get('gap_pen', -2))
        match = int(request.form.get('match', 1))
        mismatch = int(request.form.get('mismatch', -1))
        algo = request.form.get('algo', 'global')

        if file:
            if not file.filename.endswith(('.fasta', '.fa')):
                error_message = "Invalid file format. Please upload a FASTA file."
            else:
                filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
                file.save(filepath)
                records_cache = parse_fasta(filepath)

        if id1 and id2 and records_cache:
            if id1 in records_cache and id2 in records_cache:
                alignment_result = align_sequences(records_cache, id1, id2, gap_pen, match, mismatch, algo)
            else:
                error_message = "One or both sequence IDs not found in the uploaded file."

    return render_template('alignment.html', alignment=alignment_result, error_message=error_message)

@app.route('/motif_search', methods=['GET', 'POST'])
def motif_search():
    """Page for searching motifs in a FASTA file."""
    global records_cache
    motif_results = []  # Clear previous results
    error_message = None

    if request.method == 'POST':
        motif = request.form.get('motif')
        if motif and records_cache:
            motif_results = find_motifs(records_cache, motif)

    return render_template('motif_search.html', motif_results=motif_results, error_message=error_message)

if __name__ == '__main__':
    app.run(debug=True)