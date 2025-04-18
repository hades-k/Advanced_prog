from flask import Flask, render_template, request
from part3 import load_genomes, align_genomes, find_motifs, compare_to_reference, visualize_differences_bar
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
STATIC_FOLDER = 'static'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['STATIC_FOLDER'] = STATIC_FOLDER

# Ensure necessary folders exist
if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)
if not os.path.exists(STATIC_FOLDER):
    os.makedirs(STATIC_FOLDER)

# Global variable to store genomes
genomes = []

@app.route('/')
def home():
    """Home page with navigation."""
    return render_template('home.html')

@app.route('/compare', methods=['GET', 'POST'])
def compare_genomes():
    """Compare two genomes."""
    global genomes
    comparison_result = None
    summary = None
    error_message = None

    if request.method == 'POST':
        # Handle file upload
        file = request.files.get('fasta_file')
        if file and file.filename.endswith(('.fasta', '.fa')):
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
            file.save(filepath)
            genomes = load_genomes(filepath)  # Use the function from part3.py
        else:
            error_message = "Invalid file format. Please upload a FASTA file."

        # Perform genome comparison if IDs are provided
        id1 = request.form.get('id1')
        id2 = request.form.get('id2')
        if id1 and id2 and genomes:
            g1 = next((g for g in genomes if g.ID == id1), None)
            g2 = next((g for g in genomes if g.ID == id2), None)

            if g1 and g2:
                # Align genomes and generate comparison result
                aligned_seq1, comparison_line, aligned_seq2 = align_genomes(g1, g2)  # Use the function from part3.py
                comparison_result = [
                    {"block1": aligned_seq1, "comp_line": comparison_line, "block2": aligned_seq2}
                ]
                # Generate summary statistics
                matches = comparison_line.count('|')
                mismatches = comparison_line.count('*')
                gaps = aligned_seq1.count('-') + aligned_seq2.count('-')
                total = len(comparison_line)
                summary = {
                    "matches": matches,
                    "mismatches": mismatches,
                    "gaps": gaps,
                    "total": total
                }

    return render_template('compare.html', genomes=genomes, comparison_result=comparison_result, summary=summary, error_message=error_message)

@app.route('/motif_search', methods=['GET', 'POST'])
def motif_search():
    """Find conserved motifs and specific motif positions."""
    global genomes
    motifs = []  # Initialize motifs as an empty list
    conservation_matrix = None
    motif_results = []
    specific_motif = None  # Initialize specific_motif to avoid UnboundLocalError
    error_message = None

    if request.method == 'POST':
        # Handle file upload
        file = request.files.get('fasta_file')
        if file and file.filename.endswith(('.fasta', '.fa')):
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
            file.save(filepath)
            genomes = load_genomes(filepath)  # Use the function from part3.py
        else:
            error_message = "Invalid file format. Please upload a FASTA file."

        # Perform specific motif search
        specific_motif = request.form.get('specific_motif')  # User provides a single motif
        if specific_motif and genomes:
            motif_results = find_motifs(genomes, specific_motif)  # Use the function from part3.py

    return render_template(
        'motif_search.html',
        genomes=genomes,
        motifs=motifs,
        conservation_matrix=conservation_matrix,
        specific_motif=specific_motif,
        motif_results=motif_results,
        error_message=error_message
    )

@app.route('/reference', methods=['GET', 'POST'])
def reference_genome():
    """Compare genomes to a reference genome."""
    global genomes
    results = None
    reference_id = None
    error_message = None

    if request.method == 'POST':
        # Handle file upload
        if 'fasta_file' in request.files:
            file = request.files.get('fasta_file')
            if file and file.filename.endswith(('.fasta', '.fa')):
                filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
                file.save(filepath)
                genomes = load_genomes(filepath)  # Use the function from part3.py
            else:
                error_message = "Invalid file format. Please upload a FASTA file."

        # Handle reference genome selection
        reference_id = request.form.get('reference_id')
        if reference_id and genomes:
            reference = next((g for g in genomes if g.ID == reference_id), None)
            if reference:
                results = compare_to_reference(reference, [g for g in genomes if g.ID != reference_id])  # Use the function from part3.py
            else:
                error_message = f"Reference genome with ID '{reference_id}' not found."

    return render_template('reference.html', genomes=genomes, results=results, reference_id=reference_id, error_message=error_message)

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

    plot_path = os.path.join(app.config['STATIC_FOLDER'], 'gc_histogram.png')
    plt.savefig(plot_path)
    plt.close()
    return plot_path

@app.route('/statistics', methods=['GET', 'POST'])
def genome_statistics():
    """Page for viewing FASTA sequence statistics."""
    global genomes
    stats = []
    error_message = None
    gc_bar_chart_path = None
    gc_pie_chart_path = None
    gc_species_histogram_path = None

    if request.method == 'POST':
        # Handle file upload
        file = request.files.get('fasta_file')
        if file and file.filename.endswith(('.fasta', '.fa')):
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], file.filename)
            file.save(filepath)
            genomes = load_genomes(filepath)  # Use the function from part3.py

            # Calculate statistics
            stats = [
                {
                    "id": genome.ID,
                    "description": genome.description,
                    "length": len(genome.seq),
                    "gc_content": round((genome.seq.count('G') + genome.seq.count('C')) / len(genome.seq) * 100, 2)
                }
                for genome in genomes
            ]

            # Generate GC content bar chart
            gc_bar_chart_path = plot_gc(stats)

            # Generate overall GC content pie chart
            total_bases = sum(stat["length"] for stat in stats)
            total_gc_bases = sum((stat["gc_content"] / 100) * stat["length"] for stat in stats)
            total_at_bases = total_bases - total_gc_bases
            plt.figure(figsize=(8, 8))
            plt.pie(
                [total_gc_bases, total_at_bases],
                labels=["GC Content", "AT Content"],
                autopct='%1.1f%%',
                colors=["orange", "skyblue"],
                startangle=140
            )
            plt.title("Overall GC Content in File")
            gc_pie_chart_path = os.path.join(app.config['STATIC_FOLDER'], 'gc_pie_chart.png')
            plt.savefig(gc_pie_chart_path)
            plt.close()

            # Generate GC content distribution histogram
            gc_contents = [stat["gc_content"] for stat in stats]
            plt.figure(figsize=(10, 6))
            plt.hist(gc_contents, bins=10, color='teal', edgecolor='black')
            plt.title("GC Content Distribution Across Species")
            plt.xlabel("GC Content (%)")
            plt.ylabel("Frequency")
            gc_species_histogram_path = os.path.join(app.config['STATIC_FOLDER'], 'gc_species_histogram.png')
            plt.savefig(gc_species_histogram_path)
            plt.close()

        else:
            error_message = "Invalid file format. Please upload a FASTA file."

    return render_template(
        'statistics.html',
        stats=stats,
        error_message=error_message,
        gc_bar_chart_path='/static/gc_histogram.png',
        gc_pie_chart_path='/static/gc_pie_chart.png',
        gc_species_histogram_path='/static/gc_species_histogram.png'
    )

if __name__ == '__main__':
    app.run(debug=True)