import os
from flask import Flask, render_template, request, jsonify
from parser import parser
from models import SequenceAlignment, MitochondrialDNA, MotifFinder

app = Flask(__name__)

UPLOAD_FOLDER = "uploads"
if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

@app.route("/")
def home():
    return render_template("index.html")

@app.route("/upload", methods=["POST"])
def upload_fasta():
    if "file" not in request.files:
        return jsonify({"error": "No file uploaded"}), 400
    
    file = request.files["file"]
    if file.filename == "":
        return jsonify({"error": "No selected file"}), 400

    file_path = os.path.join(UPLOAD_FOLDER, file.filename)
    file.save(file_path)

    try:
        # Use the existing parser
        df = parser(file_path)
        sequences = []
        for _, row in df.iterrows():
            mtdna = MitochondrialDNA(row['seq'], row['id'], row.get('description', ''))
            seq_data = {
                "id": mtdna.ID,
                "description": mtdna.description,
                "sequence": mtdna.seq,
                "length": mtdna.get_length(),
                "gc_content": mtdna.get_GC_content()
            }
            sequences.append(seq_data)
        return jsonify({"sequences": sequences})
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    finally:
        if os.path.exists(file_path):
            os.remove(file_path)

@app.route("/search_motif", methods=["POST"])
def search_motif():
    data = request.get_json()
    sequence = data.get("sequence", "")
    motif = data.get("motif", "")
    
    if not sequence or not motif:
        return jsonify({"error": "Both sequence and motif are required"}), 400

    # Use the existing MotifFinder class
    finder = MotifFinder(motif)
    positions = finder.search_motif(sequence)
    count = finder.count_occurrences(sequence)

    return jsonify({
        "positions": positions,
        "count": count
    })

@app.route("/align", methods=["POST"])
def align():
    data = request.get_json()
    seq1 = data.get("seq1")
    seq2 = data.get("seq2")
    algo = data.get("algo", "global")
    gap = int(data.get("gap", -2))
    match = int(data.get("match", 1))
    mismatch = int(data.get("mismatch", -1))

    aligner = SequenceAlignment(seq1, seq2)
    aligned_seq1, comparison, aligned_seq2 = aligner.align_sequences(
        gap_pen=gap, 
        match=match, 
        mismatch=mismatch, 
        algo=algo
    )
    score = aligner.get_alignment_scores(
        gap_pen=gap, 
        match=match, 
        mismatch=mismatch, 
        algo=algo
    )

    return jsonify({
        "seq1": aligned_seq1,
        "seq2": aligned_seq2,
        "comparison": comparison,
        "score": score
    })

if __name__ == "__main__":
    print("Go to http://127.0.0.1:5000")
    app.run(host='127.0.0.1', port=5000)
    
    
