# Mitochondrial DNA Analysis Tool

A Python-based tool for analyzing mitochondrial DNA sequences, implementing sequence alignment, motif searching, and statistical analysis.

## Table of Contents
1. [Project Overview](#project-overview)
2. [System Design](#system-design)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Examples](#examples)
6. [Project Structure](#project-structure)
7. [Design Documentation](#design-documentation)

## Project Overview

This tool provides functionality for:
- FASTA file parsing and data management
- Mitochondrial DNA sequence analysis
- Motif searching and pattern matching
- Global and local sequence alignment
- GC content analysis and visualization

The project integrates a web interface with a command line analysis script, providing both interactive and batch processing capabilities.

## System Design

### CRC Cards

#### MitochondrialDNA
- **Class**: MitochondrialDNA
- **Responsibilities**:
  - Store DNA sequence data
  - Calculate sequence properties
  - Extract subsequences
- **Collaborators**: Parser, SequenceAlignment

#### MotifFinder
- **Class**: MotifFinder
- **Responsibilities**:
  - Search for patterns in sequences
  - Count motif occurrences
- **Collaborators**: FMIndexQuery

#### SequenceAlignment
- **Class**: SequenceAlignment
- **Responsibilities**:
  - Perform sequence alignments
  - Calculate alignment scores
- **Collaborators**: GlobalAlignment, LocalAlignment

#### WebApp
- **Class**: WebApp
- **Responsibilities**:
  - Handle user interactions
  - Integrate command line analysis
- **Collaborators**: MitochondrialDNA, SequenceAlignment, Part3

### UML Diagram

```mermaid
classDiagram
    class MitochondrialDNA {
        +str seq
        +str ID
        +str description
        +get_subsequence(start, end)
        +get_GC_content()
        +get_length()
    }
    
    class MotifFinder {
        +str motif_seq
        +count_occurrences(target_seq)
        +search_motif(target_seq)
    }
    
    class SequenceAlignment {
        +str seq1
        +str seq2
        +align_sequences()
        +get_alignment_scores()
    }
    
    class WebApp {
        +app.py
        +handle_uploads()
        +display_results()
        +integrate_analysis()
    }
    
    class AnalysisScript {
        +Part 3.py
        +run_analysis()
        +generate_plots()
    }
    
    MitochondrialDNA --> SequenceAlignment
    MotifFinder --> FMIndexQuery
    SequenceAlignment --> GlobalAlignment
    SequenceAlignment --> LocalAlignment
    WebApp --> MitochondrialDNA
    WebApp --> SequenceAlignment
    WebApp --> AnalysisScript
```

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/mitochondrial-dna-analysis.git
cd mitochondrial-dna-analysis
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Web Interface
1. Start the Flask server:
```bash
python app.py
```

2. Open your browser to `http://localhost:5000`

3. Upload a FASTA file and use the interface to:
   - View sequence statistics
   - Perform sequence alignments
   - Search for motifs
   - View visualizations
   - Run integrated analysis from `Part 3.py`

### Command Line Analysis
Run the analysis script directly:
```bash
python "Part 3"
```

## Examples

### Input Example (FASTA format)
```
>NC_012920.1 Homo sapiens mitochondrion
GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGG
```

### Output Examples

1. Sequence Statistics:
```
Sequence ID:
Length: 
GC Content: 
```

2. Motif Search:
```
Motif: 
Occurrences: 
Positions: 
```

3. Alignment Results:
```
Sequence 1: 
Sequence 2: 
Comparison: 
Score: 
```

## Project Structure

```
.
├── app.py                 # Flask web application
├── Part 3.py             # Command line analysis script
├── models.py             # Core DNA analysis classes
├── parser.py             # FASTA file parsing
├── global_alignment_algo.py  # Global sequence alignment
├── local_alignment_algo.py   # Local sequence alignment
├── fm_index_query.py     # Pattern matching
├── templates/            # HTML templates
│   ├── home.html
│   ├── statistics.html
│   ├── alignment.html
│   └── motif_search.html
├── static/              # Static files
└── uploads/             # Uploaded files
```

## Design Documentation

### Genomic Element Modeling

1. **MitochondrialDNA Class**
   - Represents individual mitochondrial DNA sequences
   - Encapsulates sequence data and analysis methods
   - Provides methods for sequence manipulation and analysis

2. **MotifFinder Class**
   - Implements pattern searching using FM-Index
   - Provides efficient motif searching capabilities
   - Returns both counts and positions of motifs

3. **SequenceAlignment Class**
   - Implements both global and local alignment algorithms
   - Provides configurable alignment parameters
   - Returns alignment results and scores

### Component Interactions

1. **Data Flow**:
   - FASTA files → Parser → MitochondrialDNA objects
   - MitochondrialDNA → SequenceAlignment for comparisons
   - MitochondrialDNA → MotifFinder for pattern searching
   - WebApp → AnalysisScript for integrated analysis

2. **Analysis Pipeline**:
   - Sequence loading and validation
   - Statistical analysis (GC content, length)
   - Motif searching
   - Sequence alignment
   - Results visualization

### Object-Oriented Principles

1. **Encapsulation**:
   - Each class encapsulates its data and methods
   - Internal implementation details are hidden
   - Clear public interfaces

2. **Abstraction**:
   - High-level sequence manipulation
   - Simplified analysis interfaces
   - Clear separation of concerns

3. **Modularity**:
   - Independent components
   - Clear interfaces
   - Easy to extend

4. **Extensibility**:
   - New analysis methods can be added
   - Support for additional file formats
   - Custom visualization options
