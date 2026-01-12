import os
import math
from typing import List, Tuple, Dict
from collections import Counter

#INPUT_FILE = "./data/alignments.fasta"
#OUTPUT_DIR = "./tests/analysis"
INPUT_FILE = "/workdir/data/alignments.fasta"
OUTPUT_DIR = "/workdir/analysis"

def parse_alignments(filepath: str) -> List[List[Tuple[str, str]]]:
    """Parse FASTA file into separate alignments"""
    alignments = []
    current_alignment = []
    current_id = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            
            if line.startswith('### ALIGNMENT_'):
                if current_id:
                    current_alignment.append((current_id, ''.join(current_seq)))
                if current_alignment:
                    alignments.append(current_alignment)
                current_alignment = []
                current_id = None
                current_seq = []
            elif line.startswith('>'):
                if current_id:
                    current_alignment.append((current_id, ''.join(current_seq)))
                current_id = line[1:]
                current_seq = []
            elif line:
                current_seq.append(line)
        
        if current_id:
            current_alignment.append((current_id, ''.join(current_seq)))
        if current_alignment:
            alignments.append(current_alignment)
    
    return alignments

def calculate_entropy(column: str) -> float:
    """Calculate Shannon entropy for a column"""
    non_gaps = [c for c in column if c != '-']
    if len(non_gaps) == 0:
        return 0.0
    
    counts = Counter(non_gaps)
    total = len(non_gaps)
    
    entropy = 0.0
    for count in counts.values():
        if count > 0:
            p = count / total
            entropy -= p * math.log2(p)
    
    return entropy

def get_conservation_category(column: str) -> str:
    """Determine conservation category for a column"""
    non_gaps = [c for c in column if c != '-']
    
    if len(non_gaps) == 0:
        return "GAP_ONLY"
    
    unique = set(non_gaps)
    if len(unique) == 1:
        return "CONSERVED"
    else:
        return "VARIABLE"

def get_most_common_residue(column: str) -> str:
    """Get most common residue in a column"""
    non_gaps = [c for c in column if c != '-']
    
    if len(non_gaps) == 0:
        return "-"
    
    counts = Counter(non_gaps)
    max_count = max(counts.values())
    most_common = sorted([aa for aa, count in counts.items() if count == max_count])
    
    return most_common[0]

def calculate_pairwise_identity(seq1: str, seq2: str) -> float:
    """Calculate percent identity between two sequences"""
    # CRITICAL FIX: Handle diagonal case for all-gap sequences
    if seq1 == seq2:
        return 100.0
    
    matching = 0
    comparable = 0
    
    for i in range(len(seq1)):
        a, b = seq1[i], seq2[i]
        
        # Skip gap-gap positions
        if a == '-' and b == '-':
            continue
        
        comparable += 1
        
        # Count matching non-gap positions
        if a == b and a != '-':
            matching += 1
    
    if comparable == 0:
        return 0.0
    
    return (matching / comparable) * 100.0

def process_alignment(sequences: List[Tuple[str, str]], aln_num: int):
    """Process a single alignment and generate output files"""
    seq_ids = [sid for sid, _ in sequences]
    seqs = [seq for _, seq in sequences]
    
    num_sequences = len(sequences)
    alignment_length = len(seqs[0])
    
    # Calculate statistics
    conserved_count = 0
    variable_count = 0
    gap_only_count = 0
    
    conservation_data = []
    
    for pos in range(alignment_length):
        column = ''.join(seq[pos] for seq in seqs)
        
        gap_fraction = column.count('-') / num_sequences
        entropy = calculate_entropy(column)
        category = get_conservation_category(column)
        most_common = get_most_common_residue(column)
        
        conservation_data.append({
            'position': pos + 1,
            'gap_fraction': gap_fraction,
            'entropy': entropy,
            'category': category,
            'most_common': most_common
        })
        
        if category == 'CONSERVED':
            conserved_count += 1
        elif category == 'VARIABLE':
            variable_count += 1
        elif category == 'GAP_ONLY':
            gap_only_count += 1
    
    # Calculate pairwise identities
    identity_matrix = {}
    total_identity = 0.0
    count = 0
    
    for i, id1 in enumerate(seq_ids):
        identity_matrix[id1] = {}
        for j, id2 in enumerate(seq_ids):
            identity = calculate_pairwise_identity(seqs[i], seqs[j])
            identity_matrix[id1][id2] = identity
            
            if i < j:
                total_identity += identity
                count += 1
    
    avg_identity = total_identity / count if count > 0 else 0.0
    
    # Write statistics file
    stats_file = os.path.join(OUTPUT_DIR, f"alignment_{aln_num}_stats.txt")
    with open(stats_file, 'w') as f:
        f.write("Alignment Statistics\n")
        f.write("====================\n")
        f.write(f"Number of sequences: {num_sequences}\n")
        f.write(f"Alignment length: {alignment_length}\n")
        f.write(f"Conserved positions: {conserved_count}\n")
        f.write(f"Variable positions: {variable_count}\n")
        f.write(f"Gap-only columns: {gap_only_count}\n")
        f.write(f"Average pairwise identity: {avg_identity:.2f}%\n")
    
    # Write conservation file
    conservation_file = os.path.join(OUTPUT_DIR, f"alignment_{aln_num}_conservation.tsv")
    with open(conservation_file, 'w') as f:
        f.write("Position\tGap_Fraction\tEntropy\tCategory\tMost_Common\n")
        for data in conservation_data:
            f.write(f"{data['position']}\t{data['gap_fraction']:.6f}\t{data['entropy']:.6f}\t{data['category']}\t{data['most_common']}\n")
    
    # Write identity matrix file
    identity_file = os.path.join(OUTPUT_DIR, f"alignment_{aln_num}_identity.tsv")
    with open(identity_file, 'w') as f:
        # Header
        f.write("\t" + "\t".join(seq_ids) + "\n")
        
        # Rows
        for id1 in seq_ids:
            row = [id1] + [f"{identity_matrix[id1][id2]:.6f}" for id2 in seq_ids]
            f.write("\t".join(row) + "\n")

def main():
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Parse alignments
    alignments = parse_alignments(INPUT_FILE)
    print(f"Found {len(alignments)} alignments")
    
    # Process each alignment
    for i, sequences in enumerate(alignments, 1):
        print(f"Processing alignment {i} with {len(sequences)} sequences...")
        process_alignment(sequences, i)
    
    print(f"Complete! Generated {len(alignments) * 3} files in {OUTPUT_DIR}/")

if __name__ == "__main__":
    main()
