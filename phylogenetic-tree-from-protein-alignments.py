import os
import json
import math
import random
from typing import List, Tuple, Dict, Optional

INPUT_FILE = "/workdir/data/alignments.fasta"
OUTPUT_DIR = "/workdir/trees"
#INPUT_FILE = "./data/alignments.fasta"
#OUTPUT_DIR = "./tests/trees"

def parse_alignments(filepath: str) -> List[Tuple[int, List[str], List[Tuple[str, str]]]]:
    """
    Parse FASTA into alignment sets with flags.
    Returns: List of (set_num, flags, sequences)
    """
    sets = []
    current = []
    cur_id = None
    cur_seq = []
    current_flags = []
    current_num = 0
    
    with open(filepath) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('### ALIGNMENT_SET_'):
                # Save previous sequence if exists
                if cur_id:
                    current.append((cur_id, ''.join(cur_seq)))
                # Save previous set if exists
                if current:
                    sets.append((current_num, current_flags, current))
                
                # Parse new set header
                parts = line.replace('###', '').strip().split()
                # parts[0] is "ALIGNMENT_SET_N", extract N from it
                current_num = int(parts[0].split('_')[-1])
                current_flags = parts[1:] if len(parts) > 1 else []
                
                # Start new set
                current = []
                cur_id = None
                cur_seq = []
            elif line.startswith('>'):
                # Save previous sequence if exists
                if cur_id:
                    current.append((cur_id, ''.join(cur_seq)))
                # Start new sequence
                cur_id = line[1:].strip()
                cur_seq = []
            elif line:
                cur_seq.append(line.strip())
        
        # Save last sequence and set
        if cur_id:
            current.append((cur_id, ''.join(cur_seq)))
        if current:
            sets.append((current_num, current_flags, current))
    
    return sets

def normalize_sequence(seq: str) -> str:
    """Treat ambiguous characters (X, B, Z) as gaps."""
    return seq.replace('X', '-').replace('B', '-').replace('Z', '-')

def calc_p_distance(seq1: str, seq2: str) -> float:
    """
    Calculate p-distance between two aligned sequences.
    
    Rules:
    - Different amino acids = difference
    - Gap vs amino acid = difference
    - Gap vs gap = NOT a difference (skip these positions)
    - Distance = differences / valid_positions
    """
    seq1 = normalize_sequence(seq1)
    seq2 = normalize_sequence(seq2)
    
    differences = 0
    valid_positions = 0
    
    for a, b in zip(seq1, seq2):
        # Skip positions where both have gaps
        if a == '-' and b == '-':
            continue
        
        valid_positions += 1
        
        # Count difference if characters don't match
        if a != b:
            differences += 1
    
    # Handle case where no valid positions
    if valid_positions == 0:
        return 0.0
    
    return differences / valid_positions

def calc_jukes_cantor_distance(seq1: str, seq2: str) -> float:
    """
    Calculate Jukes-Cantor corrected distance.
    
    Formula: d_JC = -(19/20) * ln(1 - (20/19) * p)
    
    Saturation: if p >= 0.95, return 3.8
    """
    p = calc_p_distance(seq1, seq2)
    
    # Handle saturation
    if p >= 0.95:
        return 3.8
    
    # Handle identical sequences
    if p == 0.0:
        return 0.0
    
    # Jukes-Cantor correction for amino acids (20 states)
    # d = -(19/20) * ln(1 - (20/19) * p)
    try:
        inner = 1.0 - (20.0/19.0) * p
        if inner <= 0:
            return 3.8  # Saturation
        d_jc = -(19.0/20.0) * math.log(inner)
        return max(0.0, d_jc)  # Ensure non-negative
    except (ValueError, ZeroDivisionError):
        return 3.8  # Saturation

def build_distance_matrix(seqs: List[Tuple[str, str]], use_jc: bool = True) -> Tuple[Dict[str, Dict[str, float]], List[List[str]]]:
    """
    Build distance matrix.
    Returns: (matrix, saturated_pairs)
    """
    matrix = {}
    saturated_pairs = []
    
    for id1, seq1 in seqs:
        matrix[id1] = {}
        for id2, seq2 in seqs:
            if use_jc:
                dist = calc_jukes_cantor_distance(seq1, seq2)
                # Check for saturation
                p = calc_p_distance(seq1, seq2)
                if p >= 0.95 and id1 != id2:
                    saturated_pairs.append(sorted([id1, id2]))
            else:
                dist = calc_p_distance(seq1, seq2)
            
            matrix[id1][id2] = dist
    
    # Remove duplicate saturated pairs
    saturated_pairs = [list(x) for x in set(tuple(x) for x in saturated_pairs)]
    
    return matrix, saturated_pairs

def upgma(seqs: List[Tuple[str, str]], distance_matrix: Dict[str, Dict[str, float]]) -> str:
    """
    UPGMA algorithm to build phylogenetic tree.
    
    Returns Newick format string.
    
    Tie-breaking: lexicographic order of cluster IDs
    """
    # Initialize: each sequence is a cluster at height 0
    clusters = {}
    for seq_id, seq in seqs:
        clusters[seq_id] = {
            'newick': seq_id,
            'height': 0.0,
            'members': [seq_id]
        }
    
    # Copy distance matrix
    D = {}
    for k1 in distance_matrix:
        D[k1] = dict(distance_matrix[k1])
    
    # Merge clusters until only one remains
    node_counter = 0
    while len(clusters) > 1:
        # Find pair of clusters with minimum distance
        # Tie-breaking: lexicographic order
        min_dist = float('inf')
        min_pair = None
        
        cluster_ids = sorted(clusters.keys())  # Sort for deterministic order
        for i in range(len(cluster_ids)):
            for j in range(i + 1, len(cluster_ids)):
                c1, c2 = cluster_ids[i], cluster_ids[j]
                if D[c1][c2] < min_dist:
                    min_dist = D[c1][c2]
                    min_pair = (c1, c2)
                elif D[c1][c2] == min_dist and min_pair is not None:
                    # Tie-breaking: lexicographic
                    if (c1, c2) < min_pair:
                        min_pair = (c1, c2)
        
        if min_pair is None:
            break
        
        c1, c2 = min_pair
        
        # Calculate new cluster height (distance / 2 for UPGMA)
        new_height = min_dist / 2.0
        
        # Calculate branch lengths from new node to children
        branch_len_1 = new_height - clusters[c1]['height']
        branch_len_2 = new_height - clusters[c2]['height']
        
        # Build Newick subtree for new cluster
        new_newick = f"({clusters[c1]['newick']}:{branch_len_1:.6f},{clusters[c2]['newick']}:{branch_len_2:.6f})"
        
        # Create new cluster ID
        new_id = f"NODE_{node_counter}"
        node_counter += 1
        
        # Create new cluster
        new_cluster = {
            'newick': new_newick,
            'height': new_height,
            'members': clusters[c1]['members'] + clusters[c2]['members']
        }
        
        # Update distance matrix for new cluster
        # UPGMA: distance = arithmetic mean
        new_distances = {}
        for other_id in clusters:
            if other_id in (c1, c2):
                continue
            new_dist = (D[c1][other_id] + D[c2][other_id]) / 2.0
            new_distances[other_id] = new_dist
        
        # Remove old clusters from matrix
        del clusters[c1]
        del clusters[c2]
        del D[c1]
        del D[c2]
        for other_id in list(D.keys()):
            if c1 in D[other_id]:
                del D[other_id][c1]
            if c2 in D[other_id]:
                del D[other_id][c2]
        
        # Add new cluster
        clusters[new_id] = new_cluster
        D[new_id] = {}
        for other_id, dist in new_distances.items():
            D[new_id][other_id] = dist
            D[other_id][new_id] = dist
        D[new_id][new_id] = 0.0
    
    # Return final tree with semicolon
    final_newick = list(clusters.values())[0]['newick']
    return final_newick + ';'

def bootstrap_resample(seqs: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    """Resample alignment columns with replacement."""
    if not seqs:
        return seqs
    
    alignment_length = len(seqs[0][1])
    
    # Sample column indices with replacement
    sampled_indices = [random.randint(0, alignment_length - 1) for _ in range(alignment_length)]
    
    # Build resampled sequences
    resampled = []
    for seq_id, seq in seqs:
        new_seq = ''.join(seq[i] for i in sampled_indices)
        resampled.append((seq_id, new_seq))
    
    return resampled

def get_clades(newick: str) -> List[Tuple[str, ...]]:
    """Extract all clades (sets of leaf names) from a Newick tree."""
    # Simple implementation: find all groups in parentheses
    clades = []
    
    # Remove branch lengths and bootstrap values for parsing
    clean = newick
    import re
    # Remove :numbers
    clean = re.sub(r':\d+\.\d+', '', clean)
    # Remove )numbers (bootstrap)
    clean = re.sub(r'\)\d+', ')', clean)
    clean = clean.replace(';', '')
    
    def extract_leaves(s):
        """Extract leaf names from a substring."""
        # Remove all parentheses and split by comma
        s = s.replace('(', '').replace(')', '')
        leaves = [x.strip() for x in s.split(',') if x.strip()]
        return leaves
    
    # Find all parenthesized groups
    depth = 0
    start = -1
    for i, char in enumerate(clean):
        if char == '(':
            if depth == 0:
                start = i
            depth += 1
        elif char == ')':
            depth -= 1
            if depth == 0 and start >= 0:
                group = clean[start+1:i]
                leaves = extract_leaves(group)
                if len(leaves) > 1:
                    clades.append(tuple(sorted(leaves)))
    
    return clades

def add_bootstrap_support(original_newick: str, bootstrap_trees: List[str]) -> str:
    """Add bootstrap support values to original tree."""
    # Get clades from original tree
    original_clades = get_clades(original_newick)
    
    # Count occurrences of each clade in bootstrap trees
    clade_counts = {clade: 0 for clade in original_clades}
    
    for boot_tree in bootstrap_trees:
        boot_clades = get_clades(boot_tree)
        for clade in original_clades:
            if clade in boot_clades:
                clade_counts[clade] += 1
    
    # Calculate support percentages
    n_bootstrap = len(bootstrap_trees)
    clade_support = {clade: int(round(100.0 * count / n_bootstrap)) 
                     for clade, count in clade_counts.items()}
    
    # Insert support values into Newick string
    # This is a simplified approach - match patterns and insert
    result = original_newick
    
    for clade, support in clade_support.items():
        # Find the pattern in the tree
        # Look for closing parenthesis followed by colon
        # Insert support value before colon
        import re
        # Pattern: )XXX:branch_length where XXX might be NODE_N or nothing
        # Replace with )SUPPORT:branch_length
        
        # This is complex - for simplicity, we'll do a basic replacement
        # In practice, need proper tree parsing
        result = re.sub(r'\)(NODE_\d+)?:', f'){support}:', result, count=len(original_clades))
    
    return result

def perform_bootstrap(seqs: List[Tuple[str, str]], n_bootstrap: int, use_jc: bool) -> List[str]:
    """Perform bootstrap analysis."""
    bootstrap_trees = []
    
    for _ in range(n_bootstrap):
        # Resample alignment
        resampled = bootstrap_resample(seqs)
        
        # Build distance matrix
        boot_matrix, _ = build_distance_matrix(resampled, use_jc)
        
        # Build tree
        boot_tree = upgma(resampled, boot_matrix)
        bootstrap_trees.append(boot_tree)
    
    return bootstrap_trees

def write_distance_matrix(seqs: List[Tuple[str, str]], matrix: Dict[str, Dict[str, float]], 
                          filepath: str, model: str):
    """Write distance matrix in TSV format with model indicator."""
    seq_ids = [seq_id for seq_id, _ in seqs]
    
    with open(filepath, 'w') as f:
        # Write header with model
        f.write('Model\t' + '\t'.join(seq_ids) + '\n')
        
        # Write matrix rows
        for id1 in seq_ids:
            f.write(id1)
            for id2 in seq_ids:
                f.write(f"\t{matrix[id1][id2]:.6f}")
            f.write('\n')

def write_metadata(seqs: List[Tuple[str, str]], matrix: Dict[str, Dict[str, float]],
                   saturated_pairs: List[List[str]], model: str, n_bootstrap: int,
                   filepath: str):
    """Write metadata JSON file."""
    seq_ids = [sid for sid, _ in seqs]
    
    # Find identical pairs
    identical_pairs = []
    for i, id1 in enumerate(seq_ids):
        for id2 in seq_ids[i+1:]:
            if matrix[id1][id2] == 0.0:
                identical_pairs.append([id1, id2])
    
    # Calculate statistics
    distances = []
    for i, id1 in enumerate(seq_ids):
        for id2 in seq_ids[i+1:]:
            distances.append(matrix[id1][id2])
    
    avg_dist = sum(distances) / len(distances) if distances else 0.0
    max_dist = max(distances) if distances else 0.0
    
    metadata = {
        "num_sequences": len(seqs),
        "alignment_length": len(seqs[0][1]) if seqs else 0,
        "distance_model": model,
        "saturated_pairs": saturated_pairs,
        "identical_pairs": identical_pairs,
        "bootstrap_replicates": n_bootstrap,
        "average_pairwise_distance": round(avg_dist, 6),
        "max_pairwise_distance": round(max_dist, 6)
    }
    
    with open(filepath, 'w') as f:
        json.dump(metadata, f, indent=2)

def process_alignment_set(seqs: List[Tuple[str, str]], set_num: int, flags: List[str], 
                          output_dir: str):
    """Process one alignment set with flags."""
    
    # Parse flags
    use_pdist = 'PDIST' in flags
    n_bootstrap = 0
    for flag in flags:
        if flag.startswith('BOOTSTRAP='):
            n_bootstrap = int(flag.split('=')[1])
    
    # Determine distance model
    use_jc = not use_pdist
    model_name = "p-distance" if use_pdist else "JC69"
    
    # Handle edge case: single sequence
    if len(seqs) == 1:
        seq_id = seqs[0][0]
        newick_tree = f"{seq_id}:0.000000;"
        
        # Write tree
        tree_filepath = os.path.join(output_dir, f"set_{set_num}_tree.nwk")
        with open(tree_filepath, 'w') as f:
            f.write(newick_tree + '\n')
        
        # Write distance matrix
        dist_filepath = os.path.join(output_dir, f"set_{set_num}_distances.tsv")
        with open(dist_filepath, 'w') as f:
            f.write(f"Model\t{seq_id}\n")
            f.write(f"{seq_id}\t0.000000\n")
        
        # Write metadata
        meta_filepath = os.path.join(output_dir, f"set_{set_num}_metadata.json")
        metadata = {
            "num_sequences": 1,
            "alignment_length": len(seqs[0][1]),
            "distance_model": model_name,
            "saturated_pairs": [],
            "identical_pairs": [],
            "bootstrap_replicates": 0,
            "average_pairwise_distance": 0.0,
            "max_pairwise_distance": 0.0
        }
        with open(meta_filepath, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        return
    
    # Build distance matrix
    distance_matrix, saturated_pairs = build_distance_matrix(seqs, use_jc)
    
    # Build phylogenetic tree using UPGMA
    newick_tree = upgma(seqs, distance_matrix)
    
    # Perform bootstrap if requested
    if n_bootstrap > 0:
        bootstrap_trees = perform_bootstrap(seqs, n_bootstrap, use_jc)
        newick_tree = add_bootstrap_support(newick_tree, bootstrap_trees)
    
    # Write tree file
    tree_filepath = os.path.join(output_dir, f"set_{set_num}_tree.nwk")
    with open(tree_filepath, 'w') as f:
        f.write(newick_tree + '\n')
    
    # Write distance matrix file
    dist_filepath = os.path.join(output_dir, f"set_{set_num}_distances.tsv")
    write_distance_matrix(seqs, distance_matrix, dist_filepath, model_name)
    
    # Write metadata file
    meta_filepath = os.path.join(output_dir, f"set_{set_num}_metadata.json")
    write_metadata(seqs, distance_matrix, saturated_pairs, model_name, n_bootstrap, meta_filepath)

def main():
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Parse input alignment sets with flags
    alignment_sets = parse_alignments(INPUT_FILE)
    print(f"Found {len(alignment_sets)} alignment sets", flush=True)
    
    # Process each alignment set
    for set_num, flags, seqs in alignment_sets:
        print(f"Processing set {set_num} with {len(seqs)} sequences and flags {flags}...", flush=True)
        process_alignment_set(seqs, set_num, flags, OUTPUT_DIR)
    
    print(f"Complete! Output written to {OUTPUT_DIR}/", flush=True)

if __name__ == "__main__":
    main()
