import os
import json
from typing import List, Tuple, Dict, Optional

#INPUT_FILE = "./data/sequences.fasta"
#OUTPUT_DIR = "./tests/expected_output"
#BLOSUM62_FILE = "./data/blosum62.txt"
#PAM250_FILE = "./data/pam250.txt"
INPUT_FILE = "/workdir/data/sequences.fasta"
OUTPUT_DIR = "/workdir/alignments"
BLOSUM62_FILE = "/workdir/data/blosum62.txt"
PAM250_FILE = "/workdir/data/pam250.txt"

NEG_INF = float('-inf')

def parse_matrix_file(filepath: str) -> Dict[Tuple[str, str], int]:
    """Parse substitution matrix from file."""
    matrix = {}
    
    with open(filepath) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
    
    # First line is header
    header = lines[0].split()
    
    # Parse data lines
    for line in lines[1:]:
        parts = line.split()
        aa = parts[0]
        scores = [int(x) for x in parts[1:]]
        
        for i, score in enumerate(scores):
            matrix[(aa, header[i])] = score
            matrix[(header[i], aa)] = score  # Symmetric
    
    return matrix

def parse_tasks(filepath: str) -> List[Tuple[int, Dict[str, str], List[Tuple[str, str]]]]:
    """
    Parse alignment tasks from FASTA.
    Returns: List of (task_id, parameters, [(id1, seq1), (id2, seq2)])
    """
    tasks = []
    current_seqs = []
    current_params = {}
    current_id = 0
    cur_seq_id = None
    cur_seq = []
    
    with open(filepath) as f:
        for line in f:
            line = line.rstrip('\n')
            
            if line.startswith('### TASK_'):
                # Save previous task
                if cur_seq_id:
                    current_seqs.append((cur_seq_id, ''.join(cur_seq)))
                if current_seqs:
                    tasks.append((current_id, current_params, current_seqs))
                
                # Parse new task
                parts = line.replace('###', '').strip().split()
                current_id = int(parts[0].split('_')[1])
                
                current_params = {}
                for part in parts[1:]:
                    if '=' in part:
                        key, val = part.split('=')
                        current_params[key] = val
                
                current_seqs = []
                cur_seq_id = None
                cur_seq = []
                
            elif line.startswith('>'):
                if cur_seq_id:
                    current_seqs.append((cur_seq_id, ''.join(cur_seq)))
                cur_seq_id = line[1:].strip()
                cur_seq = []
            elif line:
                cur_seq.append(line.strip())
        
        # Save last task
        if cur_seq_id:
            current_seqs.append((cur_seq_id, ''.join(cur_seq)))
        if current_seqs:
            tasks.append((current_id, current_params, current_seqs))
    
    return tasks

def get_score(matrix: Optional[Dict], a: str, b: str, is_simple: bool) -> int:
    """Get substitution score."""
    if is_simple:
        return 2 if a == b else -1
    else:
        return matrix.get((a, b), matrix.get((b, a), -4))

def needleman_wunsch_affine(seq1: str, seq2: str, matrix: Optional[Dict], 
                             gap_open: int, gap_extend: int, is_simple: bool) -> Tuple[float, str, str]:
    """
    Needleman-Wunsch with affine gap penalties.
    Returns: (score, aligned_seq1, aligned_seq2)
    """
    m, n = len(seq1), len(seq2)
    
    # Three matrices
    M = [[NEG_INF]*(n+1) for _ in range(m+1)]
    Ix = [[NEG_INF]*(n+1) for _ in range(m+1)]
    Iy = [[NEG_INF]*(n+1) for _ in range(m+1)]
    
    # Initialization
    M[0][0] = 0
    
    for i in range(1, m+1):
        Ix[i][0] = gap_open + (i-1) * gap_extend
    
    for j in range(1, n+1):
        Iy[0][j] = gap_open + (j-1) * gap_extend
    
    # Fill matrices
    for i in range(1, m+1):
        for j in range(1, n+1):
            # M matrix
            match_score = get_score(matrix, seq1[i-1], seq2[j-1], is_simple)
            M[i][j] = match_score + max(M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1])
            
            # Ix matrix (gap in seq1)
            Ix[i][j] = max(
                M[i-1][j] + gap_open,
                Ix[i-1][j] + gap_extend
            )
            
            # Iy matrix (gap in seq2)
            Iy[i][j] = max(
                M[i][j-1] + gap_open,
                Iy[i][j-1] + gap_extend
            )
    
    # Optimal score
    optimal_score = max(M[m][n], Ix[m][n], Iy[m][n])
    
    # Traceback
    align1, align2 = [], []
    i, j = m, n
    
    # Determine starting matrix
    if M[m][n] == optimal_score:
        current_matrix = 'M'
    elif Ix[m][n] == optimal_score:
        current_matrix = 'Ix'
    else:
        current_matrix = 'Iy'
    
    while i > 0 or j > 0:
        if current_matrix == 'M':
            if i == 0 or j == 0:
                break
            
            match_score = get_score(matrix, seq1[i-1], seq2[j-1], is_simple)
            
            # Priority: M > Ix > Iy
            if M[i][j] == match_score + M[i-1][j-1]:
                align1.append(seq1[i-1])
                align2.append(seq2[j-1])
                i -= 1
                j -= 1
                current_matrix = 'M'
            elif M[i][j] == match_score + Ix[i-1][j-1]:
                align1.append(seq1[i-1])
                align2.append(seq2[j-1])
                i -= 1
                j -= 1
                current_matrix = 'Ix'
            else:
                align1.append(seq1[i-1])
                align2.append(seq2[j-1])
                i -= 1
                j -= 1
                current_matrix = 'Iy'
        
        elif current_matrix == 'Ix':
            if i == 0:
                break
            
            align1.append(seq1[i-1])
            align2.append('-')
            
            # Determine previous matrix
            if Ix[i][j] == M[i-1][j] + gap_open:
                current_matrix = 'M'
            else:
                current_matrix = 'Ix'
            i -= 1
        
        else:  # Iy
            if j == 0:
                break
            
            align1.append('-')
            align2.append(seq2[j-1])
            
            if Iy[i][j] == M[i][j-1] + gap_open:
                current_matrix = 'M'
            else:
                current_matrix = 'Iy'
            j -= 1
    
    # Handle remaining
    while i > 0:
        align1.append(seq1[i-1])
        align2.append('-')
        i -= 1
    
    while j > 0:
        align1.append('-')
        align2.append(seq2[j-1])
        j -= 1
    
    return optimal_score, ''.join(reversed(align1)), ''.join(reversed(align2))

def smith_waterman_affine(seq1: str, seq2: str, matrix: Optional[Dict],
                          gap_open: int, gap_extend: int, is_simple: bool) -> Tuple[float, str, str, Tuple[int, int], Tuple[int, int]]:
    """
    Smith-Waterman with affine gaps.
    Returns: (score, aligned_seq1, aligned_seq2, start_pos, end_pos)
    """
    m, n = len(seq1), len(seq2)
    
    M = [[0]*(n+1) for _ in range(m+1)]
    Ix = [[0]*(n+1) for _ in range(m+1)]
    Iy = [[0]*(n+1) for _ in range(m+1)]

    max_score = 0
    max_i, max_j = 0, 0
    max_matrix = 'M'
    
    for i in range(1, m+1):
        for j in range(1, n+1):
            match_score = get_score(matrix, seq1[i-1], seq2[j-1], is_simple)
            M[i][j] = max(0, match_score + max(M[i-1][j-1], Ix[i-1][j-1], Iy[i-1][j-1]))
            
            Ix[i][j] = max(0, M[i-1][j] + gap_open, Ix[i-1][j] + gap_extend)
            Iy[i][j] = max(0, M[i][j-1] + gap_open, Iy[i][j-1] + gap_extend)
            
            # Track maximum (deterministic: take first occurrence)
            if M[i][j] > max_score:
                max_score = M[i][j]
                max_i, max_j = i, j
                max_matrix = 'M'
            
            if Ix[i][j] > max_score:
                max_score = Ix[i][j]
                max_i, max_j = i, j
                max_matrix = 'Ix'
            
            if Iy[i][j] > max_score:
                max_score = Iy[i][j]
                max_i, max_j = i, j
                max_matrix = 'Iy'
    
    # Traceback from max
    align1, align2 = [], []
    i, j = max_i, max_j
    current_matrix = max_matrix
    end_i, end_j = i, j
    
    while i > 0 and j > 0:
        if current_matrix == 'M':
            if M[i][j] == 0:
                break
            
            match_score = get_score(matrix, seq1[i-1], seq2[j-1], is_simple)
            
            if M[i][j] == match_score + M[i-1][j-1]:
                current_matrix = 'M'
            elif M[i][j] == match_score + Ix[i-1][j-1]:
                current_matrix = 'Ix'
            else:
                current_matrix = 'Iy'
            
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
        
        elif current_matrix == 'Ix':
            if Ix[i][j] == 0:
                break
            
            align1.append(seq1[i-1])
            align2.append('-')
            
            if Ix[i][j] == M[i-1][j] + gap_open:
                current_matrix = 'M'
            else:
                current_matrix = 'Ix'
            i -= 1
        
        else:  # Iy
            if Iy[i][j] == 0:
                break
            
            align1.append('-')
            align2.append(seq2[j-1])
            
            if Iy[i][j] == M[i][j-1] + gap_open:
                current_matrix = 'M'
            else:
                current_matrix = 'Iy'
            j -= 1
    
    start_i, start_j = i + 1, j + 1
    
    # Handle case where no alignment found (max_score = 0)
    if max_score == 0:
        # No local alignment exists - return empty alignment
        return 0.0, '', '', (0, 0), (0, 0)
    
    return max_score, ''.join(reversed(align1)), ''.join(reversed(align2)), (start_i, start_j), (end_i, end_j)

def semiglobal_alignment(seq1: str, seq2: str, matrix: Optional[Dict],
                         gap_open: int, gap_extend: int, is_simple: bool) -> Tuple[float, str, str]:
    """
    Semi-global alignment: no penalty for end gaps.
    Use standard DP but don't penalize gaps at sequence ends.
    """
    m, n = len(seq1), len(seq2)
    
    # Single DP matrix (simplified approach)
    F = [[0]*(n+1) for _ in range(m+1)]
    
    # No gap penalties for first row/column (free starting gaps)
    for i in range(m+1):
        F[i][0] = 0
    for j in range(n+1):
        F[0][j] = 0
    
    # Fill matrix with standard NW recurrence
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = F[i-1][j-1] + get_score(matrix, seq1[i-1], seq2[j-1], is_simple)
            delete = F[i-1][j] + gap_open  # Gap in seq2
            insert = F[i][j-1] + gap_open   # Gap in seq1
            F[i][j] = max(match, delete, insert)
    
    # Find best score in last row or last column (free ending gaps)
    max_score = F[m][n]
    max_i, max_j = m, n
    
    for i in range(m+1):
        if F[i][n] > max_score:
            max_score = F[i][n]
            max_i, max_j = i, n
    
    for j in range(n+1):
        if F[m][j] > max_score:
            max_score = F[m][j]
            max_i, max_j = m, j
    
    # Traceback from best end position
    align1, align2 = [], []
    i, j = max_i, max_j
    
    # Add free end gaps
    while i < m:
        align1.append(seq1[i])
        align2.append('-')
        i += 1
    
    while j < n:
        align1.append('-')
        align2.append(seq2[j])
        j += 1
    
    # Traceback to origin
    i, j = max_i, max_j
    while i > 0 and j > 0:
        score_diag = F[i-1][j-1] + get_score(matrix, seq1[i-1], seq2[j-1], is_simple)
        score_up = F[i-1][j] + gap_open
        score_left = F[i][j-1] + gap_open
        
        # Priority: diagonal > up > left
        if F[i][j] == score_diag:
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif F[i][j] == score_up:
            align1.append(seq1[i-1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(seq2[j-1])
            j -= 1
    
    # Add free start gaps
    while i > 0:
        align1.append(seq1[i-1])
        align2.append('-')
        i -= 1
    
    while j > 0:
        align1.append('-')
        align2.append(seq2[j-1])
        j -= 1
    
    return max_score, ''.join(reversed(align1)), ''.join(reversed(align2))

def calculate_statistics(seq1_id: str, seq2_id: str, seq1_orig: str, seq2_orig: str,
                         align1: str, align2: str, score: float, task_id: int,
                         algo: str, matrix_name: str, gap_open: int, gap_extend: int,
                         start_pos: Optional[Tuple[int, int]] = None,
                         end_pos: Optional[Tuple[int, int]] = None) -> Dict:
    '''Calculate alignment statistics.'''
    
    num_matches = sum(1 for a, b in zip(align1, align2) if a == b and a != '-')
    num_mismatches = sum(1 for a, b in zip(align1, align2) if a != b and a != '-' and b != '-')
    num_gaps = align1.count('-') + align2.count('-')
    
    # Count gap opens and extends
    gap_opens = 0
    gap_extends = 0
    
    in_gap = False
    for char in align1:
        if char == '-':
            if not in_gap:
                gap_opens += 1
                in_gap = True
            else:
                gap_extends += 1
        else:
            in_gap = False
    
    in_gap = False
    for char in align2:
        if char == '-':
            if not in_gap:
                gap_opens += 1
                in_gap = True
            else:
                gap_extends += 1
        else:
            in_gap = False
    
    alignment_length = len(align1)
    percent_identity = (num_matches / alignment_length * 100) if alignment_length > 0 else 0.0
    
    if start_pos is None:
        start_pos = (1, 1)
    if end_pos is None:
        end_pos = (len(seq1_orig), len(seq2_orig))
    
    return {
        "task_id": task_id,
        "algorithm": algo,
        "matrix": matrix_name,
        "gap_open": gap_open,
        "gap_extend": gap_extend,
        "seq1_id": seq1_id,
        "seq2_id": seq2_id,
        "seq1_length": len(seq1_orig),
        "seq2_length": len(seq2_orig),
        "alignment_length": alignment_length,
        "optimal_score": round(score, 2),
        "num_matches": num_matches,
        "num_mismatches": num_mismatches,
        "num_gaps": num_gaps,
        "num_gap_opens": gap_opens,
        "num_gap_extends": gap_extends,
        "percent_identity": round(percent_identity, 2),
        "alignment_start_pos": list(start_pos),
        "alignment_end_pos": list(end_pos)
    }

def process_task(task_id: int, params: Dict, seqs: List[Tuple[str, str]], 
                 blosum62: Dict, pam250: Dict, output_dir: str):
    """Process one alignment task."""
    
    algo = params['ALGO']
    matrix_name = params['MATRIX']
    gap_open = int(params.get('GAP_OPEN', -5))
    gap_extend = int(params.get('GAP_EXTEND', -1))
    
    seq1_id, seq1 = seqs[0]
    seq2_id, seq2 = seqs[1]
    
    # Select matrix
    if matrix_name == 'SIMPLE':
        matrix = None
        is_simple = True
    elif matrix_name == 'BLOSUM62':
        matrix = blosum62
        is_simple = False
    elif matrix_name == 'PAM250':
        matrix = pam250
        is_simple = False
    else:
        matrix = None
        is_simple = True
    
    # Run alignment
    if algo == 'NW':
        score, align1, align2 = needleman_wunsch_affine(seq1, seq2, matrix, gap_open, gap_extend, is_simple)
        start_pos, end_pos = (1, 1), (len(seq1), len(seq2))
    elif algo == 'SW':
        score, align1, align2, start_pos, end_pos = smith_waterman_affine(seq1, seq2, matrix, gap_open, gap_extend, is_simple)
        # Handle empty SW alignment (no local alignment found)
        if not align1:
            align1 = '-'
            align2 = '-'
    elif algo == 'SEMIGLOBAL':
        score, align1, align2 = semiglobal_alignment(seq1, seq2, matrix, gap_open, gap_extend, is_simple)
        start_pos, end_pos = (1, 1), (len(seq1), len(seq2))
    else:
        raise ValueError(f"Unknown algorithm: {algo}")
    
    # Write outputs
    # Alignment FASTA
    align_file = os.path.join(output_dir, f"alignment_{task_id}.fasta")
    with open(align_file, 'w') as f:
        f.write(f">aligned_{seq1_id}\n{align1}\n")
        f.write(f">aligned_{seq2_id}\n{align2}\n")
    
    # Score
    score_file = os.path.join(output_dir, f"alignment_{task_id}_score.txt")
    with open(score_file, 'w') as f:
        f.write(f"Score: {score:.2f}\n")
    
    # Statistics
    stats = calculate_statistics(seq1_id, seq2_id, seq1, seq2, align1, align2,
                                 score, task_id, algo, matrix_name, gap_open, gap_extend,
                                 start_pos, end_pos)
    stats_file = os.path.join(output_dir, f"alignment_{task_id}_stats.json")
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load matrices
    blosum62 = parse_matrix_file(BLOSUM62_FILE)
    pam250 = parse_matrix_file(PAM250_FILE)
    
    # Parse tasks
    tasks = parse_tasks(INPUT_FILE)
    print(f"Found {len(tasks)} alignment tasks", flush=True)
    
    # Process each task
    for task_id, params, seqs in tasks:
        print(f"Processing task {task_id}: {params['ALGO']} with {params['MATRIX']}...", flush=True)
        process_task(task_id, params, seqs, blosum62, pam250, OUTPUT_DIR)
    
    print(f"Complete! Alignments written to {OUTPUT_DIR}/", flush=True)

if __name__ == "__main__":
    main()
