import os, json, math
from typing import Dict, List, Tuple

INPUT_FILE = "/workdir/data/sequences.fna"
OUTPUT_DIR = "/workdir/predictions"

# Full codon relative adaptiveness table for E. coli
ECOLI_CODON_TABLE: Dict[str,float] = {
    "CTG":1.0,"CGC":1.0,"GGC":1.0,"GCC":1.0,"ACC":1.0,
    "GAC":1.0,"GAA":1.0,"TTC":1.0,"AAC":1.0,"CAG":1.0,
    "AAA":1.0,"AGC":1.0,"CCA":1.0,"GTC":1.0,"ATC":1.0,
    "ATG":1.0,"TGG":1.0,
    "TTG":0.5,"CGT":0.7,"GGT":0.6,"GCT":0.8,"ACT":0.7,
    "GAT":0.8,"GAG":0.6,"TTT":0.7,"AAT":0.7,"CAA":0.5,
    "AAG":0.8,"AGT":0.5,"CCG":0.7,"GTG":0.8,"ATT":0.8,
    "CTA":0.2,"CGA":0.3,"GGA":0.4,"GCA":0.5,"ACA":0.5,
    "CCC":0.3,"GTT":0.6,"ATA":0.3,
    "TTA":0.05,"AGG":0.05,"AGA":0.05,"TCG":0.1,"GTA":0.2,
    "TCA":0.5,"TCT":0.6,"TCC":0.6,"CCT":0.5,
    "GCG":0.6,"ACG":0.4,"CGG":0.3,"GGG":0.4,
    "CTC":0.8,"TGC":0.8,"TGT":0.6,"CAC":0.7,"CAT":0.6,
    "TAA":1.0,"TAG":0.5,"TGA":0.3
}

# Fill missing codons with default value
bases = ['A','C','G','T']
for a in bases:
    for b in bases:
        for c in bases:
            cod = a+b+c
            if cod not in ECOLI_CODON_TABLE:
                ECOLI_CODON_TABLE[cod] = 0.3

# Genetic code for translation
CODON_TO_AA = {
 "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
 "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
 "TCT":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S",
 "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
 "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
 "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
 "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
 "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
 "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
 "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
 "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
 "CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R",
 "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

def reverse_complement(s: str) -> str:
    comp = {'A':'T','C':'G','G':'C','T':'A'}
    return ''.join(comp.get(ch,'N') for ch in reversed(s))

def translate(dna: str)->str:
    aa=[]
    for i in range(0,len(dna)-2,3):
        aa.append(CODON_TO_AA.get(dna[i:i+3],'X'))
    return ''.join(aa)

def cai_of_dna(dna: str)->float:
    vals=[]
    for i in range(0,len(dna)-2,3):
        c = dna[i:i+3]
        w = ECOLI_CODON_TABLE.get(c,0.0)
        if w>0:
            vals.append(math.log(w))
    if not vals:
        return 0.0
    return math.exp(sum(vals)/len(vals))

def hamming(a:str,b:str)->int:
    return sum(x!=y for x,y in zip(a,b))

def score_rbs(full_seq: str, atg_pos0: int) -> Tuple[int,str,int]:
    """
    Score RBS upstream of ATG.
    Returns (score, rbs_sequence, position_relative_to_atg)
    """
    best = (0,'',0)
    consensus = "AGGAGGU"
    
    for offset in range(5,16):  # 5-15bp upstream
        start = atg_pos0 - offset - len(consensus)
        if start < 0:
            continue
        window = full_seq[start:start+len(consensus)]
        if len(window) < len(consensus):
            continue
        
        dist = hamming(window, consensus)
        if dist == 0:
            score = 7
        elif dist == 1:
            score = 5
        elif dist == 2:
            score = 3
        else:
            score = 0
        
        if score > best[0]:
            best = (score, window, -offset)
    
    return best

def find_orfs(seq:str, strand:str)->List[dict]:
    seq_len = len(seq)
    orfs=[]
    
    for frame in range(3):
        i = frame
        while i < seq_len-2:
            cod = seq[i:i+3]
            if cod == "ATG":
                # Search for stop codon
                j = i+3
                found_stop = False
                while j <= seq_len-3:
                    stop = seq[j:j+3]
                    if stop in ("TAA","TAG","TGA"):
                        found_stop = True
                        end = j+3
                        orf_seq = seq[i:end]
                        aa_len = (end - i)//3
                        
                        # Apply filters
                        if aa_len < 100:
                            break
                        
                        cai = cai_of_dna(orf_seq)
                        if cai < 0.25:
                            break
                        
                        score, rbs_seq, rbs_pos = score_rbs(seq, i)
                        if score < 5:
                            break
                        
                        gc = (orf_seq.count('G') + orf_seq.count('C')) / len(orf_seq)
                        protein = translate(orf_seq)
                        
                        orfs.append({
                            "strand": strand,
                            "frame": frame+1,
                            "start": i+1,
                            "end": end,
                            "length_nt": len(orf_seq),
                            "length_aa": aa_len,
                            "cai": round(cai,3),
                            "rbs_score": score,
                            "rbs_sequence": rbs_seq,
                            "rbs_position": rbs_pos,
                            "gc_content": round(gc,3),
                            "start_codon": "ATG",
                            "stop_codon": stop,
                            "protein": protein
                        })
                        break
                    j += 3
                
                if found_stop:
                    i = end  # Skip past this ORF
                else:
                    i += 3
            else:
                i += 3
    
    return orfs

def resolve_overlaps(orfs:List[dict])->List[dict]:
    if not orfs:
        return []
    
    orfs_sorted = sorted(orfs, key=lambda x: x['start'])
    kept=[]
    
    for o in orfs_sorted:
        should_add = True
        to_remove = []
        
        for i, k in enumerate(kept):
            overlap_start = max(o['start'], k['start'])
            overlap_end = min(o['end'], k['end'])
            overlap_len = max(0, overlap_end - overlap_start)
            
            if overlap_len > 0:
                shorter = min(o['length_nt'], k['length_nt'])
                pct = overlap_len / shorter if shorter > 0 else 0
                
                if pct > 0.5:
                    # Keep longer, or higher CAI if same length
                    if o['length_nt'] > k['length_nt']:
                        to_remove.append(i)
                    elif o['length_nt'] == k['length_nt'] and o['cai'] > k['cai']:
                        to_remove.append(i)
                    else:
                        should_add = False
                        break
        
        # Remove marked entries
        for idx in sorted(to_remove, reverse=True):
            kept.pop(idx)
        
        if should_add:
            kept.append(o)
    
    return sorted(kept, key=lambda x: x['start'])

def format_fasta(s:str, width=60)->str:
    return "\n".join(s[i:i+width] for i in range(0,len(s),width))

# Read FASTA
os.makedirs(OUTPUT_DIR, exist_ok=True)

with open(INPUT_FILE) as fh:
    lines = [l.rstrip('\n') for l in fh]

sequences=[]
cur_id=None
cur_seq=[]

for line in lines:
    if line.startswith(">"):
        if cur_id:
            sequences.append((cur_id,''.join(cur_seq)))
        cur_id = line[1:].split()[0]
        cur_seq=[]
    else:
        cur_seq.append(line.strip())

if cur_id:
    sequences.append((cur_id,''.join(cur_seq)))

all_results=[]
summary=[]
codon_data = {"reference_table": ECOLI_CODON_TABLE, "sequences": {}}

for sid, seq in sequences:
    s = seq.upper()
    
    # Find ORFs on both strands
    forward = find_orfs(s, "+")
    rev_seq = reverse_complement(s)
    reverse = find_orfs(rev_seq, "-")
    
    # Adjust reverse strand coordinates
    adjusted_rev = []
    L = len(s)
    for r in reverse:
        new_end = L - r['start'] + 1
        new_start = L - r['end'] + 1
        r['start'] = new_start
        r['end'] = new_end
        adjusted_rev.append(r)
    
    # Calculate strand scores
    fscore = sum(o['length_nt'] for o in forward) / len(s) if len(s) > 0 else 0
    rscore = sum(o['length_nt'] for o in adjusted_rev) / len(s) if len(s) > 0 else 0
    
    # Select strand
    if fscore >= rscore:
        predicted = "+"
        selected = forward
    else:
        predicted = "-"
        selected = adjusted_rev
    
    # Resolve overlaps
    final = resolve_overlaps(selected)
    
    # Write genes FASTA
    genes_file = os.path.join(OUTPUT_DIR, f"{sid}_genes.faa")
    with open(genes_file,"w") as gf:
        for i, g in enumerate(final,1):
            header = f">{sid}|gene{i}|{g['strand']}|{g['start']}..{g['end']}|{g['length_aa']}aa|CAI={g['cai']:.2f}|RBS={g['rbs_score']}"
            gf.write(header+"\n")
            gf.write(format_fasta(g['protein'])+"\n")
    
    # Write metadata JSON
    metadata = {
        "sequence_id": sid,
        "length": len(s),
        "predicted_strand": predicted,
        "strand_score": {"forward": round(fscore,4), "reverse": round(rscore,4)},
        "genes": [
            {
                "gene_id": f"gene{i+1}",
                "strand": g["strand"],
                "frame": g["frame"],
                "start": g["start"],
                "end": g["end"],
                "length_nt": g["length_nt"],
                "length_aa": g["length_aa"],
                "cai": g["cai"],
                "rbs_sequence": g["rbs_sequence"],
                "rbs_score": g["rbs_score"],
                "rbs_position": g["rbs_position"],
                "gc_content": g["gc_content"],
                "start_codon": g["start_codon"],
                "stop_codon": g["stop_codon"]
            }
            for i,g in enumerate(final)
        ]
    }
    
    with open(os.path.join(OUTPUT_DIR, f"{sid}_metadata.json"), "w") as mf:
        json.dump(metadata, mf, indent=2)
    
    all_results.append(metadata)
    
    # Summary entry
    summary.append({
        "sequence_id": sid,
        "num_genes": len(final),
        "predicted_strand": predicted,
        "avg_cai": round(sum(g['cai'] for g in final)/len(final),2) if final else 0.0,
        "avg_length": round(sum(g['length_aa'] for g in final)/len(final),1) if final else 0.0,
        "total_coding_pct": round(sum(g['length_nt'] for g in final)/len(s),3) if len(s)>0 else 0.0
    })
    
    # Codon usage: count codons in predicted genes
    observed={}
    for g in final:
        sidx = g['start'] - 1
        eidx = g['end']
        region = s[sidx:eidx]
        for i in range(0,len(region)-2,3):
            cod = region[i:i+3]
            if len(cod)==3 and all(ch in "ACGT" for ch in cod) and cod in ECOLI_CODON_TABLE:
                observed[cod] = observed.get(cod,0)+1
    
    # Fallback: if no genes, count codons in full sequence
    if not observed:
        for i in range(0,len(s)-2,3):
            cod = s[i:i+3]
            if len(cod)==3 and all(ch in "ACGT" for ch in cod) and cod in ECOLI_CODON_TABLE:
                observed[cod] = observed.get(cod,0)+1
    
    codon_data['sequences'][sid] = {
        "observed_usage": observed,
        "relative_adaptiveness": ECOLI_CODON_TABLE
    }

# Write summary TSV
with open(os.path.join(OUTPUT_DIR,"summary.tsv"),"w") as sf:
    sf.write("sequence_id\tnum_genes\tpredicted_strand\tavg_cai\tavg_length\ttotal_coding_pct\n")
    for r in summary:
        sf.write(f"{r['sequence_id']}\t{r['num_genes']}\t{r['predicted_strand']}\t{r['avg_cai']}\t{r['avg_length']}\t{r['total_coding_pct']}\n")

# Write codon usage JSON
with open(os.path.join(OUTPUT_DIR,"codon_usage.json"),"w") as cf:
    json.dump(codon_data, cf, indent=2)

print(f"Processed {len(sequences)} sequences")
print(f"Total genes predicted: {sum(x['num_genes'] for x in summary)}")
