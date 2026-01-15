import os
import json
import math
import re
from collections import defaultdict
from typing import Dict, List, Tuple, Set
from itertools import combinations

# Paths
#DATA_DIR = "./data"
DATA_DIR = "/workdir/data"
ALIGNMENTS_DIR = os.path.join(DATA_DIR, "alignments")
REFERENCE = os.path.join(DATA_DIR, "reference.fasta")
POPULATIONS_FILE = os.path.join(DATA_DIR, "populations.txt")
#RESULTS_DIR = "./tests/expected_results"
RESULTS_DIR = "/workdir/results"

# Output files
VCF_RAW = os.path.join(RESULTS_DIR, "variants.vcf")
VCF_FILTERED = os.path.join(RESULTS_DIR, "variants_filtered.vcf")
STATS_JSON = os.path.join(RESULTS_DIR, "population_stats.json")
DISTANCE_TSV = os.path.join(RESULTS_DIR, "distance_matrix.tsv")
PHYLOGENY_NWK = os.path.join(RESULTS_DIR, "phylogeny.nwk")
SUMMARY_TXT = os.path.join(RESULTS_DIR, "summary.txt")

# Create results directory
os.makedirs(RESULTS_DIR, exist_ok=True)

# ============================================================================
# PART 1: PARSE INPUT DATA
# ============================================================================

def parse_fasta(filepath):
    """Parse FASTA file and return dict of sequences"""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def parse_populations(filepath):
    """Parse population mapping file"""
    pop_map = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 2:
                    sample, population = parts[0], parts[1]
                    pop_map[sample] = population
    return pop_map

def parse_simple_bam(bam_file, min_baseq=20):
    """
    Parse simplified BAM format
    Format: chrom pos ref_base read_bases qualities
    """
    pileup = defaultdict(lambda: defaultdict(list))
    
    with open(bam_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 5:
                continue
            
            chrom = parts[0]
            pos = int(parts[1])
            ref_base = parts[2].upper()
            read_bases = parts[3].upper()
            qualities = parts[4]
            
            # Process read bases
            for read_base, qual_char in zip(read_bases, qualities):
                qual = ord(qual_char) - 33  # Phred+33
                
                if qual >= min_baseq and read_base in 'ACGT':
                    pileup[chrom][pos].append(read_base)
    
    return pileup

# ============================================================================
# PART 2: VARIANT CALLING
# ============================================================================

def call_variants(bam_files, reference_seqs, sample_names, min_depth=5):
    """Call variants from multiple BAM files"""
    
    # Parse all BAMs
    all_pileups = {}
    for bam_file, sample in zip(bam_files, sample_names):
        pileup = parse_simple_bam(bam_file)
        all_pileups[sample] = pileup
    
    # Find all positions with coverage across all samples
    all_positions = set()
    for sample_pileup in all_pileups.values():
        for chrom in sample_pileup:
            for pos in sample_pileup[chrom]:
                all_positions.add((chrom, pos))
    
    variants = []
    
    # Call variants at each position
    for chrom, pos in sorted(all_positions):
        if chrom not in reference_seqs:
            continue
            
        ref_seq = reference_seqs[chrom]
        if pos < 1 or pos > len(ref_seq):
            continue
            
        ref_base = ref_seq[pos-1]
        
        if ref_base not in 'ACGT':
            continue
        
        # Count alleles across all samples
        allele_counts = defaultdict(int)
        sample_genotypes = {}
        
        for sample in sample_names:
            if chrom in all_pileups[sample] and pos in all_pileups[sample][chrom]:
                bases = all_pileups[sample][chrom][pos]
                
                if len(bases) < min_depth:
                    sample_genotypes[sample] = ('./.',  0, 0)
                    continue
                
                # Count bases
                base_counts = defaultdict(int)
                for base in bases:
                    base_counts[base] += 1
                    allele_counts[base] += 1
                
                # Call genotype
                total = len(bases)
                ref_count = base_counts.get(ref_base, 0)
                ref_freq = ref_count / total
                
                # Find most common non-ref allele
                non_ref_alleles = [(b, c) for b, c in base_counts.items() if b != ref_base]
                
                # Simple genotype calling
                if ref_freq >= 0.9:
                    gt = '0/0'
                    gq = min(99, int(-10 * math.log10(max(0.001, 1 - ref_freq))))
                elif ref_freq <= 0.1 and non_ref_alleles:
                    gt = '1/1'
                    gq = min(99, int(-10 * math.log10(max(0.001, ref_freq))))
                else:
                    gt = '0/1'
                    deviation = abs(ref_freq - 0.5)
                    gq = min(99, int(-10 * math.log10(max(0.001, deviation * 2))))
                
                sample_genotypes[sample] = (gt, total, gq)
            else:
                sample_genotypes[sample] = ('./.',  0, 0)
        
        # Find alternate allele (most common non-ref)
        non_ref_alleles = [(allele, count) for allele, count in allele_counts.items() 
                          if allele != ref_base]
        
        if non_ref_alleles:
            alt_base = max(non_ref_alleles, key=lambda x: x[1])[0]
            
            # Only keep if at least one sample has non-reference allele
            has_alt = any(gt != '0/0' and gt != './.' 
                         for gt, _, _ in sample_genotypes.values())
            
            if has_alt:
                variants.append((chrom, pos, ref_base, alt_base, sample_genotypes))
    
    return variants

def write_vcf(variants, sample_names, output_file, reference_seqs):
    """Write variants to VCF format"""
    with open(output_file, 'w') as f:
        # Header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=CustomVariantCaller\n")
        
        # Contigs
        for chrom, seq in reference_seqs.items():
            f.write(f"##contig=<ID={chrom},length={len(seq)}>\n")
        
        # INFO fields
        f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        f.write("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count\">\n")
        f.write("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele Number\">\n")
        f.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
        
        # FORMAT fields
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
        f.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
        
        # Column header
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n")
        
        # Variants
        for chrom, pos, ref, alt, sample_gts in variants:
            # Calculate INFO fields
            total_depth = sum(dp for _, dp, _ in sample_gts.values())
            
            # Count alleles
            ac = 0  # alt allele count
            an = 0  # total allele number
            
            for sample in sample_names:
                gt, dp, gq = sample_gts[sample]
                if gt != './.':
                    alleles = gt.split('/')
                    an += 2
                    ac += sum(1 for a in alleles if a == '1')
            
            af = ac / an if an > 0 else 0.0
            
            info = f"DP={total_depth};AC={ac};AN={an};AF={af:.4f}"
            
            # Format genotypes
            format_str = "GT:DP:GQ"
            gt_strings = []
            for sample in sample_names:
                gt, dp, gq = sample_gts[sample]
                gt_strings.append(f"{gt}:{dp}:{gq}")
            
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\t{format_str}\t" + "\t".join(gt_strings) + "\n")

# ============================================================================
# PART 3: VARIANT FILTERING
# ============================================================================

def parse_vcf(vcf_file):
    """Parse VCF file"""
    variants = []
    samples = []
    
    with open(vcf_file) as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                parts = line.split('\t')
                samples = parts[9:]
            elif line:
                parts = line.split('\t')
                chrom = parts[0]
                pos = int(parts[1])
                ref = parts[3]
                alt = parts[4]
                info = parts[7]
                format_fields = parts[8].split(':')
                
                # Parse INFO
                info_dict = {}
                for item in info.split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info_dict[key] = value
                
                # Parse genotypes
                genotypes = {}
                for i, sample in enumerate(samples):
                    gt_data = parts[9 + i].split(':')
                    gt_dict = dict(zip(format_fields, gt_data))
                    genotypes[sample] = gt_dict
                
                variants.append({
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'info': info_dict,
                    'genotypes': genotypes
                })
    
    return variants, samples

def filter_variants(variants, samples):
    """Filter variants based on quality criteria"""
    filtered = []
    
    for var in variants:
        # Check DP
        dp = int(var['info'].get('DP', 0))
        if dp < 20 or dp > 500:
            continue
        
        # Check AF
        af = float(var['info'].get('AF', 0))
        if af < 0.05 or af > 0.95:
            continue
        
        # Check GQ
        high_gq_count = 0
        missing_count = 0
        
        for sample in samples:
            gt_data = var['genotypes'][sample]
            gt = gt_data.get('GT', './.')
            gq = int(gt_data.get('GQ', 0))
            
            if gt == './.':
                missing_count += 1
            elif gq >= 20:
                high_gq_count += 1
        
        # At least 80% with GQ >= 20
        if high_gq_count < 0.8 * len(samples):
            continue
        
        # At most 20% missing
        if missing_count > 0.2 * len(samples):
            continue
        
        # Check biallelic
        if ',' in var['alt']:
            continue
        
        filtered.append(var)
    
    return filtered

def write_filtered_vcf(variants, samples, output_file, header_file):
    """Write filtered variants to VCF"""
    # Copy header from original VCF
    header_lines = []
    with open(header_file) as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                break
    
    with open(output_file, 'w') as f:
        # Write header
        for line in header_lines:
            f.write(line)
        
        # Write variants
        for var in variants:
            info_str = ';'.join(f"{k}={v}" for k, v in var['info'].items())
            
            gt_strings = []
            for sample in samples:
                gt_data = var['genotypes'][sample]
                gt_str = ':'.join(str(gt_data.get(field, '.')) for field in ['GT', 'DP', 'GQ'])
                gt_strings.append(gt_str)
            
            f.write(f"{var['chrom']}\t{var['pos']}\t.\t{var['ref']}\t{var['alt']}\t.\tPASS\t{info_str}\tGT:DP:GQ\t")
            f.write('\t'.join(gt_strings) + '\n')

# ============================================================================
# PART 4: POPULATION GENETICS STATISTICS
# ============================================================================

def calculate_nucleotide_diversity(variants, samples):
    """Calculate pi (nucleotide diversity)"""
    if len(samples) < 2:
        return 0.0
    
    total_pi = 0.0
    n_sites = 0
    
    for var in variants:
        # Calculate allele frequencies
        ref_count = 0
        alt_count = 0
        total_alleles = 0
        
        for sample in samples:
            gt = var['genotypes'][sample].get('GT', './.')
            if gt != './.':
                alleles = gt.split('/')
                for allele in alleles:
                    if allele == '0':
                        ref_count += 1
                    elif allele == '1':
                        alt_count += 1
                    total_alleles += 1
        
        if total_alleles > 0:
            p = ref_count / total_alleles
            q = alt_count / total_alleles
            total_pi += 2 * p * q
            n_sites += 1
    
    return total_pi / n_sites if n_sites > 0 else 0.0

def calculate_tajimas_d(variants, samples):
    """Calculate Tajima's D"""
    n = len(samples)
    
    if n < 4:
        return 0.0
    
    # Count segregating sites
    S = len(variants)
    
    if S == 0:
        return 0.0
    
    # Calculate a_n
    a_n = sum(1.0 / i for i in range(1, n))
    
    # Calculate theta_W
    theta_w = S / a_n
    
    # Calculate pi
    pi = calculate_nucleotide_diversity(variants, samples) * S
    
    # Calculate variance components
    a_n_sq = sum(1.0 / (i * i) for i in range(1, n))
    b_n = (n + 1) / (3 * (n - 1))
    c_n = 2 * (n * a_n - 2 * (n - 1)) / ((n - 1) * (n - 2))
    e_n = c_n / a_n
    
    # Variance
    var = (e_n * S + (a_n_sq / (a_n * a_n + a_n_sq)) * (S * (S - 1)))
    
    if var <= 0:
        return 0.0
    
    # Tajima's D
    d = (pi - theta_w) / math.sqrt(var)
    
    return d

def calculate_heterozygosity(variants, samples):
    """Calculate observed and expected heterozygosity"""
    het_count = 0
    total_count = 0
    exp_het_sum = 0.0
    
    for var in variants:
        # Observed heterozygosity
        for sample in samples:
            gt = var['genotypes'][sample].get('GT', './.')
            if gt != './.':
                total_count += 1
                if '0/1' in gt or '1/0' in gt:
                    het_count += 1
        
        # Expected heterozygosity
        ref_count = 0
        alt_count = 0
        total_alleles = 0
        
        for sample in samples:
            gt = var['genotypes'][sample].get('GT', './.')
            if gt != './.':
                alleles = gt.split('/')
                for allele in alleles:
                    if allele == '0':
                        ref_count += 1
                    elif allele == '1':
                        alt_count += 1
                    total_alleles += 1
        
        if total_alleles > 0:
            p = ref_count / total_alleles
            q = alt_count / total_alleles
            exp_het_sum += 2 * p * q
    
    h_obs = het_count / total_count if total_count > 0 else 0.0
    h_exp = exp_het_sum / len(variants) if len(variants) > 0 else 0.0
    
    return h_obs, h_exp

def calculate_fst(variants, samples1, samples2):
    """Calculate Fst between two populations"""
    if len(samples1) < 1 or len(samples2) < 1:
        return 0.0
    
    fst_sum = 0.0
    n_sites = 0
    
    for var in variants:
        # Calculate allele frequencies in each population
        def get_freq(samples):
            ref_count = 0
            total = 0
            for sample in samples:
                gt = var['genotypes'][sample].get('GT', './.')
                if gt != './.':
                    alleles = gt.split('/')
                    for allele in alleles:
                        if allele == '0':
                            ref_count += 1
                        total += 1
            return ref_count / total if total > 0 else 0.0
        
        p1 = get_freq(samples1)
        p2 = get_freq(samples2)
        
        # Average frequency
        p_avg = (p1 + p2) / 2
        q_avg = 1 - p_avg
        
        # Total heterozygosity
        h_t = 2 * p_avg * q_avg
        
        # Average within-population heterozygosity
        h_s = (2 * p1 * (1 - p1) + 2 * p2 * (1 - p2)) / 2
        
        if h_t > 0:
            fst = (h_t - h_s) / h_t
            fst_sum += fst
            n_sites += 1
    
    return fst_sum / n_sites if n_sites > 0 else 0.0

# ============================================================================
# PART 5: DISTANCE MATRIX AND PHYLOGENY
# ============================================================================

def calculate_pairwise_distance(variants, sample1, sample2):
    """Calculate genetic distance between two samples"""
    differences = 0
    total_sites = 0
    
    for var in variants:
        gt1 = var['genotypes'][sample1].get('GT', './.')
        gt2 = var['genotypes'][sample2].get('GT', './.')
        
        if gt1 != './.' and gt2 != './.':
            alleles1 = [int(a) for a in gt1.split('/') if a != '.']
            alleles2 = [int(a) for a in gt2.split('/') if a != '.']
            
            # Count differences
            for a1 in alleles1:
                for a2 in alleles2:
                    if a1 != a2:
                        differences += 0.5
            
            total_sites += 1
    
    return differences / (total_sites * 2) if total_sites > 0 else 0.0

def build_distance_matrix(variants, samples):
    """Build pairwise distance matrix"""
    n = len(samples)
    matrix = [[0.0] * n for _ in range(n)]
    
    for i in range(n):
        for j in range(i + 1, n):
            dist = calculate_pairwise_distance(variants, samples[i], samples[j])
            if dist == 0.0 and i != j:
                dist = 0.000001
            matrix[i][j] = dist
            matrix[j][i] = dist
    
    return matrix

def write_distance_matrix(matrix, samples, output_file):
    """Write distance matrix to TSV file"""
    with open(output_file, 'w') as f:
        # Header row: empty cell, then sample names
        header_parts = [''] + list(samples)
        f.write('\t'.join(header_parts) + '\n')
        
        # Data rows: sample name, then distances
        for i, sample in enumerate(samples):
            row_parts = [sample] + [f"{matrix[i][j]:.6f}" for j in range(len(samples))]
            f.write('\t'.join(row_parts) + '\n')

def upgma(distance_matrix, samples):
    """UPGMA clustering - SIMPLIFIED WORKING VERSION"""
    n = len(samples)
    
    if n == 0:
        return "();"
    if n == 1:
        return f"{samples[0]}:0.0;"
    if n == 2:
        dist = distance_matrix[0][1]
        return f"({samples[0]}:{dist/2:.6f},{samples[1]}:{dist/2:.6f});"
    
    # Copy distance matrix and track active nodes
    dist = [row[:] for row in distance_matrix]
    active = list(range(n))
    node_names = {i: samples[i] for i in range(n)}
    node_heights = {i: 0.0 for i in range(n)}
    node_sizes = {i: 1 for i in range(n)}
    next_id = n
    
    while len(active) > 1:
        # Find minimum distance
        min_dist = float('inf')
        pair = None
        
        for i in range(len(active)):
            for j in range(i + 1, len(active)):
                if dist[active[i]][active[j]] < min_dist:
                    min_dist = dist[active[i]][active[j]]
                    pair = (active[i], active[j])
        
        if pair is None:
            break
        
        i, j = pair
        
        # New cluster height and branch lengths
        new_height = min_dist / 2.0
        branch_i = max(0.000001, new_height - node_heights[i])
        branch_j = max(0.000001, new_height - node_heights[j])
        
        # Create new node
        node_names[next_id] = f"({node_names[i]}:{branch_i:.6f},{node_names[j]}:{branch_j:.6f})"
        node_heights[next_id] = new_height
        node_sizes[next_id] = node_sizes[i] + node_sizes[j]
        
        # Update distances to new cluster
        for k in active:
            if k != i and k != j:
                new_dist = (dist[i][k] * node_sizes[i] + dist[j][k] * node_sizes[j]) / node_sizes[next_id]
                dist[i][k] = dist[k][i] = new_dist
        
        # Remove j, replace i with new cluster
        active.remove(j)
        active.remove(i)
        active.append(next_id)
        
        # Copy distances for new cluster ID
        if next_id >= len(dist):
            # Extend matrix if needed
            for row in dist:
                row.extend([0.0] * (next_id - len(row) + 1))
            dist.extend([[0.0] * len(dist[0]) for _ in range(next_id - len(dist) + 1)])
        
        for k in range(len(dist)):
            dist[next_id][k] = dist[i][k]
            dist[k][next_id] = dist[i][k]
        
        next_id += 1
    
    # Return the tree
    root = active[0]
    return node_names[root] + ";"

# ============================================================================
# PART 6: SUMMARY STATISTICS
# ============================================================================

def calculate_transition_transversion_ratio(variants):
    """Calculate Ts/Tv ratio"""
    transitions = 0
    transversions = 0
    
    transition_pairs = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
    
    for var in variants:
        ref = var['ref']
        alt = var['alt']
        
        if len(ref) == 1 and len(alt) == 1:
            if (ref, alt) in transition_pairs:
                transitions += 1
            else:
                transversions += 1
    
    return transitions / transversions if transversions > 0 else 0.0

def generate_summary(variants_raw, variants_filtered, samples, pop_map, stats, output_file):
    """Generate human-readable summary report"""
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("VARIANT CALLING AND POPULATION GENETICS ANALYSIS SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        # Variant calling summary
        f.write("VARIANT CALLING SUMMARY\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total variants called: {len(variants_raw)}\n")
        f.write(f"Variants after filtering: {len(variants_filtered)}\n")
        f.write(f"Filtering rate: {(1 - len(variants_filtered)/len(variants_raw))*100:.1f}%\n" if len(variants_raw) > 0 else "Filtering rate: N/A\n")
        f.write(f"Transition/Transversion ratio: {stats['global_stats']['transition_transversion_ratio']:.2f}\n")
        f.write(f"Mean depth: {stats['global_stats']['mean_depth']:.1f}\n")
        f.write("\n")
        
        # Per-sample statistics
        f.write("PER-SAMPLE STATISTICS\n")
        f.write("-" * 40 + "\n")
        f.write(f"{'Sample':<15} {'Population':<12} {'Missing%':<10} {'Het%':<10} {'AvgDP':<10}\n")
        
        for sample in samples:
            pop = pop_map.get(sample, 'Unknown')
            
            # Calculate stats
            missing = 0
            het = 0
            total_dp = 0
            count = 0
            
            for var in variants_filtered:
                gt = var['genotypes'][sample].get('GT', './.')
                dp = int(var['genotypes'][sample].get('DP', 0))
                
                if gt == './.':
                    missing += 1
                elif '0/1' in gt or '1/0' in gt:
                    het += 1
                
                total_dp += dp
                count += 1
            
            missing_pct = (missing / count * 100) if count > 0 else 0
            het_pct = (het / (count - missing) * 100) if (count - missing) > 0 else 0
            avg_dp = total_dp / count if count > 0 else 0
            
            f.write(f"{sample:<15} {pop:<12} {missing_pct:<10.1f} {het_pct:<10.1f} {avg_dp:<10.1f}\n")
        
        f.write("\n")
        
        # Population statistics
        f.write("POPULATION GENETICS STATISTICS\n")
        f.write("-" * 40 + "\n")
        
        for pop, pop_stats in stats['populations'].items():
            f.write(f"\nPopulation: {pop}\n")
            f.write(f"  Sample count: {pop_stats['sample_count']}\n")
            f.write(f"  Nucleotide diversity (π): {pop_stats['nucleotide_diversity']:.6f}\n")
            f.write(f"  Tajima's D: {pop_stats['tajimas_d']:.3f}\n")
            f.write(f"  Observed heterozygosity: {pop_stats['observed_heterozygosity']:.4f}\n")
            f.write(f"  Expected heterozygosity: {pop_stats['expected_heterozygosity']:.4f}\n")
            f.write(f"  Segregating sites: {pop_stats['segregating_sites']}\n")
        
        f.write("\n")
        
        # Fst matrix
        f.write("POPULATION DIFFERENTIATION (Fst)\n")
        f.write("-" * 40 + "\n")
        
        populations = list(stats['populations'].keys())
        for i, pop1 in enumerate(populations):
            for pop2 in populations[i+1:]:
                key = f"{pop1}_vs_{pop2}"
                if key in stats['fst_matrix']:
                    fst = stats['fst_matrix'][key]
                    f.write(f"{pop1} vs {pop2}: {fst:.4f}\n")
        
        f.write("\n")
        
        # Interpretation
        f.write("INTERPRETATION NOTES\n")
        f.write("-" * 40 + "\n")
        f.write("Tajima's D interpretation:\n")
        f.write("  D < -2: Excess of rare alleles (positive selection or population expansion)\n")
        f.write("  D ≈ 0: Neutral evolution\n")
        f.write("  D > +2: Excess of intermediate frequency alleles (balancing selection)\n")
        f.write("\n")
        f.write("Fst interpretation:\n")
        f.write("  0.00-0.05: Little differentiation\n")
        f.write("  0.05-0.15: Moderate differentiation\n")
        f.write("  0.15-0.25: Great differentiation\n")
        f.write("  >0.25: Very great differentiation\n")
        f.write("\n")

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    print("Starting variant calling pipeline...", flush=True)
    
    # Parse input data
    print("Parsing reference genome...", flush=True)
    reference_seqs = parse_fasta(REFERENCE)
    
    print("Parsing population mapping...", flush=True)
    pop_map = parse_populations(POPULATIONS_FILE)
    
    # Find BAM files
    bam_files = sorted([os.path.join(ALIGNMENTS_DIR, f) for f in os.listdir(ALIGNMENTS_DIR) if f.endswith('.bam')])
    sample_names = [os.path.basename(f).replace('.bam', '') for f in bam_files]
    
    print(f"Found {len(bam_files)} BAM files", flush=True)
    
    # PART 1: Call variants
    print("Calling variants...", flush=True)
    variants_raw = call_variants(bam_files, reference_seqs, sample_names, min_depth=5)
    print(f"Called {len(variants_raw)} variants", flush=True)
    
    # Write raw VCF
    print("Writing raw VCF...", flush=True)
    write_vcf(variants_raw, sample_names, VCF_RAW, reference_seqs)
    
    # PART 2: Filter variants
    print("Filtering variants...", flush=True)
    variants_parsed, samples_parsed = parse_vcf(VCF_RAW)
    variants_filtered = filter_variants(variants_parsed, samples_parsed)
    print(f"Retained {len(variants_filtered)} variants after filtering", flush=True)
    
    write_filtered_vcf(variants_filtered, samples_parsed, VCF_FILTERED, VCF_RAW)
    
    # PART 3: Population genetics statistics
    print("Calculating population statistics...", flush=True)
    
    # Group samples by population
    pop_samples = defaultdict(list)
    for sample in sample_names:
        pop = pop_map.get(sample, 'Unknown')
        pop_samples[pop].append(sample)
    
    pop_stats = {}
    for pop, samples in pop_samples.items():
        if len(samples) < 2:
            pop_stats[pop] = {
                'sample_count': len(samples),
                'nucleotide_diversity': 0.0,
                'tajimas_d': 0.0,
                'observed_heterozygosity': 0.0,
                'expected_heterozygosity': 0.0,
                'segregating_sites': len(variants_filtered)
            }
        else:
            pi = calculate_nucleotide_diversity(variants_filtered, samples)
            tajd = calculate_tajimas_d(variants_filtered, samples)
            h_obs, h_exp = calculate_heterozygosity(variants_filtered, samples)
            
            pop_stats[pop] = {
                'sample_count': len(samples),
                'nucleotide_diversity': pi,
                'tajimas_d': tajd,
                'observed_heterozygosity': h_obs,
                'expected_heterozygosity': h_exp,
                'segregating_sites': len(variants_filtered)
            }
    
    # Calculate Fst
    fst_matrix = {}
    populations = list(pop_samples.keys())
    for i, pop1 in enumerate(populations):
        for pop2 in populations[i+1:]:
            fst = calculate_fst(variants_filtered, pop_samples[pop1], pop_samples[pop2])
            fst_matrix[f"{pop1}_vs_{pop2}"] = fst
    
    # Global statistics
    mean_depth = sum(int(var['info'].get('DP', 0)) for var in variants_filtered) / len(variants_filtered) if variants_filtered else 0
    tstv = calculate_transition_transversion_ratio(variants_filtered)
    
    stats = {
        'populations': pop_stats,
        'fst_matrix': fst_matrix,
        'global_stats': {
            'total_variants': len(variants_raw),
            'total_filtered': len(variants_filtered),
            'mean_depth': mean_depth,
            'transition_transversion_ratio': tstv
        }
    }
    
    # Write stats JSON
    print("Writing population statistics...", flush=True)
    with open(STATS_JSON, 'w') as f:
        json.dump(stats, f, indent=2)
    
    # PART 4: Distance matrix and phylogeny
    print("Calculating genetic distances...", flush=True)
    distance_matrix = build_distance_matrix(variants_filtered, sample_names)
    
    print("Writing distance matrix...", flush=True)
    write_distance_matrix(distance_matrix, sample_names, DISTANCE_TSV)
    
    with open(DISTANCE_TSV) as f:
        for i, line in enumerate(f):
            print(f"Line {i}: {repr(line)}", flush=True)
            if i < 3:  # Just first 3 lines
                parts = line.strip().split('\t')
                print(f"  Parts: {parts}", flush=True)
                print(f"  Num parts: {len(parts)}", flush=True)
    
        print("Building phylogenetic tree...", flush=True)
        newick_tree = upgma(distance_matrix, sample_names)
    
        with open(PHYLOGENY_NWK, 'w') as f:
            f.write(newick_tree + '\n')
    
    # PART 5: Summary report
    print("Generating summary report...", flush=True)
    generate_summary(variants_raw, variants_filtered, sample_names, pop_map, stats, SUMMARY_TXT)
    
    print("Pipeline complete!", flush=True)
    print(f"Results written to {RESULTS_DIR}/", flush=True)

if __name__ == "__main__":
    main()
