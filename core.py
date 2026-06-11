from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction, molecular_weight
import primer3
import pandas as pd
import io

def extinction_coefficient(seq: str) -> int:
    s = seq.upper()
    return (
        s.count('A') * 15400 +
        s.count('C') * 7400 +
        s.count('G') * 11500 +
        s.count('T') * 8700
    )

def calculate_tm_custom(primer_seq: str, Na: float = 50, Mg: float = 0, dntps: float = 0):
    return mt.Tm_NN(primer_seq, Na=Na, Mg=Mg, dNTPs=dntps)

def has_gc_clamp(seq: str, clamp_nt: int = 1) -> bool:
    s = seq.upper()
    return (s[-clamp_nt:].count('G') + s[-clamp_nt:].count('C')) >= clamp_nt

def primer_properties(primer_seq: str):
    gc_content = round(gc_fraction(primer_seq) * 100, 2)
    mw = molecular_weight(primer_seq, seq_type='DNA')
    ext_coeff = extinction_coefficient(primer_seq)
    nmol_per_od = 1 / ext_coeff * 1e6 if ext_coeff != 0 else None
    ug_per_od = (nmol_per_od * mw) / 1e6 if nmol_per_od is not None else None
    return {
        "GC%": gc_content,
        "Molecular Weight (g/mol)": round(mw, 2),
        "Extinction Coefficient (L/mol/cm)": ext_coeff,
        "nmol/OD260": round(nmol_per_od, 2) if nmol_per_od is not None else None,
        "µg/OD260": round(ug_per_od, 2) if ug_per_od is not None else None,
    }

def check_secondary_structure(primer_seq: str):
    hairpin = primer3.bindings.calcHairpin(primer_seq)
    homodimer = primer3.bindings.calcHomodimer(primer_seq)
    return {
        "Hairpin ΔG (kcal/mol)": getattr(hairpin, "dg", None),
        "Hairpin found": getattr(hairpin, "structure_found", None),
        "Homodimer ΔG (kcal/mol)": getattr(homodimer, "dg", None),
        "Homodimer found": getattr(homodimer, "structure_found", None)
    }

def design_inner_primers_match_tm(seq_record, snp_pos, ref, alt, target_tm, tm_tolerance=1.0, na=50, mg=0, dntps=0):
    seq = str(seq_record.seq).upper()
    snp_index = int(snp_pos) - 1

    for length in range(18, 26):
        if_start = max(0, snp_index - length + 1)
        if_primer_seq = seq[if_start:snp_index] + ref
        ir_downstream = seq[snp_index + 1: snp_index + 1 + (length - 1)]
        ir_primer_seq = alt + ir_downstream
        ir_primer = str(Seq(ir_primer_seq).reverse_complement())

        try:
            if_tm = calculate_tm_custom(if_primer_seq, Na=na, Mg=mg, dntps=dntps)
            ir_tm = calculate_tm_custom(ir_primer, Na=na, Mg=mg, dntps=dntps)
        except Exception:
            continue

        if abs(if_tm - target_tm) <= tm_tolerance and abs(ir_tm - target_tm) <= tm_tolerance:
            return if_primer_seq, round(if_tm, 2), ir_primer, round(ir_tm, 2)

    return None, None, None, None

def calculate_amplicon_sizes(seq_record, snp_pos, ref, alt, of_primer, or_primer, if_primer, ir_primer):
    seq = str(seq_record.seq).upper()
    of_pos = seq.find(of_primer) if of_primer else -1
    rc_or = str(Seq(or_primer).reverse_complement()) if or_primer else None
    or_pos = seq.rfind(rc_or) if rc_or else -1
    if_pos = seq.find(if_primer) if if_primer else -1
    rc_ir = str(Seq(ir_primer).reverse_complement()) if ir_primer else None
    ir_pos = seq.rfind(rc_ir) if rc_ir else -1

    def size_between(left_pos, right_pos, right_primer_len):
        if left_pos != -1 and right_pos != -1 and right_pos > left_pos:
            return right_pos - left_pos + right_primer_len
        return None

    outer_size = size_between(of_pos, or_pos, len(or_primer) if or_primer else 0)
    wt_size = size_between(of_pos, ir_pos, len(ir_primer) if ir_primer else 0)
    mut_size = size_between(if_pos, or_pos, len(or_primer) if or_primer else 0)
    internal_size = size_between(if_pos, ir_pos, len(ir_primer) if ir_primer else 0)

    return {
        "outer_product": outer_size,
        "wild_type_product": wt_size,
        "mutant_product": mut_size,
        "internal_control": internal_size
    }

def process_primers(fasta_content: str, csv_content: str, params: dict):
    snp_df = pd.read_csv(io.StringIO(csv_content))
    required_columns = ['ID', 'pos', 'ref', 'alt', 'outer_tm']
    
    if not all(col in snp_df.columns for col in required_columns):
        raise ValueError(f"CSV file must contain these columns: {', '.join(required_columns)}")

    fasta_records = SeqIO.to_dict(SeqIO.parse(io.StringIO(fasta_content), "fasta"))
    if not fasta_records:
        raise ValueError("No valid sequences found in FASTA file.")

    na = params.get('na', 50.0)
    mg = params.get('mg', 0.0)
    dntps = params.get('dntps', 0.0)
    tm_tolerance = params.get('tm_tolerance', 1.0)
    of_primer = params.get('of_primer', "").strip()
    or_primer = params.get('or_primer', "").strip()

    results = []
    failed_snps = []
    gel_data = []

    for idx, row in snp_df.iterrows():
        seq_id = str(row['ID'])
        if seq_id not in fasta_records:
            failed_snps.append(f"Sequence '{seq_id}' not found")
            continue
        
        snp_pos = int(row['pos'])
        ref = str(row['ref']).upper()
        alt = str(row['alt']).upper()
        outer_tm = float(row['outer_tm'])

        if_primer, if_tm, ir_primer, ir_tm = design_inner_primers_match_tm(
            fasta_records[seq_id], snp_pos, ref, alt, outer_tm,
            tm_tolerance, na, mg, dntps
        )

        if if_primer and ir_primer:
            if_props = primer_properties(if_primer)
            ir_props = primer_properties(ir_primer)
            
            amplicon_sizes = {}
            if of_primer and or_primer:
                amplicon_sizes = calculate_amplicon_sizes(
                    fasta_records[seq_id], snp_pos, ref, alt, 
                    of_primer.upper(), or_primer.upper(), if_primer, ir_primer
                )
                
                gel_data.append({
                    'seq_id': seq_id,
                    'snp_pos': snp_pos,
                    'ref_alt': f"{ref}/{alt}",
                    'amplicon_sizes': amplicon_sizes
                })
            
            results.append({
                'Seq ID': seq_id,
                'SNP Position': snp_pos,
                'Ref/Alt': f"{ref}/{alt}",
                'IF Primer': if_primer,
                'IF Tm': if_tm,
                'IF GC%': if_props['GC%'],
                'IF GC Clamp': "Yes" if has_gc_clamp(if_primer) else "No",
                'IR Primer': ir_primer,
                'IR Tm': ir_tm,
                'IR GC%': ir_props['GC%'],
                'IR GC Clamp': "Yes" if has_gc_clamp(ir_primer) else "No",
                'Target Tm': outer_tm,
                **{f'Amplicon {k}': v for k, v in amplicon_sizes.items()}
            })
        else:
            failed_snps.append(f"SNP at position {snp_pos} ({ref}→{alt}) in {seq_id}")

    return {
        "results": results,
        "failed_snps": failed_snps,
        "gel_data": gel_data,
        "success_count": len(results),
        "total_count": len(snp_df)
    }
