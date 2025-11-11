import streamlit as st
from Bio import SeqIO, Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.ligation import nt_extinction_coefficient
import primer3
import pandas as pd

# ----- Utilities -----

@st.cache_data
def calculate_tm_custom(primer_seq, Na=50, Mg=0, dNTPs=0):
    return mt.Tm_NN(primer_seq, Na=Na, Mg=Mg, dNTPs=dNTPs)

def has_gc_clamp(seq, clamp_nt=1):
    return seq[-clamp_nt:].count('G') + seq[-clamp_nt:].count('C') >= clamp_nt

def primer_properties(primer_seq):
    gc_content = round(gc_fraction(primer_seq) * 100, 2)
    mw = molecular_weight(primer_seq, seq_type='DNA')
    ext_coeff = nt_extinction_coefficient(primer_seq)
    nmol_per_od = 1 / ext_coeff * 1e6
    ug_per_od = (nmol_per_od * mw) / 1e6
    return {
        "GC%": gc_content,
        "Molecular Weight (g/mol)": round(mw, 2),
        "Extinction Coefficient (L/mol/cm)": ext_coeff,
        "nmol/OD260": round(nmol_per_od, 2),
        "µg/OD260": round(ug_per_od, 2),
    }

def check_secondary_structure(primer_seq):
    hairpin = primer3.bindings.calcHairpin(primer_seq)
    homodimer = primer3.bindings.calcHomodimer(primer_seq)
    return {
        "Hairpin ΔG (kcal/mol)": hairpin.dg,
        "Hairpin found": hairpin.structure_found,
        "Homodimer ΔG (kcal/mol)": homodimer.dg,
        "Homodimer found": homodimer.structure_found
    }

# ----- Primer Design Logic -----
def design_inner_primers_match_tm(seq_record, snp_pos, ref, alt, target_tm, tm_tolerance=1.0):
    seq = str(seq_record.seq).upper()
    snp_pos -= 1

    for length in range(18, 26):
        if_start = max(0, snp_pos - length)
        if_primer = seq[if_start:snp_pos] + ref
        if_tm = calculate_tm_custom(if_primer)

        ir_end = min(len(seq), snp_pos + 1 + length)
        ir_primer_seq = alt + seq[snp_pos + 1:ir_end]
        ir_primer = str(Seq.Seq(ir_primer_seq).reverse_complement())
        ir_tm = calculate_tm_custom(ir_primer)

        if abs(if_tm - target_tm) <= tm_tolerance and abs(ir_tm - target_tm) <= tm_tolerance:
            return if_primer, round(if_tm, 2), ir_primer, round(ir_tm, 2)

    return None, None, None, None

# ----- Streamlit App -----
st.title("Tetra-Primer ARMS-PCR Designer with Primer Analysis")

fasta_file = st.file_uploader("Upload FASTA file", type=["fasta", "fa"])
snp_file = st.file_uploader("Upload SNP data (CSV: ID, pos, ref, alt, outer_tm)", type="csv")

na = st.number_input("Na+ (mM)", value=50.0)
mg = st.number_input("Mg2+ (mM)", value=0.0)
dntps = st.number_input("dNTPs (mM)", value=0.0)

if fasta_file and snp_file:
    snp_df = pd.read_csv(snp_file)
    fasta_records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    results = []

    for _, row in snp_df.iterrows():
        seq_id = str(row['ID'])
        if seq_id not in fasta_records:
            continue

        snp_pos = int(row['pos'])
        ref = row['ref']
        alt = row['alt']
        outer_tm = float(row['outer_tm'])

        if_primer, if_tm, ir_primer, ir_tm = design_inner_primers_match_tm(
            fasta_records[seq_id], snp_pos, ref, alt, outer_tm
        )

        if if_primer and ir_primer:
            if_props = primer_properties(if_primer)
            ir_props = primer_properties(ir_primer)
            if_struct = check_secondary_structure(if_primer)
            ir_struct = check_secondary_structure(ir_primer)

            results.append({
                'Seq ID': seq_id,
                'IF Primer': if_primer,
                'IF Tm': if_tm,
                'IF GC Clamp': has_gc_clamp(if_primer),
                **{f'IF {k}': v for k, v in if_props.items()},
                **{f'IF {k}': v for k, v in if_struct.items()},
                'IR Primer': ir_primer,
                'IR Tm': ir_tm,
                'IR GC Clamp': has_gc_clamp(ir_primer),
                **{f'IR {k}': v for k, v in ir_props.items()},
                **{f'IR {k}': v for k, v in ir_struct.items()},
                'Target Tm': outer_tm
            })

    result_df = pd.DataFrame(results)
    st.subheader("Inner Primer Results with Thermodynamic Analysis")
    st.dataframe(result_df)

    csv = result_df.to_csv(index=False)
    st.download_button("Download Results CSV", csv, "inner_primers_results.csv", "text/csv")

