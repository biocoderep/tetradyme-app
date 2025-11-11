import streamlit as st
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction, molecular_weight
import primer3
import pandas as pd
import io
import plotly.graph_objects as go

# Set page configuration
st.set_page_config(
    page_title="Tetradyme",
    page_icon="",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ----- Utilities -----
def extinction_coefficient(seq: str) -> int:
    s = seq.upper()
    return (
        s.count('A') * 15400 +
        s.count('C') * 7400 +
        s.count('G') * 11500 +
        s.count('T') * 8700
    )

@st.cache_data
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

def create_gel_visualization(band_sizes, lane_labels, title="Expected Gel Electrophoresis Results"):
    fig = go.Figure()
    num_lanes = len(lane_labels)
    
    fig.add_shape(type="rect", x0=-1, y0=0, x1=num_lanes, y1=100,
                  line=dict(color="black", width=2), 
                  fillcolor="rgba(245,245,245,0.95)")
    
    for lane_idx in range(num_lanes):
        fig.add_shape(type="circle", x0=lane_idx-0.2, y0=-5, x1=lane_idx+0.2, y1=-1,
                     line=dict(color="black"), fillcolor="black")

    flat_sizes = [s for lane in band_sizes for s in lane if lane]
    if flat_sizes:
        max_size = max(flat_sizes)
        min_size = min(flat_sizes)
    else:
        max_size, min_size = 100, 50

    range_span = max_size - min_size if max_size != min_size else 1.0

    colors = ['#1e40af', '#6b21a8', '#d97706', '#b91c1c', '#047857']
    
    for lane_idx, (sizes, label) in enumerate(zip(band_sizes, lane_labels)):
        if not sizes:
            fig.add_annotation(x=lane_idx, y=105, text=label, showarrow=False, 
                             font=dict(size=10))
            continue

        fig.add_annotation(x=lane_idx, y=105, text=label, showarrow=False, 
                         font=dict(size=10))
        for i, size in enumerate(sizes):
            position = 100 - ((size - min_size) / range_span) * 90
            color = colors[i % len(colors)]
            
            fig.add_trace(go.Scatter(
                x=[lane_idx-0.3, lane_idx+0.3, lane_idx+0.3, lane_idx-0.3, lane_idx-0.3],
                y=[position, position, position+4, position+4, position],
                fill="toself",
                fillcolor=color,
                line=dict(color="black", width=1),
                hoverinfo="text",
                text=f"Size: {size} bp",
                showlegend=False
            ))
            fig.add_annotation(x=lane_idx, y=position+2, text=str(size), 
                             showarrow=False, font=dict(color="white", size=9, weight="bold"))

    fig.update_layout(
        title=dict(text=title, font=dict(size=18)),
        xaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
        yaxis=dict(title="Fragment Size (bp)", range=[-10, 110], autorange="reversed",
                  gridcolor="lightgray"),
        plot_bgcolor="white",
        height=500,
        showlegend=False
    )
    return fig

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

# ----- Main App -----
col1, col2, col3 = st.columns([1, 2, 1])
with col2:
    st.title("Tetradyme")
    st.subheader("Tetra-Primer ARMS-PCR Designer")
    st.info("Design inner primers for Tetra-Primer ARMS-PCR based on SNP information and target melting temperature.")

# Sidebar
with st.sidebar:
    st.header("Parameters")
    
    st.subheader("Upload Files")
    fasta_file = st.file_uploader("FASTA File", type=["fasta", "fa"])
    snp_file = st.file_uploader("SNP Data (CSV)", type=["csv"])
    
    st.subheader("Reaction Conditions")
    na = st.number_input("Na+ (mM)", value=50.0, min_value=0.0, max_value=1000.0)
    mg = st.number_input("Mg2+ (mM)", value=0.0, min_value=0.0, max_value=100.0)
    dntps = st.number_input("dNTPs (mM)", value=0.0, min_value=0.0, max_value=10.0)
    
    st.subheader("Primer Parameters")
    tm_tolerance = st.number_input("Tm Tolerance (°C)", value=1.0, min_value=0.1, max_value=5.0)
    min_primer_len = st.number_input("Min Primer Length", value=18, min_value=15, max_value=25)
    max_primer_len = st.number_input("Max Primer Length", value=25, min_value=20, max_value=30)
    
    st.subheader("Outer Primers (Optional)")
    of_primer = st.text_input("Outer Forward Primer", value="")
    or_primer = st.text_input("Outer Reverse Primer", value="")

# CSV Format Info
with st.expander("CSV Format Requirements"):
    st.write("Your CSV file must contain the following columns:")
    st.write("- **ID**: Sequence identifier (must match FASTA header)")
    st.write("- **pos**: SNP position (1-based)")
    st.write("- **ref**: Reference allele")
    st.write("- **alt**: Alternative allele")
    st.write("- **outer_tm**: Target melting temperature")
    
    st.code("""ID,pos,ref,alt,outer_tm
sequence_1,10,A,T,60.5
sequence_1,25,G,C,62.0""")

# Main processing
if fasta_file and snp_file:
    try:
        # Read files
        snp_df = pd.read_csv(snp_file)
        required_columns = ['ID', 'pos', 'ref', 'alt', 'outer_tm']
        
        if not all(col in snp_df.columns for col in required_columns):
            st.error(f"CSV file must contain these columns: {', '.join(required_columns)}")
        else:
            fasta_text = io.TextIOWrapper(fasta_file, encoding='utf-8')
            fasta_records = SeqIO.to_dict(SeqIO.parse(fasta_text, "fasta"))
            
            if not fasta_records:
                st.error("No valid sequences found in FASTA file.")
            else:
                # Show sequence info
                st.success(f"Loaded {len(fasta_records)} sequence(s)")
                with st.expander("Sequence Information"):
                    for seq_id, record in fasta_records.items():
                        st.write(f"**{seq_id}**: {len(record.seq)} bp")
                
                # Process button
                if st.button("Design Primers", type="primary"):
                    results = []
                    failed_snps = []
                    gel_data = []
                    
                    progress_bar = st.progress(0)
                    
                    for idx, (_, row) in enumerate(snp_df.iterrows()):
                        progress = (idx + 1) / len(snp_df)
                        progress_bar.progress(progress)
                        
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
                            if_struct = check_secondary_structure(if_primer)
                            ir_struct = check_secondary_structure(ir_primer)
                            
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
                    
                    progress_bar.empty()
                    
                    # Display results
                    if results:
                        st.success(f"Successfully designed primers for {len(results)} SNP(s)")
                        
                        # Summary metrics
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("Successful", len(results))
                        with col2:
                            st.metric("Failed", len(failed_snps))
                        with col3:
                            st.metric("Success Rate", f"{len(results)/len(snp_df)*100:.1f}%")
                        
                        # Results table
                        st.subheader("Primer Design Results")
                        result_df = pd.DataFrame(results)
                        st.dataframe(result_df, use_container_width=True)
                        
                        # Detailed view
                        with st.expander("Detailed Primer Information"):
                            for i, result in enumerate(results):
                                st.write(f"**Primer Pair {i+1}: {result['Seq ID']} - Position {result['SNP Position']}**")
                                
                                col1, col2 = st.columns(2)
                                with col1:
                                    st.write("**Inner Forward (IF)**")
                                    st.code(result['IF Primer'])
                                    st.write(f"Tm: {result['IF Tm']}°C")
                                    st.write(f"GC%: {result['IF GC%']}")
                                    st.write(f"GC Clamp: {result['IF GC Clamp']}")
                                
                                with col2:
                                    st.write("**Inner Reverse (IR)**")
                                    st.code(result['IR Primer'])
                                    st.write(f"Tm: {result['IR Tm']}°C")
                                    st.write(f"GC%: {result['IR GC%']}")
                                    st.write(f"GC Clamp: {result['IR GC Clamp']}")
                                
                                st.divider()
                        
                        # Gel visualization
                        if of_primer and or_primer and gel_data:
                            st.subheader("Gel Electrophoresis Visualization")
                            
                            lane_labels = []
                            band_sizes = []
                            
                            for i, data in enumerate(gel_data):
                                lane_labels.append(f"Lane {i+1}\n{data['seq_id']}")
                                sizes = [size for size in data['amplicon_sizes'].values() if size is not None]
                                band_sizes.append(sizes)
                            
                            fig = create_gel_visualization(band_sizes, lane_labels)
                            st.plotly_chart(fig, use_container_width=True)
                            
                            st.info("**Interpretation**: Outer product appears in all samples. Wild-type and mutant products are allele-specific.")
                        
                        # Download
                        csv = result_df.to_csv(index=False)
                        st.download_button(
                            label="Download Results (CSV)",
                            data=csv,
                            file_name="tetradyme_results.csv",
                            mime="text/csv"
                        )
                    else:
                        st.warning("No suitable primers found for any SNPs.")
                    
                    if failed_snps:
                        with st.expander("Failed SNPs"):
                            for error in failed_snps:
                                st.error(error)
    
    except Exception as e:
        st.error(f"Error: {str(e)}")
else:
    st.warning("Please upload both FASTA and CSV files to begin.")

st.divider()
st.caption("Tetradyme - Tetra-Primer ARMS-PCR Designer")