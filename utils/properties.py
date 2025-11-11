from Bio.SeqUtils import gc_fraction, molecular_weight

def calculate_gc(seq):
    return round(gc_fraction(seq) * 100, 2)

def calculate_molecular_weight(seq):
    return round(molecular_weight(seq, 'DNA'), 2)

def extinction_coefficient(seq):
    # Approximate values
    return {
        'A': 15400,
        'C': 7400,
        'G': 11500,
        'T': 8700
    }.get(seq[0], 10000) * len(seq)  # crude sum

def od260_to_ug_nmole(seq):
    mw = molecular_weight(seq, 'DNA')
    return round(33 / mw * 1e6, 2), round(33 / 330, 2)  # µg/OD260 and nmol/OD260

