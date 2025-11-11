import primer3

def check_3prime_gc_clamp(primer):
    return primer[-1] in ['G', 'C']

def check_dimer_hairpin(primer):
    hairpin = primer3.bindings.calcHairpin(primer)
    homodimer = primer3.bindings.calcHomodimer(primer)
    return hairpin.dg, homodimer.dg

