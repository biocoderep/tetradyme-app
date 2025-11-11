from Bio.SeqUtils import MeltingTemp as mt

def calculate_tm(seq, Na=50, Mg=0, dNTPs=0):
    return mt.Tm_NN(seq, nn_table=mt.DNA_NN4, Na=Na, Mg=Mg, dNTPs=dNTPs)

