# Add sequences required for NEB EnGen sgRNA synthesis kit
# https://www.neb.com/protocols/2016/05/11/engen-sqrna-synthesis-kit-s-pyogenes-protocol-e3322

import pandas as pd

def check_for_5_prime_G(sequence):
    '''Add G nucleotide to 5' end of sequence if not present. Requirement
    of NEB kit.

    Args:
        sequence (str): Target sequence (no PAM site).
    '''
    if sequence.upper() != 'G':
        sequence = 'G' + sequence
    
    return sequence


def add_NEB_seqs(T7_seq, scaffold, target):
    return T7_seq + target + scaffold


def format_target_seq_for_NEB_kit(T7_seq, scaffold, target_seq):
    target_seq = check_for_5_prime_G(target_seq)
    return add_NEB_seqs(T7_seq, scaffold, target_seq)


def main():
    
    targets_df = pd.read_csv(str(snakemake.input), sep='\t')
    config = snakemake.params['config']
    print(config)
    
    T7_seq, scaffold = config['T7_PROMOTER'], config['THREE_PRIME_OVERLAP']
    targets_no_PAM = list(targets_df['target_no_pam'])
    print(targets_no_PAM)
    NEB_oligos = [
        format_target_seq_for_NEB_kit(T7_seq, scaffold, target_seq)
        for target_seq in targets_no_PAM]
    targets_df['NEB_EnGen_Oligos'] = NEB_oligos
    
    targets_df.to_csv(str(snakemake.output), sep='\t', index=False)

if __name__ == '__main__':
    main()
    