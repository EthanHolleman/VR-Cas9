# Add a guide RNA column to flashFry output. FlashFry only lists the target
# sites which include PAMS. To convert to gRNA sequences need to remove the
# PAM site (last 3 nucleotides)


import pandas as pd
from Bio.Seq import Seq

def read_df(filepath):
    
    return pd.read_csv(str(filepath), sep='\t')


def add_target_no_pam_column(df):
    
    df['target_no_pam'] = [l[:-3] for l in list(df['target'])]
    df['end_no_pam'] = df['stop'] - 3
    return df

def add_gRNA_column(df):
    no_pams_rc = [str(Seq(s).reverse_complement()) for s in list(df['target_no_pam'])]
    df['gRNA_as_DNA'] = no_pams_rc
    return df
    

def main():
    
    score_file = str(snakemake.input)
    df = read_df(score_file)
    df_no_pam = add_target_no_pam_column(df)
    df_no_pam.dropna(inplace=True)
    df_no_pam_grna = add_gRNA_column(df)
    df.to_csv(str(snakemake.output), sep='\t', index=False)


if __name__ == '__main__':
    main()
    