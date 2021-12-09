# Add pamless target site information to FlashFry output

import pandas as pd
from Bio.Seq import Seq

def read_df(filepath):
    
    return pd.read_csv(str(filepath), sep='\t')


def add_target_no_pam_column(df):
    
    pamless_columns = []
    for each_target, strand in zip(list(df['targets']), list(df['orientation'])):
        if strand == 'FWD':  # pam is last three nucs
            row = {
                'target_no_pam': each_target[:-3],
                'start_no_pam': '',
                'end_no_pam': ''
            }
        elif strand == 'REV':
            pass
        else:
            raise TypeError('Orientation must be FWD or REV')
        
        
    df['target_no_pam'] = [l[:-3] for l in list(df['target'])]
    df['end_no_pam'] = df['stop'] - 3
    return df


def main():
    
    score_file = str(snakemake.input)
    df = read_df(score_file)
    df_no_pam = add_target_no_pam_column(df)
    df_no_pam.dropna(inplace=True)
    df.to_csv(str(snakemake.output), sep='\t', index=False)


if __name__ == '__main__':
    main()
    