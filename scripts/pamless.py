# Add pamless target site information to FlashFry output

import pandas as pd
from Bio.Seq import Seq

def read_df(filepath):
    
    return pd.read_csv(str(filepath), sep='\t')


def add_target_no_pam_column(df):
    
    pamless_columns = []
    for index, row in df.iterrows():
        
        # FlashFry outputs all sequences in the forward orientation even
        # if the target site is on the reverse strand. So the PAM site will
        # always be at the end of the sequence string.
        row = {
            'target_no_pam': each_target[:-3],
            'start_no_pam': row['start'],
            'end_no_pam': row['end'] - 3
        }
        pamless_columns.append(row)
    
    # "cbind" the pamless columns to the old dataframe
    # https://stackoverflow.com/questions/33088010/pandas-column-bind-cbind-two-data-frames
    df_c = pd.concat(
        [df.reset_index(drop=True), pd.Dataframe(pamless_columns)], axis=1
        )
    return df_c
        


def main():
    
    score_file = str(snakemake.input)
    df = read_df(score_file)
    df_no_pam = add_target_no_pam_column(df)
    df_no_pam.dropna(inplace=True)
    df.to_csv(str(snakemake.output), sep='\t', index=False)


if __name__ == '__main__':
    main()
    