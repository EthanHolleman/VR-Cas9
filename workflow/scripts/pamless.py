# Add pamless target site information to FlashFry output

import pandas as pd
from Bio.Seq import Seq


def read_df(filepath):

    return pd.read_csv(str(filepath), sep="\t")


def add_target_no_pam_column(df):

    pamless = []
    for index, row in df.iterrows():

        # FlashFry outputs all sequences in the forward orientation even
        # if the target site is on the reverse strand. So the PAM site will
        # always be at the end of the sequence string.
        pamless.append(row["target"][:-3])
    
    df['target_no_pam'] = pamless
    return df



def main():

    score_file = str(snakemake.input)
    df = read_df(score_file)
    df_no_pam = add_target_no_pam_column(df)
    print(df_no_pam)
    
    assert "target_no_pam" in df_no_pam.columns
    assert len(df) == len(df_no_pam)
    
    df_no_pam.to_csv(str(snakemake.output), sep="\t", index=False)


if __name__ == "__main__":
    main()
