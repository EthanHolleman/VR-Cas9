import pandas as pd
from Bio import SeqIO
import numpy as np
from Bio import SeqIO


def remove_excluded_targets(df, excluded):
    # https://www.kite.com/python/answers/how-to-get-rows-from-a-dataframe-that-are-not-in-another-dataframe-in-python
    common = df.merge(excluded, on=["target"])
    return df[~df.target.isin(common.target)]


def filter_dangerous(df):
    """Remove sgRNA's that flashfry deems to have dangerous GC of polyT
    content.

    Args:
        df (DataFrame): Pandas df.
    """
    print(df)
    df = df.loc[df["dangerous_GC"] == "NONE"]
    df = df.loc[df["dangerous_polyT"] == "NONE"]
    return df


def filter_high_off_target(df, min_score=0.95):
    """Remove guides with predicted off targets via the 
     Doench et al. Nature Biotechnology, 2016 specificity score. Higher
     scores indicate lower off target potential. Default min value os 0.95
     as off targets should be very low for plasmid templates. 

    Args:
        df ([DataFrame]): Pandas df of flashfry produced metrics.
    """
    return df.loc[df["DoenchCFD_specificityscore"] > min_score]


def filter_strand(df, strand):

    return df.loc[df["orientation"] == strand]


def select_targets(df, num_targets, force_zero, seq_len, end_offset=0):
    """Select final Cas9 target sites, return as a pandas DataFrame.

    Args:
        df (DataFrame): Input df, should already be filtered.
        num_targets (int): Number of target sites to select.
        force_zero (bool): Force selection of closest site to start of sequence.
        seq_len (int): Length of template sequence.
        end_offset (int, optional): Length to subtract from total length. Defaults to 0.

    Returns:
        [dataFrame]: Dataframe of selected target sequences.
    """

    targets = []
    if force_zero:
        print(df)
        targets.append(df.sort_values("start").iloc[0])

    spread = int((seq_len) / num_targets)
    locs = list(np.flip(np.arange(0, seq_len - end_offset, spread)))  # flip for popping
    print(locs, seq_len - end_offset, spread)

    def select_target(search_start, search_end):
        in_region = df.loc[(df["start"] >= search_start) & (df["stop"] <= search_end)]
        if len(in_region) > 0:
            in_region_by_doench = in_region.sort_values(
                "Doench2014OnTarget", ascending=False
            )
            return in_region_by_doench.iloc[0]  # best on target score in defined region

    while locs:
        cur_loc = locs.pop()
        print(cur_loc)
        search_start, search_end = cur_loc - (spread / 1.5), cur_loc + (spread / 1.5)
        if search_start < 0:
            search_start = 0
        if search_end > seq_len:
            search_end = seq_len
        targets.append(select_target(search_start, search_end))
        targets = [t for t in targets if isinstance(t, pd.Series)]

    return pd.DataFrame(targets)


def label_select_targets(scored, selected):

    scored["selected"] = False
    scored.loc[scored.target.isin(selected.target), "selected"] = True
    return scored


def main():

    excluded_targets = pd.read_csv(snakemake.params['excluded'], sep='\t')
    scored_targets = pd.read_csv(snakemake.input["scored"], sep="\t")
    
    remove_excluded_targets(scored_targets, excluded_targets)
    
    
    # drop dangerous
    scored_targets = filter_dangerous(scored_targets)
    scored_targets = filter_high_off_target(scored_targets, snakemake.config["MINOFF"])
    scored_targets = filter_strand(scored_targets, snakemake.config["STRAND"])

    seq_len = len(SeqIO.read(snakemake.input["fasta"], format="fasta"))

    targets = select_targets(
        scored_targets,
        snakemake.config["NUM_TARGETS"],
        snakemake.config["ZERO_TARGET"],
        seq_len,
        snakemake.config["END_OFFSET"],
    )

    all_scored_targets_label = label_select_targets(scored_targets, targets)

    targets.to_csv(str(snakemake.output["selected"]), index=False, sep="\t")
    all_scored_targets_label.to_csv(
        snakemake.output["all_labeled"], index=False, sep="\t"
    )


if __name__ == "__main__":
    main()
