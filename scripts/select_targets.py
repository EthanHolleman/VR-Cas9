import pandas as pd
from Bio import SeqIO
import numpy as np


def get_ideal_points(min_dist, max_dist, number_targets, add_zero_point=True):

    assert max_dist > min_dist
    usable_length = max_dist - min_dist
    step = int(usable_length / number_targets)
    points = list(np.arange(min_dist, max_dist, step))
    if add_zero_point:  # force point at start of sequence
        points = [0] + points
    return points


def filter_targets_by_strand(df, target_strand, nick_target=False):
    # nicking occurs on strand opposite target sequence
    print(target_strand)
    if nick_target:
        return df.loc[df["orientation"] == target_strand]
    else:
        if target_strand == "FWD":
            return df.loc[df["orientation"] == "RVS"]
        elif target_strand == "RVS":
            return df.loc[df["orientation"] == "FWD"]
        else:
            raise TypeError("Target strand must be FWD or REV")


def target_picker(
    targets_df,
    ideal_point,
    filter_danger=True,
    max_dist_from_ideal=50,
    look_step=10,
    min_doench_score=0.2,
):
    # pick closest targets to each ideal point
    df = targets_df.copy()
    if filter_danger:
        df = df.loc[(df["dangerous_GC"] == "NONE") & (df["dangerous_polyT"] == "NONE")]

    # Get distance to desired target location
    df["target_dist"] = df["start"] - ideal_point
    df["target_dist"] = abs(df["target_dist"])

    # Select targets within `max_dist_from_ideal` distance plus add some
    # tolerance based on the location of the ideal point. Further point is
    # away from 0 (assumed TSS) the more we do not care as much about
    # its exact positioning.
    
    df = df.loc[df["target_dist"] <= max_dist_from_ideal + ideal_point/100]
    # Sort targets by distance to target and by doench2014ontarget score
    # (Doench et al. Nature Biotechnology, 2014 ). Higher Doench scores
    # indicate higher on-target efficiency
    df = df.sort_values(["Doench2014OnTarget"], ascending=[False])
    df = df.loc[df["Doench2014OnTarget"] >= min_doench_score]
    
    if len(df) == 0:
        return []
    return df.iloc[0]


def pick_targets(df, min_dist, max_dist, number_targets, target_strand, nicked_strand):

    ideal_points = get_ideal_points(min_dist, max_dist, number_targets)
    df_stranded = filter_targets_by_strand(df, target_strand, nicked_strand)
    target_rows = []

    for each_point in ideal_points:
        target = target_picker(df_stranded, each_point)
        if len(target) > 0:
            target_rows.append(target)

    targets = pd.DataFrame(target_rows)
    targets = targets.drop_duplicates()
    return targets


def add_plasmid_identifiers(target_df, id_dict):
    for key, value in id_dict.items():
        target_df[key] = value
    return target_df
    

def main():

    config = snakemake.params["config"]
    id_dict = snakemake.params['id_dict']
    
    all_targets = pd.read_csv(str(snakemake.input), sep="\t")
    
    targets = pick_targets(
        all_targets,
        config["MIN_DIST"],
        config["MAX_DIST"],
        config["NUM_TARGETS"],
        config["TARGET_STRAND"],
        config["NICKED_STRAND_IS_TARGET"],
    )
    targets = add_plasmid_identifiers(targets, id_dict)
    targets.to_csv(str(snakemake.output), sep="\t", index=False)


if __name__ == "__main__":
    main()
