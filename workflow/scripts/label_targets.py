from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pandas as pd


def add_cas9_targets_as_features(df, record):
    def targets_to_features(df):
        records = df.to_dict(orient="records")
        for i, each_row in enumerate(records):
            if each_row["orientation"] == "FWD":
                strand = 1
            else:
                strand = -1
            yield SeqFeature(
                FeatureLocation(each_row["start"], each_row["stop"], strand=strand),
                id=f"Cas9 Target {i+1}",
                type=f"Cas9 Target {i+1}",
                strand=strand  # snapgene uses type as
                # the displayed
                # feature name
            )

    record.features += list(targets_to_features(df))
    return record


def main():

    record = SeqIO.read(str(snakemake.input["record"]), "gb")
    targets_df = pd.read_csv(str(snakemake.input["targets"]), sep="\t")
    feature_record = add_cas9_targets_as_features(targets_df, record)
    SeqIO.write(feature_record, str(snakemake.output), "gb")


if __name__ == "__main__":
    main()
