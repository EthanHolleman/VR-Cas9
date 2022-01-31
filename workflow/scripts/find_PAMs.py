# Python script to identify potential guide RNAs for CRISPR-Cas9
# from a collection of sequence records in fasta format.

import re
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

PAM_RE = r".GG"  #  Streptococcus pyogenes PAM as default


def read_record(filepath):
    return SeqIO.parse(str(filepath), "fasta")


def find_PAMs_in_record(record, PAM_RE=PAM_RE):
    """Return a list of integers representing the start indexes of all
    PAM sites in a sequence record.

    Args:
        record (SeqRecord): SeqRecord to search
    """
    return [(m.start(0), m.end(0)) for m in re.finditer(PAM_RE, str(record.seq).upper())]


def make_gRNAs(record, PAM_RE=PAM_RE, gRNA_len=20):
    """Make list of dictionaries describing gRNAs for a given input DNA
    sequence.

    Args:
        record (SeqRecord): SeqRecord to make gRNAs for.
        PAM_RE (Regex, optional): Regex expression to find PAM sites. Defaults to PAM_RE.
        gRNA_len (int, optional): Length of gRNA. Defaults to 20.

    Returns:
        list: List of dictionaries, each describing one gRNA.
    """

    PAM_sites = find_PAMs_in_record(record)
    gRNAs = []
    for i, each_PAM in enumerate(PAM_sites):
        # TODO: Implement guide RNA for circular sequences (wrap around)
        # for now do not make guide RNAs that would cause index errors

        start = each_PAM[0] - gRNA_len
        if each_PAM[0] - gRNA_len >= 0:
            gRNA_seq = record[start : each_PAM[0]]
            # gRNA = SeqRecord(gRNA_seq, id=gRNA_seq.id, name=f'gRNA-{i+1}')
            gRNAs.append(
                {
                    "gRNA": str(gRNA_seq.seq),
                    "target_seq_id": record.id,
                    "target_seq_descrip": record.description,
                    "gRNA_num": i + 1,
                    "start": start,
                    "PAM_start": each_PAM[0],
                    "PAM_end": each_PAM[1],
                }
            )
    return gRNAs


def get_args():

    parser = ArgumentParser(
        "Design guide CRISPER-Cas9 guide RNAs for a collections \
                            of fasta formated sequences."
    )
    parser.add_argument("F", help="Input fasta file.")
    parser.add_argument(
        "--out",
        help="Output path to write gRNAs to. Default is  \
                                stout.",
        default=None,
    )
    parser.add_argument(
        "--pam",
        help='PAM site regex. Defaults to the Streptococcus \
                                        pyogenes PAM, NGG. Or as regex ".GG" \
                                        all regexes should be for capitalized \
                                        sequences.',
        default=PAM_RE,
    )
    parser.add_argument(
        "--len", help="Guide RNA length, defaults to 20 bp", default=20, type=int
    )
    return parser.parse_args()


def main():

    args = get_args()
    records = read_record(args.F)
    gRNA_dicts_list = []
    for r in records:
        gRNA_dicts_list += make_gRNAs(r, args.pam, args.len)

    df = pd.DataFrame(gRNA_dicts_list)
    if args.out:
        df.to_csv(args.out, sep="\t", index=False)
    else:
        print(df.to_csv(sep="\t", index=False))


if __name__ == "__main__":
    main()
