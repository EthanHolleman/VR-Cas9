from Bio import SeqIO


def main():

    record = SeqIO.read(str(snakemake.input), "gb")
    SeqIO.write(record, str(snakemake.output), "fasta")


if __name__ == "__main__":
    main()
