# Create tsv file that can be used for bulk IDT orders

import pandas as pd
from datetime import datetime


def make_oligo_names(name_prefix, number_oligos, start_index=0):
    def make_oligo_name(name_prefix, oligo_number):
        date_str = datetime.now().strftime("%m-%d-%Y")
        return f"{name_prefix}_{oligo_number}-{date_str}"

    for i in range(start_index, number_oligos+start_index):
        yield make_oligo_name(name_prefix, i)


def main():

    targets = pd.read_csv(str(snakemake.input), sep="\t")
    
    num_targets = len(targets)
    name_prefix = snakemake.params["name_prefix"]
    oligo_names = list(make_oligo_names(
        name_prefix, num_targets, snakemake.params['oligo_ID_start'])
                       )
    
    targets["Name"] = oligo_names
    targets["Concentration"] = "25nm"  # IDT constant for primer concentration
    targets["Desalination"] = "STD"  # IDT standard desalination
    idt_df = targets[["Name", "NEB_EnGen_Oligos", "Concentration", "Desalination"]]
    print(idt_df)
    idt_df.to_csv(str(snakemake.output), sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
