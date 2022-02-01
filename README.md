# VR-Cas9

![Workflow](https://github.com/EthanHolleman/VR-Cas9/actions/workflows/main.yml/badge.svg)

Identify potential CRISPR-Cas9 targets in plasmid sequences using [FlashFry](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-018-0545-0), select `n` approximately uniformly distributed
target sites along the plasmid sequence, design oligos for sgRNA
synthesis using [NEB enGen sgRNA synthesis kit](https://www.neb.com/products/e3322-engen-sgrna-synthesis-kit-s-pyogenes#Product%20Information) and produce an [IDT](https://www.idtdna.com/pages) oligo order sheet. 

## Dependencies

- [Conda](https://www.anaconda.com/products/individual)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)
- [Pandas](https://anaconda.org/conda-forge/pandas)

## Run

The workflow can be run most easily with `snakemake -j 1 --configfile config/config.yml --use-conda`.
from the `workflow` directory. To run the workflow on a cluster reconfigure `cluster.yml` and `run.sh` to work with your system (currently configured for SLURM (UC Davis CRICK)).

## Configuration

In order to run the workflow with your own plasmid sequences
you will need to modify a few configuration files described
below.

### `workflow/config/config.yml`

This is the snakemake config file. The parameters are commented
and explained within the file itself but those with additional
detail are described below.

#### `sample_tsv`

This should be a filepath relative to where the workflow
is executed from to a tsv file containing the following
headers.

##### `seq_name`

Unique ID for the template sequence.

##### `genbank_path`

Path to genbank file containing sequence and annotations for this template. If relative path is used should be relative to workflow execution directory.

##### `start_feature`

Plasmids are circular sequences but are usually described linearly for convenience of assigning exact coordinates to annotations. Often the arbitrary start of a plasmid within a genbank file is not the zero index we want to design sgRNAs relative to. Set this field to the “\label” field of a feature within the genbank file to use that feature as the start of the linearized plasmid.

##### `force_zero_guide`

Normally the workflow tries to evenly distribute a number of target sites throughout the provided plasmid sequence. It may be desirable to place a target site as close to the start of a sequence as possible. Setting this field to “True” forces the workflow to always include the closest target site to the start of the sequence.

##### `total_guides`

Total number of sgRNAs to select. Workflow will attempt to evenly distribute target sites along plasmid sequence.

##### `strand`

Select target sites on the + “FWD” or – “RVS” strand only.

##### `exclude_start`

Start of region of plasmid (bp index) to exclude targets which overlap. This
is relative to the feature used in `start_feature` field not the supplied
genbank file.

##### `exclude_end`

End of region of plasmid (bp index) to exclude targets which overlap. This
is relative to the feature used in `start_feature` field not the supplied
genbank file.

If you wish to include all of the plasmid for the target search set
`exclude_start` and `exclude_end` to `0`.


#### `exluded_targets_tsv`

Path to a tsv file that described target sequences (including PAM)
that should not be included in the output. Useful for ensuring
oligos with specific targets that have already been ordered do
not show up in the output again. Should contain the following
fields.


##### `oligo_name`

Unique identifier for the target sequence.

##### `target`

Target sequence to exclude from output.









