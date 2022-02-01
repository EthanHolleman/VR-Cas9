# VR-Cas9

![Workflow](https://github.com/EthanHolleman/VR-Cas9/actions/workflows/main.yml/badge.svg)

Identify potential CRISPR-Cas9 targets in plasmid sequences using [FlashFry](https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-018-0545-0).

## Dependencies

- [Conda](https://www.anaconda.com/products/individual)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)

## Run

The workflow can be run most easily with `snakemake -j 1 --configfile config/config.yml --use-conda`.
from the `workflow` directory. To run the workflow on a cluster reconfigure `cluster.yml` and `run.sh` to work with your system (currently configured for SLURM (UC Davis CRICK)).


