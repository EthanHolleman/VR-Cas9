rule reslice_genbank:
    # Slice the sequence of a circular genbank file to a linear fasta
    # file starting at a user defined feature.
    conda:
        '../envs/python.yml'
    output:
        'output/{run}/reslice/{seq_name}.reslice.fa'
    params:
        genbank=lambda wildcards: samples.loc[samples['seq_name'] == wildcards.seq_name]['genbank_path'].values[0],
        start_feature=lambda wildcards: samples.loc[samples['seq_name'] == wildcards.seq_name]['start_feature'].values[0]
    script:'../scripts/reslice_gb.py'


