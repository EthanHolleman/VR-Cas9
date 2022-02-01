rule make_IDT_order_tables:
    conda:
        '../envs/python.yml'
    input:
        'output/{run}/selectedScoredTargetsPamlessNEB/{seq_name}.targets.scored.pamless.NEB.selected.tsv'
    output:
        'output/{run}/IDTOrder/{seq_name}.oligos.tsv'
    params:
        name_prefix=lambda wildcards: f'{wildcards.seq_name}',
        oligo_ID_start=config['OLIGO_ID_START']
    script:'../scripts/IDT_order.py'