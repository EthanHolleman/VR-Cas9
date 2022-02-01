
rule add_pamless:
    conda:
        '../envs/python.yml'
    input:
       'output/{run}/flashFryScores/{seq_name}.targets.scored'
    output:
        'output/{run}/scoredTargetsPamless/{seq_name}.targets.scored.pamless.tsv'
    script:'../scripts/pamless.py'


rule append_NEB_kit_seqs:
    conda:
        '../envs/python.yml'
    input:
        'output/{run}/scoredTargetsPamless/{seq_name}.targets.scored.pamless.tsv'
    output:
        'output/{run}/scoredTargetsPamlessNEB/{seq_name}.targets.scored.pamless.NEB.tsv'
    params:
        config=config
    script:'../scripts/add_NEB_seqs.py'


# rule select_targets:
#     conda:
#         '../envs/python.yml'
#     input:
#         'output/scored-targets-NEB/{series}/{series_prefix}-{insert_number}.NEB.targets.scored.tsv'
#     output:
#         'output/selected-scored-NEB/{series}/{series_prefix}-{insert_number}.NEB.selected.targets.scored.tsv'
#     params:
#         config=config
#     script:'../scripts/select_targets.py'