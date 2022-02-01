
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


rule select_targets:
    conda:
        '../envs/python.yml'
    input:
        scored='output/{run}/scoredTargetsPamlessNEB/{seq_name}.targets.scored.pamless.NEB.tsv',
        fasta='output/{run}/reslice/{seq_name}.reslice.fa'
    output:
        selected='output/{run}/selectedScoredTargetsPamlessNEB/{seq_name}.targets.scored.pamless.NEB.selected.tsv',
        all_labeled='output/{run}/selectedScoredTargetsPamlessNEB/{seq_name}.targets.scored.pamless.NEB.label.tsv'
    script:'../scripts/select_targets.py'


rule plot_targets:
    conda:
        '../envs/R.yml'
    input:
        'output/{run}/selectedScoredTargetsPamlessNEB/{seq_name}.targets.scored.pamless.NEB.label.tsv'
    output:
        'output/{run}/plots/{seq_name}.targets.scored.pamless.NEB.label.png'
    script:'../scripts/plotOnTargetScores.R'