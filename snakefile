
NUM_INSERTS = 31


include: 'rules/flashfry.smk'
include: 'rules/get_seqs.smk'


rule all:
    input:
        expand(
            'output/selected-scored-NEB/VR-{num_inserts}.NEB.selected.targets.scored.tsv',
            num_inserts=list(range(1, NUM_INSERTS+1))
        )
