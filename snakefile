
NUM_INSERTS = 31


include: 'rules/flashfry.smk'
include: 'rules/get_seqs.smk'


rule all:
    input:
        expand(
            'output/scored-targets-NEB/VR-{num_inserts}.NEB.targets.scored.tsv',
            num_inserts=list(range(1, NUM_INSERTS+1))
        )
