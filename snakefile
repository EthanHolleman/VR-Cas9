
NUM_INSERTS = 31


include: 'rules/flashfry.smk'
include: 'rules/get_seqs.smk'


rule all:
    input:
        expand(
            'output/scored-targets-gRNA/VR-{num_inserts}.gRNA.targets.scored',
            num_inserts=list(range(1, NUM_INSERTS+1))
        )
