
NUM_INSERTS = 31


include: 'rules/flashfry.smk'
include: 'rules/get_seqs.smk'




rule all:
    input:
        expand(
            'output/plasmid-maps-with-targets/{series}/{series_prefix}-{insert_number}.label.targets.png',
            insert_number=list(range(1, NUM_INSERTS+1)), series=[config['SERIES_NAME']],
            series_prefix=config['SERIES_PREFIX']
            )
        