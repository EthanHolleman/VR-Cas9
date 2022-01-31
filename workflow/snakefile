
NUM_INSERTS = 31


include: 'rules/flashfry.smk'
include: 'rules/get_seqs.smk'

if config['RUN_BACKUP'] == True:
    backup = ['.backup.done']
else:
    backup = []

rule all:
    input:
        expand(
            'output/plasmid-maps-with-targets/{series}/{series_prefix}-{insert_number}.label.targets.png',
            insert_number=list(range(1, NUM_INSERTS+1)), series=[config['SERIES_NAME']],
            series_prefix=config['SERIES_PREFIX']
            ),
        expand(
            'output/IDT-order/{series}/{series_prefix}-{insert_number}.NEB.selected.targets.scored.IDT.tsv',
            insert_number=list(range(1, NUM_INSERTS+1)), series=[config['SERIES_NAME']],
            series_prefix=config['SERIES_PREFIX']
        ),
        backup
        