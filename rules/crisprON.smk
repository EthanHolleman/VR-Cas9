rule download_crisprON:
    output:
        'output/software/crispron/README.md'
    shell:'''
    mkdir -p output/software
    cd output/software
    git clone https://github.com/RTH-tools/crispron.git
    '''

rule download_crisprOFF:
    output:
        'output/software/crisproff-1.1.2/README.md'
    shell:'''
    mkdir -p output/software
    cd output/software
    wget https://rth.dk/resources/crispr/crisproff/downloads/crisproff-1.1.2.tar.gz
    tar -xf crisproff-1.1.2.tar.gz
    '''


rule copy_crisprOFF_files:
    input:
        crisproff='output/software/crisproff-1.1.2/README.md',
        crispron='output/software/crispron/README.md'
    output:
        pipeline='output/software/crispron/bin/CRISPRspec_CRISPRoff_pipeline.py',
        energy='output/software/crispron/data/model/energy_dics.pkl'
    params:
        crisproff='output/software/crisproff-1.1.2',
        crispron='output/software/crispron'
    shell:'''
    cp {params.crisproff}/CRISPRspec_CRISPRoff_pipeline.py "{output.pipeline}"
    cp {params.crisproff}/energy_dics.pkl "{output.energy}"
    '''

rule CRISPRON:
    conda:
        '../envs/crisprON.yml'
    input:
        fasta='output/sequences/fasta/{series}/{series_prefix}-{insert_number}.fa',
        crispron='output/software/crispron/README.md',
        energy='output/software/crispron/data/model/energy_dics.pkl'
    output:
        directory('output/crisprON/{series}/{series_prefix}-{insert_number}-predictions')
    shell:'''
    output/software/crispron/bin/CRISPRon.sh {input.fasta} {output}
    '''

rule CRISPRON_ALL:
    input:
        expand(
            'output/crisprON/{series}/{series_prefix}-{insert_number}-predictions',
            insert_number=list(range(1, NUM_INSERTS+1)), series=[config['SERIES_NAME']],
            series_prefix=config['SERIES_PREFIX']
        )
    output:
        'output/crisprON/done.txt'
    shell:'''
    touch {output}
    '''



