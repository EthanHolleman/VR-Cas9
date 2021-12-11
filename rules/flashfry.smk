

rule download_flashfry:
    conda:
        '../envs/java.yml'
    output:
        'software/FlashFry-assembly-1.12.jar'
    shell:'''
    rm -r software
    mkdir software
    cd software
    wget https://github.com/mckennalab/FlashFry/releases/download/1.12/FlashFry-assembly-1.12.jar
    '''


rule build_database:
    conda:
        '../envs/java.yml'
    input:
        ff='software/FlashFry-assembly-1.12.jar',
        target='output/sequences/fasta/{series}/{series_prefix}-{insert_number}.fa'
    output:
        'output/databases/{series}/{series_prefix}-{insert_number}'
    params:
        enzyme='spcas9ngg'
    shell:'''
    mkdir -p tmp
    java -Xmx4g -jar {input.ff} index --tmpLocation ./tmp \
    --database {output} --reference {input.target} \
    --enzyme {params.enzyme}
    '''


rule discover_targets:
    conda:
        '../envs/java.yml'
    input:
        database='output/databases/{series}/{series_prefix}-{insert_number}',
        target='output/sequences/fasta/{series}/{series_prefix}-{insert_number}.fa',
        ff='software/FlashFry-assembly-1.12.jar'
    output:
        'output/targets/{series}/{series_prefix}-{insert_number}.targets'
    shell:'''
    java -Xmx4g -jar {input.ff} discover --database {input.database} \
    --fasta {input.target} --output {output}
    '''


rule score_targets:
    conda:
        '../envs/java.yml'
    input:
        ff='software/FlashFry-assembly-1.12.jar',
        discovered_targets='output/targets/{series}/{series_prefix}-{insert_number}.targets',
        database='output/databases/{series}/{series_prefix}-{insert_number}'
    output:
        'output/scored-targets/{series}/{series_prefix}-{insert_number}.targets.scored'
    shell:'''
     java -Xmx4g -jar {input.ff} score --input {input.discovered_targets} \
     --output {output} --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot,rank \
     --database {input.database}
    '''

rule add_pamless:
    conda:
        '../envs/python.yml'
    input:
       'output/scored-targets/{series}/{series_prefix}-{insert_number}.targets.scored'
    output:
        'output/scored-targets-pamless/{series}/{series_prefix}-{insert_number}.pamless.targets.scored.tsv'
    script:'../scripts/pamless.py'


rule append_NEB_kit_seqs:
    conda:
        '../envs/python.yml'
    input:
        'output/scored-targets-pamless/{series}/{series_prefix}-{insert_number}.pamless.targets.scored.tsv'
    output:
        'output/scored-targets-NEB/{series}/{series_prefix}-{insert_number}.NEB.targets.scored.tsv'
    params:
        config=config
    script:'../scripts/add_NEB_seqs.py'


rule select_targets:
    conda:
        '../envs/python.yml'
    input:
        'output/scored-targets-NEB/{series}/{series_prefix}-{insert_number}.NEB.targets.scored.tsv'
    output:
        'output/selected-scored-NEB/{series}/{series_prefix}-{insert_number}.NEB.selected.targets.scored.tsv'
    params:
        config=config
    script:'../scripts/select_targets.py'


rule label_genbank_with_targets:
    conda:
        '../envs/python.yml'
    input:
        record='output/sequences/constructs/{series}/{series_prefix}-{insert_number}.gb',
        targets='output/selected-scored-NEB/{series}/{series_prefix}-{insert_number}.NEB.selected.targets.scored.tsv'
    output:
        'output/genbank-labeled-with-targets/{series}/{series_prefix}-{insert_number}.label.targets.gb'
    script:'../scripts/label_targets.py'


rule draw_plasmids_with_targets:
    conda:
        '../envs/bioDraw.yml'
    input:
        record='output/genbank-labeled-with-targets/{series}/{series_prefix}-{insert_number}.label.targets.gb',
        relabel_dict='label-scheme.yml'
    output:
        'output/plasmid-maps-with-targets/{series}/{series_prefix}-{insert_number}.label.targets.png'
    script:'../scripts/plot_targets.py'







