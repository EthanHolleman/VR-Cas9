

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
        target='output/sequences/split/VR-{num_inserts}.fa'
    output:
        'output/databases/VR-{num_inserts}'
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
        database='output/databases/VR-{num_inserts}',
        target='output/sequences/split/VR-{num_inserts}.fa',
        ff='software/FlashFry-assembly-1.12.jar'
    output:
        'output/targets/VR-{num_inserts}.targets'
    shell:'''
    java -Xmx4g -jar {input.ff} discover --database {input.database} \
    --fasta {input.target} --output {output}
    '''


rule score_targets:
    conda:
        '../envs/java.yml'
    input:
        ff='software/FlashFry-assembly-1.12.jar',
        discovered_targets='output/targets/VR-{num_inserts}.targets',
        database='output/databases/VR-{num_inserts}'
    output:
        'output/scored-targets/VR-{num_inserts}.targets.scored'
    shell:'''
     java -Xmx4g -jar {input.ff} score --input {input.discovered_targets} \
     --output {output} --scoringMetrics doench2014ontarget,doench2016cfd,dangerous,hsu2013,minot \
     --database {input.database}
    '''


