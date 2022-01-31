rule download_chopchop:
    output:
        'software/chopchop/README.md'
    shell:'''
    cd software
    git clone https://bitbucket.org/valenlab/chopchop.git
    '''

rule build_bowtie_index:
    conda:
        '../envs/chopchop.yml'
    input:
        fasta='output/sequences/fasta/{series}/{series_prefix}-{insert_number}.fa'
    output:
        directory('output/bowtie2/databases/{series}/{series_prefix}-{insert_number}')
    shell:'''
    bowtie-build {input.fasta} {output}
    '''




