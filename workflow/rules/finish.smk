rule tar_output:
    input:
        expand(
            'output/{run}/IDTOrder/{seq_name}.oligos.tsv',
            run=RUN_NAME, seq_name=list(samples['seq_name'])
        ),
        expand(
            'output/{run}/plots/{seq_name}.targets.scored.pamless.NEB.label.png',
            run=RUN_NAME, seq_name=list(samples['seq_name'])
        )
    output:
       'output/tars/{run}.tar.gz'
    params:
        output_dir = f'output/{RUN_NAME}'
    shell:'''
    tar -zcvf {output} {params.output_dir}
    '''

