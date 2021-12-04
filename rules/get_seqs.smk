

rule download_sliced_sequences:
    # Download R-loop prediction substrate sequences. These are VR-insert
    # sequences that have been sliced relative to the TSS.
    output:
        'output/sequences/prediction_substrates_2021-11-14.fa'
    shell:'''
    cd output/sequences
    wget https://raw.githubusercontent.com/EthanHolleman/Quon-R-loop-prediction-substrates/main/prediction_substrates_2021-11-14.fa
    '''


rule split_select_init_seqs:
    conda:
        '../envs/python.yml'
    input:
        'output/sequences/prediction_substrates_2021-11-14.fa'
    output:
        expand(
            'output/sequences/split/VR-{num_inserts}.fa', 
            num_inserts=list(range(1, NUM_INSERTS+1))
        )
    params:
        VR_regex='VR-\d+',
        output_dir='output/sequences/split'
    shell:'''
    python scripts/split_fasta.py {input} {params.output_dir} --rename "{params.VR_regex}" \
    --filter "init"
    '''



