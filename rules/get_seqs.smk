

rule download_sliced_sequences:
    output:
        expand(
            'output/sequences/constructs/{series}/{series_prefix}-{insert_number}.gb',
            insert_number=list(range(1, NUM_INSERTS+1)), series=['T7_initiation_series'],
            series_prefix=['T7_init_VR']
        )

    shell:'''
    rm -r output/sequences
    mkdir output/sequences
    cd output/sequences
    wget https://github.com/EthanHolleman/plasmid-VR-design/releases/download/v2.0/constructs.zip
    unzip -o constructs.zip
    '''

rule make_fasta_files:
    conda:
        '../envs/python.yml'
    input:
        'output/sequences/constructs/{series}/{series_prefix}-{insert_number}.gb'
    output:
        'output/sequences/fasta/{series}/{series_prefix}-{insert_number}.fa'
    script:'../scripts/genbank2fasta.py'


# rule split_select_init_seqs:
#     conda:
#         '../envs/python.yml'
#     input:
#         'output/sequences/prediction_substrates_2021-11-14.fa'
#     output:
#         expand(
#             'output/sequences/split/VR-{num_inserts}.fa', 
#             num_inserts=list(range(1, NUM_INSERTS+1))
#         )
#     params:
#         VR_regex='VR-\d+',
#         output_dir='output/sequences/split'
#     shell:'''
#     python scripts/split_fasta.py {input} {params.output_dir} --rename "{params.VR_regex}" \
#     --filter "init"
#     '''



