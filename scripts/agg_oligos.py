import pandas as pd


def read_dataframes(filelist):
    frames = []
    for each_frame in filelist:
        frames.append(pd.read_csv(each_frame, sep='\t'))
    concat_frame = pd.concat(frames)
    # add column names
    # colnames = ['oligo', 'sequence', 'size', 'prep']
    # concat_frame.columns = colnames
    return concat_frame

def main():
    
    frames = snakemake.input
    agg_frame = read_dataframes(frames)
    unique_oligos = agg_frame.drop_duplicates(subset=['context'])
    unique_oligos.to_csv(str(snakemake.output), sep='\t', index=False)

if __name__ == '__main__':
    main()
        