
import pandas as pd
from Bio import SeqIO



def get_ideal_points(min_dist, max_dist, number_targets):
    assert max_dist > min_dist
    usable_length = max_dist - min_dist
    step = int(usable_length / number_targets)
    ideal_points = []
    return [i for i in range(min_dist, max_dist, step)]



def filter_targets_by_strand(df, target_strand, nick_target=False):
# nicking occurs on strand opposite target sequence
    print(target_strand)
    if nick_target:
        return df.loc[df['orientation'] == target_strand]
    else:
        if target_strand == 'FWD':
            return df.loc[df['orientation'] == 'RVS']
        elif target_strand == 'RVS':
            return df.loc[df['orientation'] == 'FWD']
        else:
            raise TypeError('Target strand must be FWD or REV')
    

def target_picker(df, ideal_point, filter_danger=True):
    # pick closest targets to each ideal point

    if filter_danger:
        df = df.loc[(df['dangerous_GC'] == 'NONE') & (df['dangerous_polyT'] == 'NONE')]
    
    return df.iloc[(df['start']-ideal_point).abs().argsort()[:1]]


def pick_targets(df, min_dist, max_dist, number_targets, target_strand, nicked_strand):
    
    ideal_points = get_ideal_points(min_dist, max_dist, number_targets)
    df_stranded = filter_targets_by_strand(df, target_strand, nicked_strand)
    target_rows = []
    for each_point in ideal_points:
        target_rows.append(target_picker(df_stranded, each_point))
    
    return pd.concat(target_rows)


def main():
    
    config = snakemake.params['config']
    all_targets = pd.read_csv(str(snakemake.input), sep='\t')
    targets = pick_targets(
        all_targets, config['MIN_DIST'], config['MAX_DIST'], config['NUM_TARGETS'],
        config['TARGET_STRAND'], config['NICKED_STRAND_IS_TARGET']
    )
    targets.to_csv(str(snakemake.output), sep='\t', index=False)


if __name__ == '__main__':
    main()


    

    
    
        
