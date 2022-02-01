from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def reslice_genbank(gb_path, start_feature_label, sliced_output_path):
    gb = SeqIO.read(gb_path, 'gb')
    
    def make_feature_dict():
        '''Make dictionary of features for genbank file mapping "label"
        qualifier to the feature

        Returns:
            dict: Dict of labels to SeqFeatures
        '''
        fd = {}
        for each_feat in gb.features:
            if 'label' in each_feat.qualifiers:
                label = each_feat.qualifiers['label']
                if isinstance(label, list) or isinstance(label, tuple):
                    label = label[0]
                fd[label] = each_feat
        
        return fd
    
    
    
    feat_dict = make_feature_dict()
    
    if start_feature_label in feat_dict:
        # get start location of start feature
        start_feat = feat_dict[start_feature_label]
        start = int(start_feat.location.start)
        
        # slice the sequence so start feature is at start
        sliced_seq = gb.seq[start:] + gb.seq[:start]
        
        assert len(sliced_seq) == len(gb.seq)
        
        # write sliced sequences as fasta to output
        sliced_record = SeqRecord(sliced_seq, id=gb.name)
        SeqIO.write([sliced_record], sliced_output_path, format='fasta')
        
        
    
    else:
        raise TypeError('Start feature label is not present in provided gb file.')


def main():
    
    gb_path = snakemake.params['genbank']
    start_feature_label = snakemake.params['start_feature']
    output = str(snakemake.output)
    
    reslice_genbank(gb_path, start_feature_label, output)


if __name__ == '__main__':
    main()
                
                
    