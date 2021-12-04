from Bio import SeqIO
from argparse import ArgumentParser
from pathlib import Path
import re


def get_args():
    
    parser = ArgumentParser('Split fasta file into individual records.')
    parser.add_argument('record', help='Path to fasta to split.')
    parser.add_argument('output', help='Output directory to store split records.')
    parser.add_argument('--filter', default=None)
    parser.add_argument('--rename', help='', default=None)
    
    return parser.parse_args()


def read_record(filepath):
    
    return SeqIO.parse(str(filepath), 'fasta')


def write_records(records, outdir, rename_re=None, filter_re=None):
    
    for each_record in records:
        if rename_re:
            try:
                name = re.search(rename_re, each_record.description)[0]
            except TypeError as e:
                print(f'Could not parse {each_record.description}')
                continue
        else:
            name = each_record.description
        
        # clean for filename
        filename = Path(name.replace(' ', '_')).with_suffix('.fa')
        filepath = Path(outdir).joinpath(filename)
        
        if filter_re:
            filter_match = re.search(filter_re, each_record.description)

            if filter_match:
                SeqIO.write(each_record, str(filepath), 'fasta')
        else:
            SeqIO.write(each_record, str(filepath), 'fasta')


def main():
    
    args = get_args()
    records = read_record(args.record)
    write_records(
        records, outdir=args.output, rename_re=args.rename, 
        filter_re=args.filter)


if __name__ == '__main__':
    main()
    


        

        
        
        
        
        

