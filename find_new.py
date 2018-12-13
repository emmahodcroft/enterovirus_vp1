import pandas as pd
from dateutil.parser import parse
import re, os
import shutil


if __name__ == '__main__':
    import argparse

    parser = parser = argparse.ArgumentParser(description='find new sequences',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input-new', help="new meta file")
    parser.add_argument('--input-old', help="old meta file")
    parser.add_argument('--input-swedish', help="swedish meta info")
    parser.add_argument('--output', help="new meta lines")
    args = parser.parse_args()


    new_meta = pd.read_csv(args.input_new, sep="\t", index_col=False)
    old_meta = pd.read_csv(args.input_old, sep="\t", index_col=False)
    swed_meta = pd.read_csv(args.input_swedish, sep=",", index_col=False)

    newM = set(new_meta.accession)-set(old_meta.accession)-set(swed_meta.accession)

    new_actual_meta = new_meta.loc[new_meta['accession'].isin(newM)]

    print("{} new accession numbers were found and will be downloaded.".format(len(newM)))

    if not newM:
        print("No new accession numbers were added. Run will not proceed!")
    else:
        new_actual_meta.to_csv(args.output, sep='\t', index=False)
        


