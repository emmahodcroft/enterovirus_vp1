import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from augur.parse import forbidden_characters

if __name__ == '__main__':

    import argparse
    parser = parser = argparse.ArgumentParser(description='parse metadata',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--blast', help="blast result file")
    parser.add_argument('--meta', help="input meta from genbank")
    parser.add_argument('--seqs', help="input sequences from genbank")
    parser.add_argument('--out_seqs', help="output seqs")
    parser.add_argument('--out_meta', help="output meta")
    args = parser.parse_args()

    #read in blast results and give them a header
    blast = pd.read_csv(args.blast, index_col=False,
        names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs"])

    #pick the ones we want - matches of at least 700 bp. 
    # In initial tests, 600bp added only 47 sequences more, 800bp lost 289 sequences.
    keep_blast = list(blast.loc[blast.send - blast.sstart >= 700, "qseqid"])

    print("{} of the {} sequences have at least 700bp in VP1 and will be included".format(len(keep_blast), 
        len(set(blast.qseqid))))

    # These meta strain names are still uncorrected for forbidden chars, but the sequence names have been corrected...
    # so we need to make them comparable by removing the forbidden chars before we can find them in the meta...
    meta = pd.read_csv(args.meta, sep='\t', index_col=False)
    keep_meta_forbid_names = []
    for f in meta.strain:
        name = f
        for c,r in forbidden_characters:
            name=name.replace(c,r)
        if name in keep_blast:
            keep_meta_forbid_names.append(f)
    keep_meta = meta[meta['strain'].isin(keep_meta_forbid_names)]
    keep_meta.to_csv(args.out_meta, sep='\t', index=False)

    # Now find them in the fasta, but only take the matching area of the sequence
    fasta_sequences = SeqIO.parse(open(args.seqs), 'fasta')
    keepSeqs = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if name in keep_blast: 
            strt = int(blast.loc[blast.qseqid == name].qstart-1)
            en = int(blast.loc[blast.qseqid == name].qend)
            vp1_seq = sequence[strt:en]
            keepSeqs.append(SeqRecord(Seq(vp1_seq), id=name, name=name, description=""))
    SeqIO.write(keepSeqs, args.out_seqs, "fasta")




