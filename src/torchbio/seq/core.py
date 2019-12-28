import torch

from Bio import SeqIO
from . import const

class GenomeSeqEncoder(object):
    def __init__(self, genome_fp, filetype='fasta', verbose=False):
        self.genome_fp = genome_fp
        self.verbose = verbose
        self.filetype = filetype

        self.genome_dict = SeqIO.to_dict(SeqIO.parse(genome_fp, filetype))

    def get_seq(self, chrom, start, end):
        return self.genome_dict[chrom][start:end].seq.upper()

    def get_encoded_seq(self, chrom, start, end):
        seq = self.get_seq(chrom, start, end)
        
        for base in seq:
            if base not in const.SEQ2INT:
                raise ValueError(f'Invalid base: {base}')

        return torch.eye(len(seq))[[const.SEQ2INT[base] for base in seq]]

