import torch

from Bio import SeqIO
import const


class GenomeSeqEncoder(object):
    def __init__(self, genome_fp, filetype='fasta', verbose=False):
        self.genome_fp = genome_fp
        self.verbose = verbose
        self.filetype = filetype

        self.genome_dict = SeqIO.to_dict(SeqIO.parse(genome_fp, filetype))

    def get_seq(self, chrom, start, end):
        return self.genome_dict[chrom][start:end].seq.upper()

    def encode_region(self, chrom, start, end, max_len=None):
        seq = self.get_seq(chrom, start, end)

        if max_len and len(seq) > max_len:
            raise ValueError(f'Length of the sequence should not exceed max_len. {len(seq)} > {max_len}')

        for base in seq:
            if base not in const.BASE2INT:
                raise ValueError(f'Invalid base: {base}')

        if not max_len:
            return self.encode_seq(seq)
        else:
            return torch.cat([self.encode_seq(seq), torch.zeros(4, max_len - len(seq))], axis=1)

    def encode_seq(self, seq):
        """Encode DNA sequence into 4 x (sequence_length) tensor,
        where the order of nucleotides are A, C, G, and T.
        """
        return torch.eye(len(const.BASE2INT))[[const.BASE2INT[base] for base in seq]][:, :4].T

    def decode(self, encoded):
        indices = torch.argmax(encoded, axis=0).numpy()
        return ''.join([const.INT2BASE[i] for i in indices])


if __name__ == '__main__':
    encoder = GenomeSeqEncoder('/data/project/dohoon/deepmeth/reference/hg38.fa')
    encoded = encoder.encode_region('chr1', 10000000, 10000010)
    decoded = encoder.decode(encoded)

    print(encoded)
    print(encoded.shape)
    print(decoded)
    print(decoded)
