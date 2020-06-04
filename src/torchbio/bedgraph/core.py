import pysam
import numpy as np

class Bedgraph(object):

    def __init__(self, f):
        self.f = f

    def fetch(self, *args, **kwargs):
        return pysam.TabixFile(self.f).fetch(*args, **kwargs, parser=pysam.asTuple())

    def vectorize(self, chrom, start, end):
        bdg_rows = list(self.fetch(chrom, start, end))

        vec = []
        for i, row in enumerate(bdg_rows):
            # First row.
            if i == 0:
                # If the whole interval is included in the first row,
                # just return list of scores with length (end - start).
                if int(row[1]) <= start and end < int(row[2]):
                    vec += [float(row[3])] * (end - start)
                else:
                    vec += [float(row[3])] * (int(row[2]) - start)
            # Last row.
            elif i == len(bdg_rows) - 1:
                vec += [float(row[3])] * (end - int(row[1]))
            else:
                vec += [float(row[3])] * (int(row[2]) - int(row[1]))

        return np.array(vec)
