#!/usr/bin/env python3

# Alessandro Lussana <alussana@ebi.ac.uk>

from phosx.phosx import compute_kinase_activity
from os import path

if __name__ == '__main__':
    compute_kinase_activity(
        seqrnk_file=str(path.join(path.dirname(__file__), 'seqrnk/koksal2018_log2.fold.change.8min.seqrnk')),
        pssms_h5_file=str(path.join(path.dirname(__file__), '../phosx/data/pssm.h5')),
        n_perm=10,
    )