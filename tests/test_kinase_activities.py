#!/usr/bin/env python3

from phosx.pssm_enrichments import kinase_activities
from os import path

if __name__ == '__main__':
    
    kinase_activities(
    
        seqrnk_file=str(path.join(path.dirname(__file__), 'seqrnk/koksal2018_log2.fold.change.8min.seqrnk')),
        pssms_h5_file=str(path.join(path.dirname(__file__), '../phosx/data/pssm.h5')),
        n_perm=10,
    
    )