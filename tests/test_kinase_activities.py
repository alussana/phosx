#!/usr/bin/env python3

from phosx.pssm_enrichments import kinase_activities
from os import path

def test_kinase_activities_1core():
    
    kinase_activities(
    
        seqrnk_file=str(path.join(path.dirname(__file__), 'seqrnk/koksal2018_log2.fold.change.8min.seqrnk')),
        pssms_h5_file=str(path.join(path.dirname(__file__), '../phosx/data/pssm.h5')),
        n_perm=10,
        plot_figures=True
    
    )
    
def test_kinase_activities_2cores():
    
    kinase_activities(
    
        seqrnk_file=str(path.join(path.dirname(__file__), 'seqrnk/koksal2018_log2.fold.change.8min.seqrnk')),
        pssms_h5_file=str(path.join(path.dirname(__file__), '../phosx/data/pssm.h5')),
        n_perm=10,
        n_proc=2,
        plot_figures=True
    
    )
    
def test_all():
    test_kinase_activities_1core()
    test_kinase_activities_2cores()

if __name__ == '__main__':
    
    test_all()