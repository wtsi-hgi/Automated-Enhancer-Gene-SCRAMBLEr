import pyranges as pr

import data_initialisation as di

def find_fasta(plateau_regions):
    
    #find_fasta takes coordinates of plateaus and returns FASTA sequence from reference genome
    
    print("Finding FASTA sequences from intervals...")
    
    plateau_regions_pr = pr.PyRanges(plateau_regions)
    
    seq = pr.get_sequence(plateau_regions_pr, di.REFERENCE_GENOME)
    plateau_regions_pr.seq = seq
    plateau_regions = plateau_regions_pr.df
    
    return plateau_regions