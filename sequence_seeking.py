import pyranges as pr

import data_initialisation as di

def find_fasta(gene_data):
    
    print("Finding FASTA sequences from intervals...")
    
    for index, gene in gene_data.iterrows():
        
        sequences_df = gene.loc(:, ["Chromesome", "Plateau_starts", "Plateau_ends", "Strand"])
    
        pr.get_fasta.get_sequence()