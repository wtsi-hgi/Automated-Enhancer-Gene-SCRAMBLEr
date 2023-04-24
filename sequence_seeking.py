import pyranges as pr
import pandas as pd

import data_initialisation as di

def find_fasta(plateaus):
    
    #find_fasta takes coordinates of plateaus and returns FASTA sequence from reference genome
    
    print("Finding FASTA sequences from intervals...")
    
    plateaus_pr = pr.PyRanges(plateaus)
    
    seq = pr.get_sequence(plateaus_pr, di.REFERENCE_GENOME)
    plateaus_pr.seq = seq
    plateaus = plateaus_pr.df
    plateaus = plateaus.rename(columns = {"seq" : "Sequence"})
    
    return plateaus
    
def find_insertion_prefixes_and_suffixes(plateaus):

    #Within sequences of regions, finds prefixes and suffixes of sequence that is to be inserted

    print("Finding possible insertion sites within each plateau...")

    sequences_for_pridict = pd.DataFrame(columns = ["Sequence_name", "Sequence"])

    #plateaus = 

    for _, plateau in plateaus.iterrows():

        for amount_unfound in range(1, len(di.INSERTED_SEQUENCE)):
            
            for checking_insertion_position in range(1, len(di.INSERTED_SEQUENCE)):
                
                missing_insertion = di.INSERTED_SEQUENCE[checking_insertion_position:(checking_insertion_position + amount_unfound)]
                found_insertion = di.INSERTED_SEQUENCE[:checking_insertion_position] + di.INSERTED_SEQUENCE[(checking_insertion_position + amount_unfound):]
                    
                for checking_plateau_position in range(len(plateau)):
                
                    if plateau[checking_plateau_position:].str.startswith(found_insertion):
                        
                        sequence_with_insertion = plateau[:(checking_plateau_position + checking_insertion_position)] + "(+" + missing_insertion + ")" + plateau[(checking_plateau_position + checking_insertion_position):]
                        new_row = {"Sequence_name" : (plateau["Gene_name"] + " " + plateau["Chromosome"] + " " + plateau["Strand"]), "Sequence" : sequence_with_insertion}
                        
                        sequences_for_pridict = sequences_for_pridict.append(new_row, ignore_index = True)
                    
    return sequences_for_pridict
        