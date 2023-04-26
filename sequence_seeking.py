import pyranges as pr
import pandas as pd
import re

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

def iterate_prefix_suffix_search(plateaus):
    
    plateaus.apply(generate_insertion_prefixes_and_suffixes, axis = 1)

def generate_insertion_prefixes_and_suffixes(plateau):
    
    di.INSERTED_SEQUENCE
    
    plateau_specific_suggested_insertion_sites = pd.DataFrame(columns = ["Plateau", "Insertion_sequence", "Insertion_location"])
    
    for number_of_bases_absent in range(0, len(di.INSERTED_SEQUENCE)):
    
        for insertion_sequence_position_being_checked in range(0, (len(di.INSERTED_SEQUENCE) - number_of_bases_absent)):
            
            absent_sequence = di.INSERTED_SEQUENCE[insertion_sequence_position_being_checked:(insertion_sequence_position_being_checked + number_of_bases_absent)]
            present_sequence = di.INSERTED_SEQUENCE[:insertion_sequence_position_being_checked] + di.INSERTED_SEQUENCE[(insertion_sequence_position_being_checked + number_of_bases_absent):]
    
            insertion_positions = find_prefix_suffix_in_plateau(plateau, present_sequence)
            
            for position in insertion_positions:
            
                new_row = pd.Series({"Plateau" : (plateau["Gene_name"] + " " + str(plateau["Start"]) + "-" + str(plateau["End"])), "Insertion_sequence" : absent_sequence, "Insertion_location" : position})
                #plateau_specific_suggested_insertion_sites = plateau_specific_suggested_insertion_sites.append(new_row, ignore_index = True)
                new_df = pd.DataFrame([new_row])
                plateau_specific_suggested_insertion_sites = pd.concat([plateau_specific_suggested_insertion_sites, new_df], axis=0, ignore_index=True)
                
                if len(plateau_specific_suggested_insertion_sites.index) > 100:
                    
                    break
                
    print(plateau_specific_suggested_insertion_sites)
    #plateau_specific_suggested_insertion_sites.to_csv((di.RESULTS_DIRECTORY + (plateau["Gene_name"] + " " + str(plateau["Start"]) + "-" + str(plateau["End"])) + "sequences_for_pridict.csv"), index = False, columns = ["Sequence_name", "Sequence"], mode = "w", header = False)
    
def find_prefix_suffix_in_plateau(plateau, present_sequence):
    
    #insertion_positions = [index for index in range(len(plateau["Sequence"])) if plateau["Sequence"].startswith(present_sequence, index)]
    #print(insertion_positions)
    
    insertions = re.finditer(pattern = present_sequence, string = plateau["Sequence"])
    insertion_positions = [index.start() for index in insertions]
    return insertion_positions
    
def find_insertion_prefixes_and_suffixes2(plateaus):

    #Within sequences of regions, finds prefixes and suffixes of sequence that is to be inserted

    print("Finding possible insertion sites within each plateau...")

    sequences_for_pridict = pd.DataFrame(columns = ["Sequence_name", "Sequence"])

    for _, plateau in plateaus.iterrows():

        print(plateau)

        for amount_unfound in range(1, len(di.INSERTED_SEQUENCE)):
            
            print("Amount unfound: " + str(amount_unfound))
            
            for checking_insertion_position in range(1, len(di.INSERTED_SEQUENCE)):
                
                print("Position within insertion being checked: " + str(checking_insertion_position))
                
                missing_insertion = di.INSERTED_SEQUENCE[checking_insertion_position:(checking_insertion_position + amount_unfound)]
                found_insertion = di.INSERTED_SEQUENCE[:checking_insertion_position] + di.INSERTED_SEQUENCE[(checking_insertion_position + amount_unfound):]
                    
                print("Missing insertion: " + missing_insertion)
                print("Found insertion: " + found_insertion)    
                    
                for checking_plateau_position in range(len(plateau["Sequence"])):
                
                    print("Position within plateau: " + str(checking_plateau_position))
                
                    if "".join(str(plateau[checking_plateau_position:])).startswith(found_insertion):
                        
                        sequence_with_insertion = str(plateau[:(checking_plateau_position + checking_insertion_position)]) + "(+" + missing_insertion + ")" + str(plateau[(checking_plateau_position + checking_insertion_position):])
                        print("Sequence with insertion: " + str(sequence_with_insertion))
                        new_row = {"Sequence_name" : (plateau["Gene_name"] + " " + plateau["Chromosome"] + " " + plateau["Strand"]), "Sequence" : sequence_with_insertion}
                        
                        sequences_for_pridict = sequences_for_pridict.append(new_row, ignore_index = True)
                    
    return sequences_for_pridict
        