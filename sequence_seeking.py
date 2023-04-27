import pyranges as pr
import pandas as pd
import re

import data_initialisation as di

global possible_plateau_insertions

possible_plateau_insertions = pd.DataFrame(columns = ["Sequence_name", "Insertion_sequence", "Insertion_location", "Plateau_sequence"])

def find_fasta(plateaus):
    
    #find_fasta takes coordinates of plateaus and returns FASTA sequence from reference genome
    
    print("Finding FASTA sequences from intervals...")
    
    plateaus_pr = pr.PyRanges(plateaus)
    
    seq = pr.get_sequence(plateaus_pr, di.REFERENCE_GENOME)
    plateaus_pr.seq = seq
    plateaus = plateaus_pr.df
    plateaus = plateaus.rename(columns = {"seq" : "Sequence"})
    
    return plateaus

def generate_pridict_input(plateaus):
    
    global possible_plateau_insertions
    
    plateaus.apply(generate_insertion_prefixes_and_suffixes, axis = 1)
    
    possible_plateau_insertions["Sequence"] = possible_plateau_insertions.apply(insert_insertion_sequence, axis = 1)
    
    #possible_plateau_insertions["Sequence"] = possible_plateau_insertions.apply(lambda insertion : insertion["Plateau_sequence"][:insertion["Insertion_location"]] + "(+" + insertion["Insertion_sequence"] + ")" + insertion["Plateau_sequence"][insertion["Insertion_location"]:])
    
    print(possible_plateau_insertions)
    
    possible_plateau_insertions.to_csv((di.RESULTS_DIRECTORY + "sequences_for_pridict.csv"), index = False, columns = ["Sequence_name", "Sequence"], mode = "w", header = False)

def generate_insertion_prefixes_and_suffixes(plateau):
    
    global possible_plateau_insertions
    
    plateau_specific_suggested_insertion_sites = pd.DataFrame(columns = ["Sequence_name", "Insertion_sequence", "Insertion_location", "Plateau_sequence"])
    
    for number_of_bases_absent in range(0, len(di.INSERTED_SEQUENCE)):
    
        for insertion_sequence_position_being_checked in range(0, (len(di.INSERTED_SEQUENCE) - number_of_bases_absent)):
            
            absent_sequence = di.INSERTED_SEQUENCE[insertion_sequence_position_being_checked:(insertion_sequence_position_being_checked + number_of_bases_absent)]
            present_sequence = di.INSERTED_SEQUENCE[:insertion_sequence_position_being_checked] + di.INSERTED_SEQUENCE[(insertion_sequence_position_being_checked + number_of_bases_absent):]
    
            insertion_positions = find_prefix_suffix_in_plateau(plateau, present_sequence)
            
            for position in insertion_positions:
            
                new_row = pd.Series({"Sequence_name" : (plateau["Gene_name"] + " " + plateau["Chromosome"] + " " + plateau["Strand"] + " "  + str(plateau["Start"]) + "-" + str(plateau["End"])), "Insertion_sequence" : absent_sequence, "Insertion_location" : (position + insertion_sequence_position_being_checked), "Plateau_sequence" : plateau["Sequence"]})
                new_df = pd.DataFrame([new_row])
                plateau_specific_suggested_insertion_sites = pd.concat([plateau_specific_suggested_insertion_sites, new_df], axis = 0, ignore_index = True)
                
                if len(plateau_specific_suggested_insertion_sites.index) > di.PARTIAL_INSERTIONS_PER_REGION:
                    
                    possible_plateau_insertions = pd.concat([possible_plateau_insertions, plateau_specific_suggested_insertion_sites], axis = 0, ignore_index = True)
                    
                    break
    
def find_prefix_suffix_in_plateau(plateau, present_sequence):
    
    insertions = re.finditer(pattern = present_sequence, string = plateau["Sequence"])
    insertion_positions = [index.start() for index in insertions]
    
    return insertion_positions

def insert_insertion_sequence(row):
    
    return (row["Plateau_sequence"][:row["Insertion_location"]] + 
                "(+" +
                row["Insertion_sequence"] +
                ")" + 
                row["Plateau_sequence"][row["Insertion_location"]:])
        