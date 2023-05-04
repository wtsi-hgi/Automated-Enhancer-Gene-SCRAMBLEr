import sys
import hashlib
import pandas as pd

import data_initialisation as di
import sequence_seeking as ss

interesting_features = ["Std", "Anomalous_score", "Enhancer_count",
                        "Enhancer_proportion", "Specific_gene_expression",
                        "Gene_size", "Symmetry_ratio"]

def export_gene_scores_report(gene_data):
    
    #Idealy this will not read from file but from passed argument
    
    #Md5 checksum of config file is generated. Gene prioritisation report file
    #is created and checksum is included in name to differentiate different
    #configs. Report saved in given location.
    
    print("Exporting gene prioritisation report...")
    
    checksum = generate_config_checksum()
    
    with open(sys.argv[1], "r") as config:
        
        report_name = \
            "gene_prioritisation_report_" + checksum.hexdigest() + ".txt"
        report = \
            open((di.GENE_PRIORITISATION_REPORT_DIRECTORY + report_name), "w")
        report.write(config.read() + "\n")
        report.close()
        report = \
            open((di.GENE_PRIORITISATION_REPORT_DIRECTORY + report_name), "a")
        gene_data.loc[:, (["Gene_name"] +
                          ["Interest_score"] + 
                          interesting_features +
                          ["Scaled_std",
                           "Scaled_anomalous_score",
                           "Scaled_enhancer_count",
                           "Scaled_enhancer_proportion",
                           "Scaled_specific_gene_expression",
                           "Scaled_gene_size",
                           "Scaled_symmetry_ratio",
                           "Z-Std",
                           "Z-Anomalous_score",
                           "Z-Enhancer_count",
                           "Z-Enhancer_proportion",
                           "Z-Specific_gene_expression",
                           "Z-Gene_size",
                           "Z-Symmetry_ratio"])].to_csv(
            (di.GENE_PRIORITISATION_REPORT_DIRECTORY + report_name),
            sep = "\t", index = True, mode = "a")            
        report.close()
        
def generate_config_checksum():

    checksum = hashlib.md5()
    
    with open(sys.argv[1], "rb") as config:
        
        for chunk in iter(lambda: config.read(4096), b""):
            
            checksum.update(chunk)
        
    return checksum

def export_plateaus(gene_data):
    
    # export_plateaus saves plateaus associated with each gene as a bed file
    
    print("Exporting plateaus to bed file...")
    
    with open((di.RESULTS_DIRECTORY + "plateaus.bed"), "w") as f:
            f.write("Chromosome Start	End	Gene_name")
            f.write("\n")

    for index, gene in gene_data.head(di.CONVOLUTION_LIMIT).iterrows():
        
        plateaus = pd.DataFrame(
            {"Start" : gene_data.loc[index, "Plateau_starts"], 
             "End" : gene_data.loc[index, "Plateau_ends"]
            }
        )
        plateaus["Gene_name"] = gene["Gene_name"]
        plateaus["Chromosome"] = "chr" + gene["Chromosome"]
        plateaus["Strand"] = gene["Strand"]
        
        plateaus = ss.find_fasta(plateaus)
        ss.generate_pridict_input(plateaus)
        
def export_convolutions(gene_data):
    
    # Coordinates of convolutions are exported to wig file, for each gene
    
    print("Exporting enhancer density convolutions to wig file...")
    
    for index, gene in gene_data.head(di.CONVOLUTION_LIMIT).iterrows():
        with open((di.RESULTS_DIRECTORY +
            gene["Gene_name"] + "_convolutions.wiggle"), "w") as f:
            
            f.write("fixedStep chrom=chr" +
                gene["Chromosome"] + " start=" +
                    str(gene["Enhancer_convolution_x"][0]) + " step=1")
            f.write("\n")
    
        convolution_signal = pd.DataFrame({"Convolution_signal" : \
            gene_data.loc[index, "Enhancer_convolution_y"]})
        convolution_signal.to_csv(
            (di.RESULTS_DIRECTORY + gene["Gene_name"] + "_convolutions.wig"), 
            sep = "\t", 
            index = False, 
            mode = "a", 
            header = False)