import matplotlib.pyplot as plt
import numpy as np

import data_initialisation as di

def data_exploration(expression):
    print("Visualising data...")
    global relative_expression, overlaps, genes
    relative_expression = expression.sort_values(by = ["Std"], ascending = False)
    relative_expression["Difference"] = relative_expression[di.CELL_LINE_OF_INTEREST] - relative_expression["Mean"]

    figure, axis = plt.subplots(1, 3, figsize = (18.5, 10.5))
    axis[0].scatter(relative_expression["Std"], relative_expression[di.CELL_LINE_OF_INTEREST], s = 0.3, c = "red")
    axis[0].scatter(relative_expression["Std"], relative_expression["Mean"], s = 0.3, c = "blue")
    axis[0].set(xlabel = "Standard deviation of expression within cell lines", ylabel = "Normalised expression")

    axis[1].scatter(relative_expression["Difference"], relative_expression[di.CELL_LINE_OF_INTEREST], s = 0.3, c = "red")
    axis[1].scatter(relative_expression["Difference"], relative_expression["Mean"], s = 0.3, c = "blue")
    axis[1].set(xlabel = "Difference between expression within " + di.CELL_LINE_OF_INTEREST + " and mean expression", ylabel = "Normalised expression")

    axis[2].scatter(relative_expression["Std"], relative_expression["Difference"], s = 0.3, c = "green")
    axis[2].set(xlabel = "Difference between expression within " + di.CELL_LINE_OF_INTEREST + " and mean expression", ylabel = "Normalised expression")

    plt.savefig(di.RESULTS_DIRECTORY + "Significant_" + di.CELL_LINE_OF_INTEREST + "_expression.png")
    plt.close(figure)
    
    figure, axis = plt.subplots(2, 2, figsize = (18.5, 10.5))
    axis[0][0].scatter(genes["Enhancer_count"], genes["Std"], s = 0.3, c = "purple")
    axis[0][0].set(xlabel = "Number of enhancers within search area", ylabel = "Standard deviation of expression within cell lines")
    
    axis[1][0].scatter(genes["Enhancer_count"], genes[di.CELL_LINE_OF_INTEREST], s = 0.3, c = "red")
    axis[1][0].scatter(genes["Enhancer_count"] + 0.5, genes["Mean"], s = 0.3, c = "blue")
    axis[1][0].set(xlabel = "Number of enhancers within search area", ylabel = "Normalised expression")
    
    axis[0][1].scatter(genes["Enhancer_proportion"], genes["Std"], s = 0.3, c = "purple")
    axis[0][1].set(xlabel = "Proportion of search area which is enhancer", ylabel = "Standard deviation of expression within cell lines")
    
    axis[1][1].scatter(genes["Enhancer_proportion"], genes[di.CELL_LINE_OF_INTEREST], s = 0.3, c = "red")
    axis[1][1].scatter(genes["Enhancer_proportion"], genes["Mean"], s = 0.3, c = "blue")
    axis[1][1].set(xlabel = "Proportion of search area which is enhancer", ylabel = "Normalised expression")
    
    plt.savefig(di.RESULTS_DIRECTORY + "Enhancers_near_genes.png")
    plt.close(figure)
    
def plot_convolutions(gene, step_x, gene_basewise, conv_x, gene_convolution):
    
    figure, axis = plt.subplots(2, 1, figsize = (18.5, 10.5))
    axis[0].plot(step_x, gene_basewise, c = "green")
    axis[0].set(xlabel = "Coordinate on chromosome " + gene["Chromosome"])
    axis[0].plot(conv_x, gene_convolution, c = "orange")
    
    plt.savefig(di.RESULTS_DIRECTORY + gene["Gene_name"] + "_enhancer_convolution")
    plt.close()
    
    
def compare_metrics(gene_data, title, filename):
    
    print("Graphing metric comparisons...")
    
    figure, axis = plt.subplots(6, 6, figsize = (15, 15))
    figure.suptitle(title, fontsize = 20)
    
    axis[0][0].scatter(gene_data["Std"], gene_data["Std"], s = 0.2, c = "green")
    axis[1][0].scatter(gene_data["Std"], gene_data["Anomalous_score"], s = 0.2, c = "green")
    axis[2][0].scatter(gene_data["Std"], gene_data["Specific_gene_expression"], s = 0.2, c = "green")
    axis[3][0].scatter(gene_data["Std"], gene_data["Enhancer_count"], s = 0.2, c = "green")
    axis[4][0].scatter(gene_data["Std"], gene_data["Enhancer_proportion"], s = 0.2, c = "green")
    axis[5][0].scatter(gene_data["Std"], gene_data["Gene_size"], s = 0.2, c = "green")
    
    axis[0][1].scatter(gene_data["Anomalous_score"], gene_data["Std"], s = 0.2, c = "green")
    axis[1][1].scatter(gene_data["Anomalous_score"], gene_data["Anomalous_score"], s = 0.2, c = "green")
    axis[2][1].scatter(gene_data["Anomalous_score"], gene_data["Specific_gene_expression"], s = 0.2, c = "green")
    axis[3][1].scatter(gene_data["Anomalous_score"], gene_data["Enhancer_count"], s = 0.2, c = "green")
    axis[4][1].scatter(gene_data["Anomalous_score"], gene_data["Enhancer_proportion"], s = 0.2, c = "green")
    axis[5][1].scatter(gene_data["Anomalous_score"], gene_data["Gene_size"], s = 0.2, c = "green")
    
    axis[0][2].scatter(gene_data["Specific_gene_expression"], gene_data["Std"], s = 0.2, c = "green")
    axis[1][2].scatter(gene_data["Specific_gene_expression"], gene_data["Anomalous_score"], s = 0.2, c = "green")
    axis[2][2].scatter(gene_data["Specific_gene_expression"], gene_data["Specific_gene_expression"], s = 0.2, c = "green")
    axis[3][2].scatter(gene_data["Specific_gene_expression"], gene_data["Enhancer_count"], s = 0.2, c = "green")
    axis[4][2].scatter(gene_data["Specific_gene_expression"], gene_data["Enhancer_proportion"], s = 0.2, c = "green")
    axis[5][2].scatter(gene_data["Specific_gene_expression"], gene_data["Gene_size"], s = 0.2, c = "green")
    
    axis[0][3].scatter(gene_data["Enhancer_count"], gene_data["Std"], s = 0.2, c = "green")
    axis[1][3].scatter(gene_data["Enhancer_count"], gene_data["Anomalous_score"], s = 0.2, c = "green")
    axis[2][3].scatter(gene_data["Enhancer_count"], gene_data["Specific_gene_expression"], s = 0.2, c = "green")
    axis[3][3].scatter(gene_data["Enhancer_count"], gene_data["Enhancer_count"], s = 0.2, c = "green")
    axis[4][3].scatter(gene_data["Enhancer_count"], gene_data["Enhancer_proportion"], s = 0.2, c = "green")
    axis[5][3].scatter(gene_data["Enhancer_count"], gene_data["Gene_size"], s = 0.2, c = "green")
    
    axis[0][4].scatter(gene_data["Enhancer_proportion"], gene_data["Std"], s = 0.2, c = "green")
    axis[1][4].scatter(gene_data["Enhancer_proportion"], gene_data["Anomalous_score"], s = 0.2, c = "green")
    axis[2][4].scatter(gene_data["Enhancer_proportion"], gene_data["Specific_gene_expression"], s = 0.2, c = "green")
    axis[3][4].scatter(gene_data["Enhancer_proportion"], gene_data["Enhancer_count"], s = 0.2, c = "green")
    axis[4][4].scatter(gene_data["Enhancer_proportion"], gene_data["Enhancer_proportion"], s = 0.2, c = "green")
    axis[5][4].scatter(gene_data["Enhancer_proportion"], gene_data["Gene_size"], s = 0.2, c = "green")
    
    axis[0][5].scatter(gene_data["Gene_size"], gene_data["Std"], s = 0.2, c = "green")
    axis[1][5].scatter(gene_data["Gene_size"], gene_data["Anomalous_score"], s = 0.2, c = "green")
    axis[2][5].scatter(gene_data["Gene_size"], gene_data["Specific_gene_expression"], s = 0.2, c = "green")
    axis[3][5].scatter(gene_data["Gene_size"], gene_data["Enhancer_count"], s = 0.2, c = "green")
    axis[3][5].scatter(gene_data["Gene_size"], gene_data["Enhancer_proportion"], s = 0.2, c = "green")
    axis[5][5].scatter(gene_data["Gene_size"], gene_data["Gene_size"], s = 0.2, c = "green")
    
    axis[0][0].set(ylabel = "Standard Deviation")
    axis[1][0].set(ylabel = "Anomalous Score")
    axis[2][0].set(ylabel = "Expression of Gene of Interest")
    axis[3][0].set(ylabel = "Number of Enhancers")
    axis[4][0].set(ylabel = "Denisty of Enhancers")
    axis[5][0].set(xlabel = "Standard Deviation", ylabel = "Size of Gene")
    axis[5][1].set(xlabel = "Anomalous Score")
    axis[5][2].set(xlabel = "Expression of Gene of Interest")
    axis[5][3].set(xlabel = "Number of Enhancers")
    axis[5][4].set(xlabel = "Denisty of Enhancers")
    axis[5][5].set(xlabel = "Size of Gene")
            
    plt.savefig(di.RESULTS_DIRECTORY + filename)
    plt.close()
    
def gene_report(gene_data):
    
    gene_data = gene_data.sort_values("Interest_score", ascending = False)
    
    for index, gene in gene_data.iterrows():
    
        print("Creating gene report for " + gene["Gene_name"] + "...")
        
        print(gene)
        
        step_x = np.frombuffer(gene["Enhancer_searched_coordinates"], dtype = float)
        step_y = np.frombuffer(gene["Enhancer_step_function"], dtype = float)
        convolved_x = np.frombuffer(gene["Enhancer_convolved_coordinates"], dtype = float)
        convolved_y = np.frombuffer(gene["Enhancer_convolution"], dtype = float)
        
        print(step_x)
        print(step_y)
        print(convolved_x)
        print(convolved_y)
        print(gene["Plateau_starts"])
        print(gene["Plateau_ends"])
        
        
        figure, axis = plt.subplots(4, 1, figsize = (12, 12))
        
        axis[0].axvline(x = gene["Gene_start"], c = "green")
        axis[0].axvline(x = gene["Gene_end"], c = "green")
        axis[0].axvline(x = gene["Search_window_end"], c = "red")
        axis[0].axvline(x = gene["Search_window_end"], c = "red")
        
        axis[1].plot(step_x, step_y, c = "green")
        axis[1].plot(convolved_x, convolved_y, c = "orange")
        axis[1].axhline(y = di.PLATEAU_THRESHOLD, c = "red")
        
        for plateau_start in gene["Plateau_starts"]:
            
            axis[2].axvline(x = plateau_start, c = "red")
        
        for plateau_end in gene["Plateau_ends"]:
            
            axis[2].axvline(x = plateau_end, c = "orange")
        
        axis[3].axvline(x = gene["Gene_start"], c = "green")
        axis[3].axvline(x = gene["Gene_end"], c = "green")
        axis[3].axvline(x = gene["Search_window_end"], c = "red")
        axis[3].axvline(x = gene["Search_window_end"], c = "red")
        
        axis[3].plot(step_x, step_y, c = "green")
        axis[3].plot(convolved_x, convolved_y, c = "orange")
        axis[3].axhline(y = di.PLATEAU_THRESHOLD, c = "red")
        
        for plateau_start in gene["Plateau_starts"]:
            
            axis[3].axvline(x = plateau_start, c = "red")
        
        for plateau_end in gene["Plateau_ends"]:
            
            axis[3].axvline(x = plateau_end, c = "orange")
        
        plt.savefig(di.RESULTS_DIRECTORY + gene["Gene_name"] + "_enhancer_convolution")
        plt.close()
        
        if index > di.ENHANCER_CONVOLUTION:
            break