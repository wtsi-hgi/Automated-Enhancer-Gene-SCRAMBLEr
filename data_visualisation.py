import matplotlib.pyplot as plt

import data_initialisation as di

def data_exploration(expression):
    print("Visualising data...")
    global relative_expression, overlaps, genes
    relative_expression = expression.sort_values(by = ["Std"], ascending = False)
    #relative_expression = relative_expression[relative_expression["mean"] < relative_expression["HAP1"]]
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
    
def plot_convolutions(gene, step_x, gene_basewise, conv_x, gene_convolution, peaks):
    
    figure, axis = plt.subplots(2, 1, figsize = (18.5, 10.5))
    axis[0].plot(step_x, gene_basewise, c = "green")
    axis[0].set(xlabel = "Coordinate on chromosome " + gene["Chromosome"])
    axis[0].plot(conv_x, gene_convolution, c = "orange")
    axis[0].plot((peaks + (gene.Start - ((len(conv_x) - len(step_x)) / 2))), gene_convolution[peaks], "x")
    #axis[0].plot(conv_x, peaks, "x")
    #axis[0].set(xlabel = "Coordinate on chromosome " + gene["Chromosome"])
    
    plt.savefig(di.RESULTS_DIRECTORY + gene["Gene_name"] + "_enhancer_convolution")
    plt.close()