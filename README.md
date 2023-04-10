# AEG SCRAMBLE
A tool to help automate the Genome Scramble project

Jonas Koeppel has started a project with Leopold Parts which will attempt to introduce controlled breaks between clusters of enhancers near genes and recombine the piece in a random manner to investigate phenotypic changes.
After identifying one gene (OTX2) which was suitable, Jonas has approached HGI to look at automating the process, which Leo broke into steps, each dependent on the success of the last:
  - Prioritise the genes of interest
  - Identify regulatory elements around the gene
  - Finding pegRNA insertion sites
  - Designing pegRNA sequences
  
## Prioritise the genes of interest
For this step, a number of metrics need to found for each gene, so that they can be compared and filtered to find the genes of interest.
* Whether the gene is a housekeeping gene
  * As a first approach for finding this, I will find the standard deviation of expression between cell lines. Jonas has provided a table of expression of many genes in many cell lines.
* Expression within HAP1
  * For now Jonas is interested in any level of expression within HAP1 cells, but it is something to take note of for the future.
* Number of nearby enhancers
  * Given a search area upstream and downstream of each gene, how many enhancers are present; and what proportion of this area do they take?
