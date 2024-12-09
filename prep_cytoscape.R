# prep_cytoscape.R
# Author : Mathieu G.
# Date : 3/12/2024
# Description : Generates a sif file that can be imported to cytoscape
#               from double mutants amino acid sequences
#               with hamming distance of 1.

# Load the hamming distance calculator function.
source(file="hamming.R")

# CHANGE THIS VARIABLE TO ONE OF THESE:
# {"anidulafungin", "caspofungin", "micafungin"}
COMPOUND <- "micafungin"


df_doubles <- read.csv("double_aa_classification.csv")

# Keep only mutants that are resistant to a specific antifungal.
compound_df_doubles <- df_doubles[df_doubles$compound == COMPOUND,]
comp_res_df_doubles <- compound_df_doubles[compound_df_doubles$sensres == "resistant",]

Lines <- c()

# Compare each pair of sequences.
# If the hamming distance == 1, add that pair to the graph.
for (i in comp_res_df_doubles$aa_seq) {
    for (j in comp_res_df_doubles$aa_seq){
        if (hamming(i,j) == 1){
            line_to_add <- paste(i, "hamming1", j, sep = " ")
            Lines <- c(Lines, line_to_add)
        }
    }
    
}

writeLines(Lines, paste(COMPOUND, "_graph.sif", sep = ""))