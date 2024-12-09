# epistasis_analysis.R
# Author : Mathieu Giguere
# Date : 4/12/2024
# Description : analyses double mutants that lead to positive epistasis.

library(plotly)

## Choose antifungal compound 
# {"anidulafungin", "caspofungin", "micafungin"}
COMPOUND <- "anidulafungin"

epis_data <- read.csv(paste(COMPOUND, "_epistasis.csv", sep = ""))

epis_data <- epis_data[epis_data$pos_epistasis == TRUE,]

### For each position pair,
### how many times does that position pair lead to positive epistasis ?

epistasis <- c()

counter <- 0

for (i in 1:9){
    for (j in 1:9){
        counter <- nrow(epis_data[epis_data$mut_pos_1 == i & epis_data$mut_pos_2 == j,])
        epistasis <- c(epistasis, counter)
    }
}

epistasis_matrix <- matrix(data=epistasis, nrow = 9, ncol = 9)

text_matrix <- ifelse(epistasis_matrix >= 10, epistasis_matrix, "")

fig <- plot_ly(
    x = c("F639", "L640", "V641", "L642", "S643", "L644", "R645", "D646", "P647"),
    y = c("F639", "L640", "V641", "L642", "S643", "L644", "R645", "D646", "P647"),
    z = epistasis_matrix, type = "heatmap", zmin = 0, zmax = 20,
    text = text_matrix,
    texttemplate = "%{text}"
)%>%
    layout(title = paste(COMPOUND, " epistasis specific connectivity", sep = ""))

fig
# Save figure manually

### For a given position, how many mutants of that position lead to 
### positive epistasis ?

epistasis_per_position <- c()

for (k in 1:9){
    epistasis_per_position <- c(epistasis_per_position,
                               nrow(epis_data[epis_data$mut_pos_1 == k | epis_data$mut_pos_2 == k,]))
}

write(epistasis_per_position,
      file = paste("results/", COMPOUND, "_general_connectivity.txt", sep = ""), ncolumns = 10)
