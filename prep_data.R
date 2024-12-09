# prep_data.R
# Author : Mathieu G.
# Date : 03/12/2024
# Description : Prepares the dataframes required for further analysis.

## Choose antifungal compound 
# {"anidulafungin", "caspofungin", "micafungin"}
COMPOUND <- 'anidulafungin'

## Import data.
dfdouble <- read.csv("double_aa_classification.csv")
dfsingle <- read.csv("single_aa_classification.csv")

## Separate single mutant and double mutant.

singles <- dfsingle[dfsingle$seq_type == 'single',]

keeps <- c("X", "compound", "aa_seq", "s")

singles <- singles[keeps]
doubles <- dfdouble[keeps]

## Separate data according to the chosen compound.

c_singles <- singles[singles$compound == COMPOUND,]
c_doubles <- doubles[doubles$compound == COMPOUND,]


## Standardize the fitness scores in each subset dataframe.
standardize <- function(input_df) {
    input_df['s'] <- as.data.frame(scale(input_df['s']))
    input_df
}
c_singles <- standardize(c_singles)
c_doubles <- standardize(c_doubles)

## Create df to plot single and double mutants.
## Find mutation positions and alt_aa for all mutants.

WT <- 'FLVLSLRDP'

c_singles['mut_pos_1'] <- mapply(function(x, y) which(x != y)[1], 
                                strsplit(c_singles$aa_seq, ""),
                                strsplit('FLVLSLRDP', ""))

c_singles$alt_aa_1 <- mapply(function(str, idx) substr(str, idx, idx),
                               c_singles$aa_seq, c_singles$mut_pos_1)

c_doubles['mut_pos_1'] <- mapply(function(x, y) which(x != y)[1], 
                                   strsplit(c_doubles$aa_seq, ""),
                                   strsplit('FLVLSLRDP', ""))

c_doubles['mut_pos_2'] <- mapply(function(x, y) which(x != y)[2], 
                                   strsplit(c_doubles$aa_seq, ""),
                                   strsplit('FLVLSLRDP', ""))

c_doubles$alt_aa_1 <- mapply(function(str, idx) substr(str, idx, idx),
                               c_doubles$aa_seq, c_doubles$mut_pos_1)

c_doubles$alt_aa_2 <- mapply(function(str, idx) substr(str, idx, idx),
                               c_doubles$aa_seq, c_doubles$mut_pos_2)

# rename column s to s_double
names(c_doubles)[names(c_doubles) == 's'] <- 's_double'
names(c_doubles)[names(c_doubles) == 'aa_seq'] <- 'aa_seq_double'

merged_1 <- merge(c_doubles, c_singles[c('mut_pos_1', 'alt_aa_1', 's')],
                      by.x = c('mut_pos_1', 'alt_aa_1'),
                      by.y= c('mut_pos_1', 'alt_aa_1'))

names(merged_1)[names(merged_1) == 's'] <- 's_1'

merged_2 <- merge(merged_1, c_singles[c('mut_pos_1', 'alt_aa_1', 's')],
                      by.x = c('mut_pos_2', 'alt_aa_2'),
                      by.y= c('mut_pos_1', 'alt_aa_1'))

names(merged_2)[names(merged_2) == 's'] <- 's_2'


merged_final <- merged_2[, c('compound', 'aa_seq_double', 'mut_pos_1',
                                   'alt_aa_1', 'mut_pos_2', 'alt_aa_2',
                                   's_double', 's_1', 's_2')]

# Remvoe stop codons
merged_final <- merged_final[merged_final$alt_aa_1 != '*',]
merged_final <- merged_final[merged_final$alt_aa_2 != '*',]

write.csv(merged_final, paste(COMPOUND, "_merged_final.csv", sep = ""))

### 3D plot

library(plotly)

fig <- plot_ly(merged_final, x = ~s_1, y = ~s_2, z = ~s_double,
               marker = list(color = ~s_double, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE))
fig <- fig %>% add_markers()
fig <- fig %>% layout(title = paste(COMPOUND, " fitness scores", sep = ""),
                      scene = list(xaxis = list(title = 's_1'),
                                   yaxis = list(title = 's_2'),
                                   zaxis = list(title = 's_double')),
                      annotations = list(
                          x = 1.13,
                          y = 1.05,
                          text = 'fitness score',
                          xref = 'fitness score',
                          yref = 'fitness score',
                          showarrow = FALSE
                      ))
fig