# find_positive_epistasis.R
# Author : Mathieu Giguere
# Date : 4/12/2024
# Description : finds double mutants that lead to positive epistasis.

library(data.table)
source(file = "surface_at_gridpoint.R")
source(file = "plot_fitness_surface.R")
source(file = "quantile_fitness_surface_adaptive.R")

## Choose antifungal compound 
# {"anidulafungin", "caspofungin", "micafungin"}
COMPOUND <- "anidulafungin"

df <- read.csv(paste(COMPOUND, "_merged_final.csv", sep = ""))

L <- loess(s_double ~ s_1 + s_2, data = df, span = 0.2)

df2 <- data.table(df)

df2[,F_fit_loess := predict(L,newdata = .SD),,.SDcols = c("s_1","s_2")]

# calculate a loess-corrected fitness
df2[,fitness_norm := s_double-F_fit_loess]

### correct the loess approximated fitness surface by quantile surfaces
Nq = 100 # grid points along each axis
Nv = 500000 # max # of variants used in estimation of surface
span = max(c(0.01,500/nrow(df2))) # fraction of nearest neighbours to use for median calculation

## calculate quantile fitness surfaces
List = quantile_fitness_surface_adaptive(df2,Nq,Nv,span,Q=0.05)

df2[,F_fit_median := predict(List$F_median_fit,newdata = .SD) + F_fit_loess,,.SDcols = c("s_1","s_2")]
df2[,F_fit_lower := predict(List$F_lower_fit,newdata = .SD) + F_fit_loess,,.SDcols = c("s_1","s_2")]
df2[,F_fit_upper := predict(List$F_upper_fit,newdata = .SD) + F_fit_loess,,.SDcols = c("s_1","s_2")]

### PLOT FIG 1B

plot_fitness_surface(df2, L, List, dataset_dir = paste(getwd(), '/', sep = ''), prefix = paste(COMPOUND, '_', sep = ''))

### Now define data subset for positive epistasis.
## mark variants for positive epistasis analysis

#background_cutoff = df2[(s_1 + s_2) < lower_bound_F,
#                                quantile(s_double,probs=0.95,na.rm=T)]

df2[,pos_epistasis := FALSE]

df2[s_double > F_fit_upper, pos_epistasis := TRUE]

sum(df2$pos_epistasis, na.rm = TRUE)

## Save file

write.csv(df2, paste(COMPOUND, "_epistasis.csv", sep = ""))