# Large files can be downloaded from
# https://github.com/sandrejev/growth_curves/releases/download/v1.0/data.zip

install.packages(c("ggplot2", "plyr", "reshape2", "gplots", "dplyr", "tidyr", "readr", "beeswarm", "XLConnect", "RSQLite", "VennDiagram", "RColorBrewer", "pheatmap", "gridExtra", "XLConnectJars", "ape"))
dir.create("../report", showWarnings = FALSE)
dir.create("../data/supplementary", showWarnings = FALSE)
dir.create("../data/generated", showWarnings = FALSE)

source("analysis.media_species_selection.R")
unzip("../data/kegg.zip", exdir="../data")
piechart.media_summary() # !
piechart.species_summary()
cummulative_plot.enzymatic_coverage()

source("analysis.overview.R")
table.growth_matrix()
table.replicates_number() # Generate replicates table
heatmap.obsulute() # !
heatmap.relative() # !
#heatmap.curves()
itol.phylogenetic_tree()
plots.media_preferrence()

source("analysis.abundance.R")
boxplots.prevalence_correlation()

source("analyze.pH.R")
boxplots_and_scater.pH_overview()

source("analysis.phylogenetic.R")
scatter.abundance_correlations()
taxonomicrank_growth_barplots() # does not exist

source("analysis.mucin.R")
plot.mucin_validation_data()
