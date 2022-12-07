library(vegan)
library(devtools)
library(omicsArt)
library(mgcv)
library(nlme)
library(labdsv)
library(Maaslin2)
library(readxl)

setwd("/Users/matteosanchez-dahl/Desktop/F2022/PUBH 6885/group/")

BC_metabolites  <- read.delim(
  "data_MTBLS326/dataNMR_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

BC_metadata  <- read.delim(
  "data_MTBLS326/metadata_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

BC_metabolites_normalized  <- read.delim(
  "data_MTBLS326/dataNMR_normalized_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

BC_metabolites_normalized_high <- read.delim(
  "data_MTBLS326/dataNMR_normalized_high_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

BC_metabolites_normalized_low <- read.delim(
  "data_MTBLS326/dataNMR_normalized_low_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

BC_general_metadata  <- read.delim(
  "data_MTBLS326/metadata_general_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F
  #,row.names = 1
)

BC_high_metadata <- read.delim(
  "data_MTBLS326/metadata_high_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

BC_low_metadata <- read.delim(
  "data_MTBLS326/metadata_low_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

Maaslin2::Maaslin2(
  input_data  = BC_metabolites,
  input_metadata =  BC_metadata,
  output =  'analysis/maaslin2_output',
  max_significance= 0.1, # target FDR
  analysis_method = 'LM', # linear model
  standardize = FALSE,
  transform = 'LOG',
  normalization = 'NONE',
  heatmap_first_n = 20,
  reference = c("Group,control")) 

Maaslin2::Maaslin2(
  input_data  = BC_metabolites,
  input_metadata =  BC_general_metadata,
  output =  'analysis/maaslin2__general_output',
  max_significance= 0.1, # target FDR
  analysis_method = 'LM', # linear model
  standardize = FALSE,
  transform = 'LOG',
  normalization = 'NONE',
  heatmap_first_n = 20,
  reference = c("Group,control")) 

Maaslin2::Maaslin2(
  input_data  = BC_metabolites_normalized,
  input_metadata =  BC_general_metadata,
  output =  'analysis/maaslin2_normalized_output',
  max_significance= 0.1, # target FDR
  analysis_method = 'LM', # linear model
  standardize = FALSE,
  transform = 'LOG',
  normalization = 'NONE',
  heatmap_first_n = 20,
  reference = c("Group,control")) 


Maaslin2::Maaslin2(
  input_data  = BC_metabolites_normalized_high,
  input_metadata =  BC_high_metadata,
  output =  'analysis/maaslin2_normalized_high_output',
  max_significance= 0.1, # target FDR
  analysis_method = 'LM', # linear model
  standardize = FALSE,
  transform = 'LOG',
  normalization = 'NONE',
  heatmap_first_n = 20,
  reference = c("Group,control")) 


Tweedieverse::Tweedieverse(
  BC_metabolites_normalized_high,
  BC_high_metadata,
  output =  'analysis/tweedieverse_output',
  heatmap_first_n = 20,
  reference = c("Group,control")) 



library(devtools)

install_github("vqv/ggbiplot")

library(ggbiplot)

# total cohort w norm
pca_norm <- prcomp(BC_metabolites_normalized, center = TRUE, scale. = TRUE)

# case vs. control
ggbiplot(pca_norm, ellipse=TRUE, groups=BC_general_metadata$Group)

# high, low, control
ggbiplot(pca_norm, ellipse=TRUE, groups=BC_metadata$Group)

# high ip3r 
pca_norm_high <- prcomp(BC_metabolites_normalized_high, center = TRUE, scale. = TRUE)
ggbiplot(pca_norm_high, ellipse=TRUE, groups=BC_high_metadata$Group)

# low ip3r
pca_norm_low <- prcomp(BC_metabolites_normalized_low, center = TRUE, scale. = TRUE)
ggbiplot(pca_norm_low, ellipse=TRUE, groups=BC_low_metadata$Group)

# all groups 
pca <- prcomp(BC_metabolites, center = TRUE, scale. = TRUE)
ggbiplot(pca)

pcoa_plots <- omicsArt::ordplots(BC_metabolites_normalized, BC_general_metadata, output = 'analysis/', outputname = 'pcoa_general', method = 'pcoa')











