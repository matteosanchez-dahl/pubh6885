################################################################################
### PUBH 6885 
### Fall 2022
### Team Project
### MetaboLights Data Analysis
### MTBLS326
### Matteo Sanchez-Dahl Gonzalez
### Zachary Siegfried 
################################################################################


################################################################################
### Set Up                                                                   ###
################################################################################

##########################################
### Packages                           ###
##########################################

library(vegan)
library(devtools)
library(omicsArt)
library(mgcv)
library(nlme)
library(labdsv)
library(Maaslin2)
library(readxl)
#install_github("vqv/ggbiplot")
library(ggbiplot)

##########################################
### Work Space                         ###
##########################################

setwd("/Users/matteosanchez-dahl/Desktop/F2022/PUBH 6885/team_project_MTBLS326/MTBLS326")


################################################################################
### Data                                                                     ###
################################################################################

##########################################
### Metabolomic Data                   ###
##########################################

### metabolomic from MTBLS326
BC_metabolites  <- read.delim(
  "data_MTBLS326/dataNMR_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

### normalized metabolomics data 
BC_metabolites_normalized  <- read.delim(
  "data_MTBLS326/dataNMR_normalized_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

### patients with high ip3r expression & controls
BC_metabolites_normalized_high <- read.delim(
  "data_MTBLS326/dataNMR_normalized_high_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

### patients with high ip3r expression & controls
BC_metabolites_normalized_low <- read.delim(
  "data_MTBLS326/dataNMR_normalized_low_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)

##########################################
### Meta Data                          ###
##########################################

### metadata from MTBLS326
BC_metadata  <- read.delim(
  "data_MTBLS326/metadata_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F
  #,row.names = 1
)

### metadata of whole cohort with generalized labels (n=42)
### "breast cancer", "control "
BC_general_metadata  <- read.delim(
  "data_MTBLS326/metadata_general_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F
  #,row.names = 1
)

### metadata of high ip3r expression subset 
### "high", "control"
BC_high_metadata <- read.delim(
  "data_MTBLS326/metadata_high_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F
  #, row.names = 1
)

### metadata of low ip3 expression subset 
### "low", "high"
BC_low_metadata <- read.delim(
  "data_MTBLS326/metadata_low_BC_IP3R.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
  row.names = 1
)


################################################################################
### Work Flow                                                                ###
################################################################################

##########################################
### Maaslin2                           ###
##########################################

Maaslin2::Maaslin2(
  input_data  = BC_metabolites_normalized_high,
  input_metadata =  BC_high_metadata,
  output =  'analysis_MTBLS326/maaslin2_normalized_high_output',
  max_significance= 0.1, # target FDR
  analysis_method = 'LM', # linear model
  standardize = FALSE,
  transform = 'LOG',
  normalization = 'NONE',
  heatmap_first_n = 20,
  reference = c("Group,control")) 

##########################################
### Tweedieverse                       ###
##########################################

Tweedieverse::Tweedieverse(
  BC_metabolites_normalized_high,
  BC_high_metadata,
  output =  'analysis_MTBLS326/tweedieverse_output',
  heatmap_first_n = 20,
  reference = c("Group,control")) 

##########################################
### Principal Component Analysis       ###
##########################################

### all groups -- not normalized 
pca <- prcomp(BC_metabolites, 
              center = TRUE, 
              scale. = FALSE)
ggbiplot(pca)

### total cohort w norm
### no scaling to do PCA on covariance
pca_norm <- prcomp(BC_metabolites_normalized, 
                   center = TRUE, 
                   scale. = FALSE)

### case vs. control
pca_norm_p_bc <- ggbiplot(pca_norm, 
                          ellipse=TRUE, 
                          groups=BC_general_metadata$Group)
# view
pca_norm_p_bc

# download png
png(file="analysis_MTBLS326/PCA_BC_control.png",
    width=600, height=350)
    pca_norm_p_bc

### high, low, control
pca_norm_p_groups <- ggbiplot(pca_norm, 
                              ellipse=TRUE, 
                              groups=BC_metadata$Group)
# view
pca_norm_p_groups

# download png
png(file="analysis_MTBLS326/PCA_groups_control.png",
    width=600, height=350) 
    pca_norm_p_groups

### high ip3r vs. control
pca_norm_high <- prcomp(BC_metabolites_normalized_high, 
                        center = TRUE, 
                        scale. = FALSE)
# plot
pca_norm_p_high <- ggbiplot(pca_norm_high, 
                            ellipse=TRUE, 
                            groups=BC_high_metadata$Group)
# view
pca_norm_p_high

# download png
png(file="analysis_MTBLS326/PCA_high_control.png",
    width=600, height=350) 
    pca_norm_p_high

### low ip3r vs. control
pca_norm_low <- prcomp(BC_metabolites_normalized_low, 
                       center = TRUE, 
                       scale. = FALSE)
# plot
pca_norm_p_low <- ggbiplot(pca_norm_low, 
                           ellipse=TRUE, 
                           groups=BC_low_metadata$Group)
# view
pca_norm_p_low

# download png
png(file="analysis_MTBLS326/PCA_low_control.png",
    width=600, height=350)
    pca_norm_p_low

################################################################################
### EOF
################################################################################