################################################################################
### PUBH 6885 
### Fall 2022
### Team Project
### MTBLS326
### NMR Analysis Pipeline
### Matteo Sanchez-Dahl Gonzalez
### Zachary Siegfried 
################################################################################


################################################################################
### Set Up                                                                   ###
################################################################################

##########################################
### Packages                           ###
##########################################

### Bioconductor
# if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

### AlpsNMR
# BiocManager::install("AlpsNMR")
library(AlpsNMR)

# install.packages("devtools")
library(devtools)
install_github('sipss/AlpsNMR', force = TRUE)

### pracma
install.packages("pracma")
library(pracma)

##########################################
### Work Space                         ###
##########################################

### set working directory
setwd("/Users/matteosanchez-dahl/Desktop/F2022/PUBH 6885/team_project_MTBLS326/nmr/")

### breast cancer files
bc_dir <- "/Users/matteosanchez-dahl/Desktop/F2022/PUBH 6885/team_project_MTBLS326/nmr/data/bc"
bc_files <- fs::dir_ls(bc_dir)

### breast cancer directory
control_dir <- "/Users/matteosanchez-dahl/Desktop/F2022/PUBH 6885/team_project_MTBLS326/nmr/data/control"
control_files <- fs::dir_ls(control_dir)


################################################################################
### Data                                                                     ###
################################################################################

##########################################
### NMR Data                           ###
##########################################

### load nmr data
dataset <- nmr_read_samples(sample_names = c(bc_files, control_files))

##########################################
### Meta Data                          ###
##########################################

### load metadata
metadata_cohort  <- read.delim(
  "data/s_BC_IP3R_mod.txt",
  sep = '\t',
  header = T,
  fill = F,
  comment.char = "" ,
  check.names = F,
)

### clean metadata
metadata_cohort$NMRExperiment <- as.character(metadata_cohort$NMRExperiment)

### merge nmr spectra with metadata
dataset <- nmr_meta_add(
  dataset, 
  metadata = metadata_cohort, 
  by = "NMRExperiment"
)


################################################################################
### Functions                                                                ###
################################################################################

##########################################
### get_subset_NMRExps                ###
##########################################
### DESCRIPTION:
### function for obtaining list of the NMRExperiments corresponding to a 
### specified subset of samples in the study (e.g. "breast cancer"/"control")
###
### @param  subset_name  := string with sample status label (e.g. "control")
### @param  metadata     := metadata for study cohort that meets requirements
### @return experiment numbers corresponding to specified subset, if none, 
###         return all experiment names 
###
### REQUIREMENTS: 
### column labeled "Source_Name" with case vs. control status specified
### column labeled "NMRExperiment" with experiment number of sample specified 
### currently only support subsets "breast cancer" and "control "
###
### NOTE:
### currently only supports data and data from MTBLS326 dataset

get_subset_NMRExps <- function (metadata, subset_name = "all")
{
  if (strcmp(subset_name, "all")) {
    # group
    return (metadata$NMRExperiment)
  }
  samples <- c(which(metadata$Source_Name == subset_name))
  NMRExps <- metadata$NMRExperiment[samples]
  return (NMRExps)
}
##########################################


##########################################
### playful_ppm_plot                   ###
##########################################
### DESCRIPTION:
### function for obtaining list of the NMRExperiments corresponding to a 
### specified subset of samples in the study (e.g. "breast cancer"/"control")
###
### @param  data        := an AlpsNMR dataset of type [nmr_dataset_1D]
### @param  metabolite  := a metabolite of interest (e.g. "lactate)
### @param  subset      := subset to plot currently only support subsets 
###                        "breast cancer" and "control ", if none, plot all
### @return plot of NMR spectra for region of spectrum corresponding to a 
### metabolite of interest (e.g. "lactate"/"glucose"), if none, full spectrum
###
### REQUIREMENTS: 
### input dataset must be an AlpsNMR dataset of type [nmr_dataset_1D]
### 
### NOTE:
### currently only supports data and data from MTBLS326 dataset

playful_ppm_plot <- function(data, metabolite = NULL, subset = "all"){
  
  cs <- c(0.0,4.5)
  
  if (strcmp(metabolite, "TSP")){
    cs <- c(0.0,0.5)
  }
  else if (strcmp(metabolite, "lactate")){
    cs <- c(1.0,1.5)
  }
  else if (strcmp(metabolite, "glucose")){
    cs <- c(3.0,4.0)
  } 
  
  g <- plot(data, 
       interactive = TRUE,
       NMRExperiment = get_subset_NMRExps(metadata_cohort, subset),
       chemshift_range = cs)
  return(g)
} 
##########################################


################################################################################
### Work Flow                                                                ###
################################################################################

##########################################
### Subsets                            ###
##########################################

### control group
control_NMRExps <- get_subset_NMRExps(metadata_cohort, "control ")

### bc group
bc_NMRExps <- get_subset_NMRExps(metadata_cohort, "breast cancer")

### all
NMRExps <- get_subset_NMRExps(metadata_cohort, "all")

##########################################
### Interpolation                      ###
##########################################

### interpolation
dataset_int <- nmr_interpolate_1D(dataset, axis = c(min= 0.0, max = 4.5))

##########################################

### plot 

# controls
control_plot <- plot(dataset_int, 
                     NMRExperiment = control_NMRExps,
                     chemshift_range = c(0.0,4.5))
control_plot

# save png file
png(file="analysis/control_subset_spectra.png")
control_plot

# playful
playful_ppm_plot(dataset_int, metabolite = "", subset = "control ")

##########################################

# bc
bc_plot <- plot(dataset_int, 
                NMRExperiment = bc_NMRExps,
                chemshift_range = c(0.0,4.5))
bc_plot

# save png file
png(file="analysis/BC_subset_spectra.png")
bc_plot

# playful
playful_ppm_plot(dataset_int, metabolite = "", subset = "breast cancer")

##########################################
### Outlier Detection                  ###
##########################################

### outlier detection
pca_outliers_rob <- nmr_pca_outliers_robust(dataset_int, ncomp = 3)
outliers <- nmr_pca_outliers_plot(dataset_int, pca_outliers_rob)
outliers

# save png
png(file="analysis/outliers.png")
outliers

##########################################
### Peak Detection                     ###
##########################################

### peak detection
# calculate baseline threshold
baselineThresh <- nmr_baseline_threshold(dataset_int, 
                                         range_without_peaks = c(0.25,0.75))

# detect peaks
# SNR set to -1 for program to calculate it automatically
peak_list_initial <- nmr_detect_peaks(
  dataset_int,
  nDivRange_ppm = 0.1,
  scales = seq(1,16,2),
  baselineThresh = baselineThresh,
  SNR.Th = -1)

### plot
# detected peaks plot on sample
detect_peaks_plot <- nmr_detect_peaks_plot(dataset_int, 
                                           peak_list_initial, 
                                           NMRExperiment = "101",
                                           chemshift_range = c(3.0,4.0))

detect_peaks_plot

# save png
png(file="analysis/detect_peaks_plot_example.png")
detect_peaks_plot

##########################################
### Spectral Alignment                 ###
##########################################

### spectral alignment 
# find reference sample for alignment
ref_NMRExp <- nmr_align_find_ref(dataset_int, peak_list_initial)

# align spectra, allowing to for max shift of 0.05 ppm 
dataset_align <- nmr_align(
  nmr_dataset = dataset_int,
  peak_data = peak_list_initial,
  NMRExp_ref = ref_NMRExp,
  maxShift_ppm = 0.05,
  acceptLostPeak = TRUE
)
# maxShift_ppm = 0.0015

##########################################

### plot

## lactate peak

# controls
control_lac_spectra <-plot(dataset_int,
                           interactive = TRUE,
                           NMRExperiment = control_NMRExps,
                           chemshift_range = c(1.0,1.5))
control_lac_spectra

# save png
png(file="analysis/control_subset_lactate_spectra.png")
control_lac_spectra

# playful
playful_ppm_plot(dataset_int, 
                 metabolite = "lactate", 
                 subset = "control ")

##########################################

# bc
BC_lac_spectra <-plot(dataset_int,
                      interactive = TRUE,
                      NMRExperiment = bc_NMRExps,
                      chemshift_range = c(1.0,1.5))
BC_lac_spectra

# save png
png(file="analysis/BC_subset_lactate_spectra.png")
BC_lac_spectra

# playful
playful_ppm_plot(dataset_int, 
                 metabolite = "lactate", 
                 subset = "breast cancer")

##########################################

## TSP
# controls
control_TSP_spectra <-plot(dataset_int,
                           interactive = TRUE,
                           NMRExperiment = control_NMRExps,
                           chemshift_range = c(0.0,0.25))
control_TSP_spectra

# save png
png(file="analysis/control_subset_TSP_spectra.png")
control_TSP_spectra

# playful
playful_ppm_plot(dataset_int, 
                 metabolite = "TSP", 
                 subset = "control ")

##########################################

# bc
BC_TSP_spectra <-plot(dataset_int,
                      interactive = TRUE,
                      NMRExperiment = bc_NMRExps,
                      chemshift_range = c(0.0,0.25))
BC_TSP_spectra

# save png
png(file="analysis/BC_subset_TSP_spectra.png")
BC_TSP_spectra

# playful
playful_ppm_plot(dataset_int, 
                 metabolite = "TSP", 
                 subset = "breast cancer")

##########################################

### plot aligned 

## lactate peak

# controls
control_lac_aligned <-plot(dataset_align,
                           interactive = TRUE,
                           NMRExperiment = control_NMRExps,
                           chemshift_range = c(1.0,1.5))
control_lac_aligned

# save png
png(file="analysis/control_subset_lactate_aligned.png")
control_lac_aligned

# playful
playful_ppm_plot(dataset_align, 
                 metabolite = "lactate", 
                 subset = "control ")

##########################################

# bc
BC_lac_aligned <-plot(dataset_align,
                      interactive = TRUE,
                      NMRExperiment = bc_NMRExps,
                      chemshift_range = c(1.0,1.5))
BC_lac_aligned

# save png
png(file="analysis/BC_subset_lactate_aligned.png")
BC_lac_aligned

# playful
playful_ppm_plot(dataset_align, 
                 metabolite = "lactate", 
                 subset = "breast cancer")

##########################################

## TSP
# controls
control_TSP_aligned <-plot(dataset_align,
                           interactive = TRUE,
                           NMRExperiment = control_NMRExps,
                           chemshift_range = c(0.0,0.25))
control_TSP_aligned

# save png
png(file="analysis/control_subset_TSP_aligned.png")
control_TSP_aligned

# playful
playful_ppm_plot(dataset_align, 
                 metabolite = "TSP", 
                 subset = "control ")

##########################################

# bc
BC_TSP_aligned <-plot(dataset_align,
                      interactive = TRUE,
                      NMRExperiment = bc_NMRExps,
                      chemshift_range = c(0.0,0.25))
BC_TSP_aligned

# save png
png(file="analysis/BC_subset_TSP_aligned.png")
BC_TSP_aligned

# playful
playful_ppm_plot(dataset_align, 
                 metabolite = "TSP", 
                 subset = "breast cancer")

##########################################
### Data Normalization                 ###
##########################################

# area normalization
dataset_norm <- nmr_normalize(dataset_align,
                              method = "area")

##########################################
### Output NMR Data Matrix             ###
##########################################

### extract peak table 
peak_table <- nmr_data(dataset_norm)

### write to output directory
write.table(peak_table, file ="analysis/normalized_peakdata.txt", sep = "\t")

################################################################################
### EOF
################################################################################