# **PUBH6885 Team Project** 🧬 🩸 🧲 📊
This repo contains functions to analyze 1H NMR spectral data from the [MTBLS326](https://www.ebi.ac.uk/metabolights/MTBLS326/files) study on MetaboLights EMBL-EBI database.

📖[Singh, A. et al. (2017) 1H NMR Metabolomics Reveals Association of High Expression of Inositol 1, 4, 5 Trisphosphate Receptor and Metabolites in Breast Cancer Patients. PLOS ONE, 12, e0169330.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169330#pone-0169330-t001)

## **Functions**

### **get_subset_NMRExps**
```
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
```

### **playful_ppm_plot**
```
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

playful_ppm_plot <- function(data, metabolite = NULL, subset = "all")
```



