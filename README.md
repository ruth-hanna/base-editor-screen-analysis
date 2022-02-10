# base-editor-screen-analysis

This repo contains Jupyter notebooks used to generate the figures in "Massively parallel assessment of human variants with base editor screens." 

The starting files are Tables S1, S3, S6, S7, S9, and S10, provided with the manuscript. To generate filtered log-fold change (LFC) files, run the Screening_Data_Preprocessing notebook; this notebook reads in the Excel files and outputs the filtered LFC files used for subsequent analyses.

Additional annotation files for the ClinVar analysis and filtered starting files are provided here:
https://figshare.com/articles/dataset/Filtered_LFC_files/19140497
https://figshare.com/articles/dataset/Annotations/19140482

The notebooks required to analyze the validation experiments are provided in the "validation-experiments" folder. These data are deposited in the Sequencing Reads Archive at the following link:
https://www.ncbi.nlm.nih.gov/sra/PRJNA638656

In order to design and annotate sgRNAs for base editing screens, please use the base editor design tool located here: 
https://github.com/mhegde/base-editor-design-tool
