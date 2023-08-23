# p53-prime-editing-sensor

## Created by Sam Gould (samgould@mit.edu) | Sánchez-Rivera Lab | Koch Institute (MIT)

Code for analysis of the p53 prime editing sensor screening dataset.

This repository is broken down into the following sections:
1. Jupyter Notebooks for recreating all of the figures that appear in the paper, included extended data figures.
    - These notebooks are located in the main directory and are labelled according to the figure #. All extended data figures are generated in a single notebook.
    - All source data is referenced with relative paths to allow for complete reproducibility (you should be able to download the repository and directly run the ) 

2. Analysis scripts used to process raw fastq files (**see folder: analysis_scripts**):
    - **split_fastq.py**: python script for deconvolution of fastq file into separate files by sample.
    - **sensor_extraction.py**: python script for filtering reads, counting pegRNAs and recombination events, extracting sensor reads and placing them in separate fastq files for each pegRNA.
    - **crispresso_analysis.py**: python script for using crispresso2 to perform quantification of editing outcomes for each pegRNA at each time-point/condition/replicate.
        - **crispresso_environment.yml**: for generating a virtual environment that allows users to run crispresso via python commands (this script will break otherwise).
    - **crispresso_analysis_aggregation.py**: as the name suggests, a python script for aggregating the results of crispresso runs (from the 28,000+ folders per sample) and formatting them into readable csv files.
    - Other files in this folder are referenced by some (or all scripts).
    - The .sh scripts with corresponding names are the scripts that were used to parallelize and run these jobs on the Luria Computing Cluster at MIT.

3. Source data for generating figures
    - These are contained in the other folders, and are referenced (as relative paths) by the figure generation python notebooks.

##
The jupyter notebooks for recreating the figures that appear in the paper should run in most python environments (no particularly exotic packages are used). For troubleshooting purposes, here is a list of packages and version numbers used:

python==3.9.12

matplotlib==3.6.2

numpy==1.23.5

pandas==1.5.2

regex==2022.10.31

scikit-learn==0.24.2

scipy==1.10.0

seaborn==0.12.2

statannot==0.2.3