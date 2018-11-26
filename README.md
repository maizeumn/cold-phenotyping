https://doi.org/10.5281/zenodo.1553411
[![DOI](https://zenodo.org/badge/151628212.svg)](https://zenodo.org/badge/latestdoi/151628212)

# cold-phenotyping

### 1. QR code sheet generation

"generate_qrcode_R_script.pl" creates an R script to generate .pdf versions of our QR code sheets to embed metadata in filenames and track days after sowing.

### 2. Image acquisition

Shell and matlab scripts in this directory were used to acquire images.

Image acquisition is initiated by running "takeImages.sh". This script calls functions and other scripts included in this directory, and includes target values for camera angle, focalTarget, focalThresh, fNumberTarget, fNumberThresh, ExposureTarget, ExposureThresh.

### 3. Image analysis with PlantCV

"PlantCV_pipeline.py" is the python code used with PlantCV v3.0dev2 to analyze images. It was run using the plantcv-pipeline.py function, and the included Probability Density Function file ("pdfs.txt") for pixel classification.

Latest PlantCV documentation is here: https://plantcv.readthedocs.io/en/latest/

Link to PlantCV v3.0dev2 github page here: https://github.com/danforthcenter/plantcv/releases/tag/v3.0.dev2

### 4. Output data from PlantCV analysis

"PlantCV_output.zip" is a compressed version of the merged csv files from the PlantCV analysis for images here: [cyverse doi].
It contains 2 .csv files:
   1) Output_genotype_survey.csv
   2) Output_stress_duration_optimization.csv

### 5. Data analysis with R

"Figures_Endersetal2018.R" is an R script used to create figures analyzing the data in the csv files in "PlantCV_output.zip".
