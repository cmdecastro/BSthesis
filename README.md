# Extracting Philippine voting patterns through hyperspectral unmixing

This repository contains all the codes and some data I used for my undergraduate thesis on Philippine Senate election patterns from 2013-2019 using Hyperspectral Unmixing.
The following is a short decription of what each folder contains:

### 1. MATLAB 
This folder contains the codes for MVSA, VCA, and SuNSAL provided in the paper:

J. Li, A. Agathos, D. Zaharie, J. M. Bioucas-Dias, A. Plaza, and X. Li. Minimum volume simplex analysis: A fast algorithm for linear hyperspectral unmixing. 
IEEE Transactions on Geoscience and Remote Sensing, 53(9):5067-5082, Sep. 2015.

I used these codes to perform hyperspectral unmixing.

### 2. archetypes
This folder contains the codes to read, analyze, and plot the outputs from the MATLAB codes after hyperspectral unmixing. 
The Jupyter Notebooks in this folder also contain the codes I used to clean my datasets.

### 3. election_data
This folder contains some of the senate election data I used in my thesis. I wasn't able to include some because the file sizes were too big. 
They can be found in this link: http://elections.org.ph/.

### 4. kmeans_clustering
This folder contains all the code I used to perform K-mean Clustering on the election data.

### 5. province_shapefile
This folder contains the shapefiles I used to plot a Philippine map with provincial boundaries.

### 6. region_shapefile
This folder contains the shapefiles I used to plot a Philippine map with regional boundaries.

### 7. unmixing_results
This folder contains the output .csv files for archetypes and weights obtained after hyperspectral unmixing.

### 8. DE CASTRO, CMM 2.19-20_Extracting Philippine voting patterns through hyperspectral unmixing.pdf
This is my full submitted manuscript.
