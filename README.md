# Recovering Spatially-Varying Cell-Specific Gene Co-Expression Networks for Single-Cell Spatial Expression Data.

# Jinge Yu

## Data

### Abstract 

We use MERFISH mouse hypothalamus data and MERFISH U-2 OS cell line datasets for our analysis, refer to Moffitt JR, Bambah-Mukku D, Eichhorn SW, Vaughn E, Shekhar K, Perez JD, et al. Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region. Science 362, 2018 and Xia C, Fan J, Emanuel G, et al. Spatial transcriptome profiling by MERFISH reveals subcellular RNA compartmentalization and cell cycle-dependent gene expression[J]. Proceedings of the National Academy of Sciences, 2019, 116(39): 19490-19499 for details.

### Availability 

-	The MERFISH mouse hypothalamus data is publicly available for download via the online data portal at https://datadryad.org/stash/dataset/doi:10.5061/dryad.8t8s248. No registration is required!
-	The MERFISH U-2 OS cell line datasets are publicly available for download via the online data portal at https: //www.pnas.org/content/116/39/19490/tab-figures-data. No registration is required.



## Code

### Abstract

All of the data preprocessing and analysis in this paper were completed using R. The code is provided to conduct preprocessing on the raw data, apply proposed model via two-step algorithm, WGCNA, CTS, CSN-joint and CSN-separate, and generate descriptive plots.

### Description

All of the R scripts used in the report are available in a public repository on GitHub [https://github.com/jingeyu/CSSN_data_code]. The MIT license applies to all code, and no permissions are required to access the code.

### Optional Information

R version 4.0.3 (2020-10-10, “Bunny-Wunnies Freak Out”)
Copyright (C) 2020 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin17.0 (64-bit) were used in the project. And The necessary R libraries for the code used for data processing and analysis are:

-	dplyr, version 0.8.5 (https://CRAN.R-project.org/package=dplyr)
-	ggplot2, version 3.3.2 (https://CRAN.R-project.org/package=ggplot2)
-	pheatmap, version 1.0.12 (https://CRAN.R-project.org/package=pheatmap)
-	R.matlab, version 3.6.2 (https://CRAN.R-project.org/package=R.matlab)
-	dplyr, version 1.0.2 (https://CRAN.R-project.org/package=dplyr)
-	MASS, version 7.3-53 (https://CRAN.R-project.org/package=MASS)
-	Matrix, version 1.2-18 (https://CRAN.R-project.org/package=Matrix)
-	CholWishart, version 1.1.0 (https://CRAN.R-project.org/package=CholWishart)
-	stats, version 4.0.3 (https://CRAN.R-project.org/package=stats)
-	Seurat, version 4.0.0 (https://CRAN.R-project.org/package=Seurat)
- WGCNA, version 1.69 
(http://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/)



Computer information:
- Operating system: macOS 10.15.5
- CPU: Quad-Core Intel Core i5 2.4G Hz. 


## Instructions for Use

### Reproducibility

All data preprocessing and analysis as well as Table 1, Figure 2—Figure 9 in the manuscript can be reproduced.

Detailed workflow information is contained in the "Reproducibility in the Simulation.docx" in "Simulation" and "Reproducibility in the Real Application.docx" in "Real Application" directories. 

The general steps in the simulation are:
 1. Generate the data and apply the proposed model to the data.
2. Generate Figure 2—Figure 4 and Table 1 in the manuscript.

The general steps in the real application are:
 1. Conduct data preprocessing.
 2. Apply the proposed model to the preprocessed data. 
 3. Generate the Figure 5—Figure 9 in the paper.
![Uploading image.png…]()

