# Self-organized fractal architectures driven by motility-dependent chemotactic feedback

This repository contains the code and data needed to reproduce the particle and chemical configurations reported in the "Self-organized fractal architectures driven by motility-dependent chemotactic feedback" article.
## Overview

The project studies self-organization patterns including fractal-like structures arising from motility-dependent chemotactic feedback in an agent based model. The code in this repository implements the simulation and clustering analysis. This produces the time series data of the particle and chemical configurations. This further is used for quantitative analysis.  

## Contents

Repository contents
- [Active_Chem_WCA.c](https://github.com/SSK-science/Self-organized-fractal-architectures-driven-by-motility-dependent-chemotactic-feedback/blob/main/Active_Chem_WCA.c)
  - Main simulation source code (C). Compile and run to generate simulation outputs.
- [Fig1.xlsx](https://github.com/SSK-science/Self-organized-fractal-architectures-driven-by-motility-dependent-chemotactic-feedback/blob/main/Fig1.xlsx)  - This contains the particle position data (x,y) and the chemical data plotted in Figure 1(b) and (c) respectively.  
- [Fig2.xlsx](https://github.com/SSK-science/Self-organized-fractal-architectures-driven-by-motility-dependent-chemotactic-feedback/blob/main/Fig2.xlsx)  - Contains data to obtain the phase diagram and fractal dimension (Fig 2(a)), the clustering data to obtain the cluster size distibution(column 3 of each block) in Fig 2(b,c), the mean cluster size and standard deviation Fig (2d) and the system span (Fig2e) are provided.
- [Fig3.xlsx](https://github.com/SSK-science/Self-organized-fractal-architectures-driven-by-motility-dependent-chemotactic-feedback/blob/main/Fig3.xlsx)  - The chemical data and the grid points to generate plots in Fig3(a-d, f-i) are given and the particle density and chemical density data to produce Fig3(e,j) are stored in this file.
- [Fig4.xlsx](https://github.com/SSK-science/Self-organized-fractal-architectures-driven-by-motility-dependent-chemotactic-feedback/blob/main/Fig4.xlsx)  - The mean squared displacement data in Fig 4(a-b) and the effective diffusion coeeficients scaledwith bare diffusion coefficient(without chemical) is provided here.
-README.md This file

Build and run (suggested)
- Compile (example):
  - gcc -o active_chem_WCA Active_Chem_WCA.c -lm
- Run (example):
  - ./active_chem_WCA  
- `README.md` — This file.

## Contact

For questions about the code or reproducing figures, contact the repository owner:
- GitHub: [SSK-science](https://github.com/SSK-science)
- Or open an issue in this repository.

## Cite
Please cite  our work ["Self-organized fractal architectures driven by motility-dependent chemotactic feedback"] (https://doi.org/10.48550/arXiv.2504.16539) if you use our code or the data provided in this repository.
