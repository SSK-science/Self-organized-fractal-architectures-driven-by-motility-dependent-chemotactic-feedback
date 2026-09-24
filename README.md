# Self-organized fractal architectures driven by motility-dependent chemotactic feedback

This repository contains the code and data needed to reproduce the particle and chemical configurations reported in the "Self-organized fractal architectures driven by motility-dependent chemotactic feedback" article.
## Overview

The project studies self-organization patterns including fractal-like structures arising from motility-dependent chemotactic feedback in an agent based model. The code in this repository implements the simulation and clustering analysis. This produces the time series data of the particle and chemical configurations. This further is used for quantitative analysis.  

## Contents

Repository contents
- [Active_Chem_WCA.c](https://github.com/SSK-science/Self-organized-fractal-architectures-driven-by-motility-dependent-chemotactic-feedback/blob/main/code/Active_Chem_WCA.c)
  - Main simulation source code (C). Compile and run to generate simulation outputs.
- [Fig1.xlsx](https://github.com/SSK-science/Self-organized-fractal-architectures-driven-by-motility-dependent-chemotactic-feedback/blob/main/data/Fig1.xlsx)  - This contains the particle position data (x,y) and the chemical data plotted in Figure 1(b) and (c) respectively.
- [Fig2.xlsx](https://github.com/SSK-science/Self-organized-fractal-architectures-driven-by-motility-dependent-chemotactic-feedback/blob/main/data/Fig2.xlsx)  - Data for the correlation integral,c(r), as a function of distance, r, for DN, SC, HP, and RHP phases for different N, plotted in Fig 2, can be found here.    
- [Fig3.xlsx](https://github.com/SSK-science/Self-organized-fractal-architectures-driven-by-motility-dependent-chemotactic-feedback/blob/main/data/Fig3.xlsx)  - Contains data to obtain the phase diagram and fractal dimension (Fig 3(a)), the clustering data to obtain the cluster size distibution(column 3 of each block) in Fig 2(b,c), the mean cluster size and standard deviation Fig (3d) and the system span (Fig3e) are provided.
- [Fig4.xlsx](https://github.com/SSK-science/Self-organized-fractal-architectures-driven-by-motility-dependent-chemotactic-feedback/blob/main/data/Fig4.xlsx)  - The chemical data and the grid points to generate plots in Fig4(a-d, f-i) are given and the particle density and chemical density data to produce Fig4(e,j) are stored in this file.
- [Fig5.xlsx](https://github.com/SSK-science/Self-organized-fractal-architectures-driven-by-motility-dependent-chemotactic-feedback/blob/main/data/Fig5.xlsx)  - The mean squared displacement data in Fig 5(a-b) and the effective diffusion coeeficients scaledwith bare diffusion coefficient(without chemical) is provided here.
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
If you use this code in your research, please cite:

**Khuntia, S. S., Chaudhuri, D., & Chaudhuri, A. (2025).**
*Self-organized fractal architectures driven by motility-dependent chemotactic feedback.*
arXiv:2504.16539.
https://doi.org/10.48550/arXiv.2504.16539

The specific version of the simulation code and data associated with this work
is archived on Zenodo:

https://doi.org/10.5281/zenodo.22869440
