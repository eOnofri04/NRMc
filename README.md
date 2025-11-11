# A General Multi-Stratum Model for a Nanofunctionalized Releasing Capsule  
**Code Repository for the Paper:**  
[![DOI](https://img.shields.io/badge/DOI-10.3934/mbe.2025xxx-blue)](https://doi.org/10.3934/mbe.2025xxx)
[![Arxiv](http://img.shields.io/badge/arXiv-2506.17078-b31a1a)](https://arxiv.org/abs/2506.17078)  
*E. Onofri, E. Cristiani, A. Martelli, P. Gentile, J. Girón Hernández, G. Pontrelli*  
**"A general multi-stratum model for a nanofunctionalized releasing capsule: An experiment-driven computational study"**  
*Mathematical Biosciences and Engineering (MBE), 2025*  

---

## Overview

This repository contains the MATLAB implementation of the numerical framework described in the above paper.  
The core script, `diffusionMultiscaleSpheric.m`, implements the finite-volume solver for a multi-stratum spherical diffusion model, supporting spatial and temporal adaptivity, radial asymmetry, and multi-layer coupling.

The repository also provides the configuration files, experimental data, and plotting routines used to generate the figures and results presented in the publication.
Experimental results are pre-evaluated for the reader's convenience, yet can be produced again and compared with the one published via `compareMatFiles.m` below.

_Note._ Please notice that some results (particularly the ones regarding sensitivity analysis) where created using a previous version of the code, hence the structure of the results might vary; yet numerical values are consistent across versions.

---

## Repository Structure

```
├── diffusionMultiscaleSpheric.m   # Main simulation script
├── pathMBE.m                      # Adds repository path to MATLAB
├── utils/
│   └── compareMatFiles.m          # Utility to compare two .mat files (AI-generated)
│
├── test_academic/
│   ├── *.m, *.mat                 # Configurations and results for academic experiments
│   └── Plots.mlx                  # Script to generate plots from academic tests
│
├── test_real/
│   ├── exp_2412.m/.mat            # Real-data experiment with radial asymmetry
│   ├── exp_2412_SD.m/.mat         # Real-data experiment without radial asymmetry
│   └── Plots_real.m/_nb.mlx       # Plotting script for real-data tests
│
├── sensitivity_real/
│   ├── sensitivityPlot.m/_nb.mlx  # Sensitivity analysis plotting script
│   └── *.mat                      # Results of sensitivity tests
│
├── img/                           # Output folder for generated images
│
├── LICENSE                        # CC-BY 4.0 license
└── README.md                      # This file
```

---

## Usage

The main script `diffusionMultiscaleSpheric.m` can be executed directly in MATLAB after adding the repository to the MATLAB path via:

```matlab
pathMBE
```

The simulation behaviour is controlled by a set of variables that can be defined before running the script:

| Variable        | Description                                                     | Default                |
|-----------------|-----------------------------------------------------------------|------------------------|
| `exp_config_fn` | Path to the configuration file for the simulation               | `'test_real/exp_2412'` |
| `do_plot`       | Plot the simulation during execution (slower)                   | `false`                |
| `do_mov`        | Save an animation of the simulation (requires `do_plot = true`) | `false`                |
| `do_pause`      | Pause after each plotting step                                  | `false`                |
| `do_plot_finer` | Update plot after every computational step                      | `false`                |
| `do_export`     | Export each frame as a `.png` image                             | `false`                |
| `do_use_exp`    | Include experimental data as reference in plots                 | `false`                |
| `do_save_U`     | Save full time evolution of the concentration field             | `false`                |

After the simulation, results are automatically saved to a `.mat` file matching the configuration name, with suffix `_U` if `do_save_U` is enabled.

---

## Experiments
 - **test_academic/**
 Contains benchmark and validation tests (single/multi-layer, coarse/fine resolution, varying time steps).
 Results are provided as `.mat` files and can be visualised using `Plots.mlx`.
 - **test_real/**
 Includes simulations based on real experimental data with and without radial asymmetry.
 Plots can be reproduced via `Plots_real.mlx`.
 - **sensitivity_real/**
 Provides sensitivity analysis with respect to selected parameters.
 Corresponding plots are available through `sensitivityPlot.mlx`.

---

## License

This project is released under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
You are free to share and adapt the material, provided that appropriate credit is given to the authors and source publication.

See the LICENSE￼file for full terms.

---

## Citation

If you use this code or any part of it in your work, please cite:
```
E. Onofri, E. Cristiani, A. Martelli, P. Gentile, J. Girón Hernández, G. Pontrelli,
A general multi-stratum model for a nanofunctionalized releasing capsule: An experiment-driven computational study,
Mathematical Biosciences and Engineering (MBE), 2025.
DOI: 10.3934/mbe.2025xxx
```

---

## Disclaimer

All source code contained in this repository was developed by Elia Onofri.  
For any enquiries, requests for clarification, or technical issues, please contact me at [**elia.onofri@cnr.it**](mailto:elia.onofri@cnr.it).
