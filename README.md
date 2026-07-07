# radio-ecle
Code for the analysis of radio (and other multi-wavelength) observations of extreme coronal line emitters. If you
make use of this code please cite Franz et al. (2026), "The Radio Properties of Extreme Coronal Line Emitters: 
Constraints on the Sub-parsec Environment" (to be posted on arxiv any day now!). A brief description of the files and 
directories in this repo are below.

The code in this repo is very dependent on the OTTER data ecosystem and you will need to install the OTTER API:
[OTTER](https://github.com/astro-otter). To run the fitting code, you will also need to install 
and [syncfit](https://github.com/alexander-group/syncfit).

## File & Directory Descriptions
### Data & Other Directories

1. `data`: The raw datasets used in this work. This includes scripts to do forced photometry in CASA on the radio images. (Only available on zenodo)
2. `private-data`: The OTTER formatted CSV files of the photometry after we've cleaned it up.
3. `convergence-plots`: Plots generated in the `synchrotron-modeling.ipynb` notebook automatically to assess convergence and reliability of fit results.

### Python Code
* compile-data.ipynb: Converts and compiles the raw, multiwavelength datasets into the OTTER format. 
* frac.ipynb: Calculates some redshift limits for different RMSs for our "fraction of radio bright ECLE" section
* luminosity.ipynb: Generates fig 1 and 2 in the paper
* plot-all-lcs.ipynb: Plots the multi-wavelength light curves, this was for exploratory purposes.
* radio-sfr.ipynb: Calculates some radio-derived SFRs from our observations to compare to spectral lines.
* targeted-radio-followup-table.ipynb: Generates a table of targeted radio observations for the paper appendix
* Modeling
  * synchrotron-modeling-free-smoothing.ipynb: A test of our synchrotron SED model if we free the smoothing parameter.
  * synchrotron-modeling.ipynb: The code actually used in the paper for synchrotron SED modeling.
  * thermal-synchrotron-modeling.ipynb: Code to try thermal synchrotron modeling of the SEDs.
  * sed-plot.ipynb: Plot the radio SEDs with the best fit parameters
* Equipartition
  * equipartition.py: Functions with the equipartition code
  * equipartition-analysis.ipynb: Equipartition analysis code to generate fig 4 and 5. Calls equipartition.py. 

### Other data files
* `*chains.json`: The output chaings from the MCMC fitting
* `all-photometry.json`: A CSV file in the OTTER format of all of the photometry
* `ecle-metadata.csv`: The metadata for the ECLEs in a machine readable format
* `ECLE-metadata-latex.csv`: The LaTeX table with the metadata for the ECLEs
* `lit_sed_params.csv`, `ASASSN19bt_fit_params.json`, `density_profile_litTDEs.csv`, `bh_masses_host.txt`: Literature values used for comparison
