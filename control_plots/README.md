# Control Plots

This subfolder contains the macros and scripts used to generate control plots as the first step of the p–Pb master analysis. 
These plots were intended as quality assurance checks, validating reconstruction algorithms and ensuring data consistency before performing the main physics analysis.

## Contents
- To be updated. Example:
- `macro_control_plots.C` — ROOT macro to generate basic histograms.
- `plot_config.json` — configuration file for histogram styling.
- `results/` — directory where output plots are stored (not tracked in Git).

## Usage
To run the macro inside ROOT:
```bash
root -l -q macro_control_plots.C
```

### Notes
- Data input paths may need to be updated depending on the environment.

### To-do
- Add efficiency plots.
- Automate legend and axis labeling.