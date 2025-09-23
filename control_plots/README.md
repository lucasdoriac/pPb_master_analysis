# Control Plots

This subfolder contains the macros and scripts used to generate control plots as the first step of the p–Pb master analysis. 
These plots were intended as quality assurance checks, validating reconstruction algorithms and ensuring data consistency before performing the main physics analysis.

## Contents
- To be updated. Example:
- `macro_control_plots.C` -> ROOT macro to generate basic histograms.
- `results/` -> local directory where output plots are stored.

## Usage
To run the macro from bash:
```bash
root -l -q control_plots.C
```

### Notes
- Data input paths may need to be updated depending on the environment.

### To-do
- Update main code.
- Add efficiency plots.
- Automate legend and axis labeling.
