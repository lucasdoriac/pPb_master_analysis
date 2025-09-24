# Control Plots

This subfolder contains the macros and scripts used to generate control plots as the first step of the p–Pb master analysis. 
These plots were intended as quality assurance checks, validating reconstruction algorithms and ensuring data consistency before performing the main physics analysis.

## Contents
- To be updated. Example:
- `control_plots.C` -> ROOT macro to generate basic histograms.
- `print_contents.C` -> Print contents and classes of all files inside `.root` datafile.
- `results/` -> local directory where output plots are stored.

## Usage
To run the macro from bash:
```bash
root -l -q control_plots.C
```

### Notes
- Path to `.root` datafiles may need to be updated depending on the environment.

### To-do
- Update `control_plots.C` code.
- Update `Pb_vs_p.C` code.
- Add updated presentation to repo.
- At the very end, when the remote repo will be completed, we need to update the path/to/datafile in all macros, since the datafile is going to be stored in a directory inside the pPb_master_analysis/ folder.
