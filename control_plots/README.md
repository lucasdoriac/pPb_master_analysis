# Control Plots

This subfolder contains the ROOT/C++ macros used to generate control plots as the first step of the p–Pb master analysis. 
These plots were intended as quality-assurance checks, validating reconstruction algorithms and ensuring data consistency before performing the main physics analysis.

## Contents
- To be updated. Example:
- `control_plots.C` -> Plots quality-assurance histograms from two `.root` datafiles.
- `plot_hfSumEt_sNN.C` -> Plots the sum of the transverse component of the total collision energy deposited on both HF calorimeters, for both datasets.
- `print_contents.C` -> Print contents and classes of all files inside specific `.root` datafile.
- `results/` -> Still empty. To be updated. Contains plot examples from `control_plots.C` and `plot_hfSumEt_sNN.C`.

## Usage
To run the macros from bash:
```bash
root -l -q control_plots.C
root -l -q plot_hfSumEt_sNN.C
root -l -q print_contents.C
```

### Notes
- Path to the `.root` datafiles may need to be updated depending on the environment.

### To-do
- Update `control_plots.C` code.
- Add updated presentation to repo.
- At the very end, when the remote repo will be completed, we need to update the path/to/datafile in all macros, since the datafile is going to be stored in a directory inside the pPb_master_analysis/ folder.
- I have decided it's not worth it to add the creation of ```results/``` folder when saving the output files. This will be done when the project is finished.