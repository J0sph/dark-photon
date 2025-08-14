# Fit of the π⁰ → e⁺e⁻γ channel

## Fit_pizero_peak.py

Fit to the e⁺e⁻γ final state, including the π⁰ and η signal peaks and the background
- Run via: `lb-conda default python Fit_pizero_peak.py`
- `filenumber = -1` : Create a flder for each block containing the ROOT file obtained by merging all the files in `/eos`, or use the existing file in that path.
- `filenumber ≠ -1`+ : Create a folder for each block containing the ROOT file obtained by merging the specified number of files in the current path, or use the existing file in that path.

## blocks.py
Create a merged ROOT file for each block in `/eos`, split into subblocks containing the number of files specified by the variable `block_size`.
- Run via: `lb-conda default python blocks.py -r X` <br>
  Where X takes values a, b, c, d, e, …  <br>
  For example, if `block_size = 100`, the `-r a` flag creates a ROOT file by combining the first 100 files, the `-r b` flag creates a ROOT file by combining the next 100 files, and so on for the other flags.
- Then, all the subblocks can be combined into a final ROOT file using the `hadd` command, in a directory with enough space, for example `/tmp`.

## reduce_root_file.cpp
It can take the final ROOT file as input and reduce it, keeping only the variables needed for the analysis to decrease its size and store it in another directory.
- Compile with: `g++ -o reduce_root_file reduce_root_file.cpp $(root-config --cflags --libs)`
- Run via: ` ./reduce_root_file [input.root] [reduced.root]`

## Plot_pi0_fit_values_per_data_block.py
Generates a summary of the fit results for each of the blocks.
- `lb-conda default python Plot_pi0_fit_values_per_data_block.py`
