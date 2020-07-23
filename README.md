# UEP

### Welcome to UEP!
<p align="justify">
UEP is a tool for predicting the impact of mutations in a protein-protein complex. For a given PDB structure, it predicts automatically the effects of all possible mutations in the highly-packed interface positions. UEP predictions are competitive with the state-of-the art methods in the field with the advantage of being open-source code and extremelly fast in comparison, since it will take less than a second!
</p>

### UEP dependencies

- [argparse](https://docs.python.org/3/library/argparse.html) - user input
- [prody](https://github.com/prody/ProDy) - three-dimensional searches
- [itertools](https://docs.python.org/3/library/itertools.html) - combinatorial calculations
- [os](https://docs.python.org/3/library/os.html) - pathing files
- [numpy](https://numpy.org/) - calculations
- [compress_pickle](https://pypi.org/project/compress-pickle/) - reading UEP contact matrix
- [pandas](https://pandas.pydata.org/) - exporting results as csv file

### Running UEP

1. Clone UEP repository in your computer.
2. Define a PDB file and an interface to work on following this scheme.
```
python3 UEP.py --pdb=PDB.pdb --interface=A,BC
```
3. Results will be displayed and saved (csv file) in the same folder than your PDB path.

### Understanding UEP results

1. Open the results file (csv file).
2. First column represents analyzed positions on the highly packed interface.
3. Other columns represent mutations into the different residues.
4. Numerical values represent the predicted ΔΔG.
5. NaN values represent mutations that could not be scored because: i) mutation is the same residue than the wild type, ii) mutation has less than 2 predicted contacts with the other chains.
6. Negative ΔΔG values are mutations predicted to improve the binding affinity of the PPI.
7. Positive ΔΔG values are mutations predicted to decrease the binding affinity of the PPI.

### What makes UEP different from the state-of-the art?

<p align="justify">
Current state-of-the art methods for predicting the impact of mutations in a protein-protein complex rely on the description of physical energies, statistical potentials, conservation, shape complementarity, and more recently, machine learning-based approaches.

UEP moves appart from the state-of-the art and it is based on the interactions observed in the interactome data (https://interactome3d.irbbarcelona.org/). It follows a three-body contact scheme of the highly-packed positions, where one residue of one protein must be in contact with at least two residues of the other protein. We have observed that such highly-packed positions exert larger differences in the experimental ΔΔG, and therefore they i) are easier to be predicted, and ii) they are more interesting for protein-protein design campaings.
</p>

<p align="center">
<img src="images/uep_scheme.png" width="400">
</p>

<p align="justify">
Once you run UEP, it will find the highly-packed residues of your PDB, and it will examine the contacts of your protein-protein interface. Then, it will predict a ΔΔG based on the wild type and the mutation counts observed in the interactome data, without the need of generating mutation files. This feature makes UEP really fast!
</p>


