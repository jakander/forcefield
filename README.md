# forcefield
Script pulls copper-water coordinates from a QM calculation using B3LYP/cc-PVTZ computing a potential energy scan of one water
incrementing from long distance (~4 angstroms) to short distance (~1.5 angstroms). In order to use the script on multiple scans,
I attempt to leave the code as general as possible. Will likely only work for Gaussian log files. Script then computes two 
body interactions (Lennard Jones and Coulombic) and plots this potential as well as the potentail energy scan computed in the 
QM calculation. 
