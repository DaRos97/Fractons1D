# Fractons in 1D

Description of the code used for the simulations in the master's thesis [Fractons: analytical description and simulations of their thermal behavior](https://amslaurea.unibo.it/23919/).

[![Python 3.9.6](https://img.shields.io/badge/python-3.9.6-blue.svg)](https://www.python.org/downloads/release/python-396/)

## General description

In this repository are constructed spin-1 **dipole** conserving Hamiltonians using the library [QuSpin](https://weinbe58.github.io/QuSpin/) .
Are considered two Hamiltonians: with *three* and *four* site interactions.  
Various observables are considered to analyze their time-evolution and thermal behavior:

- Autocorrelator <img src="https://render.githubusercontent.com/render/math?math=\langle S_z(t)S_z(0)\rangle">
- Inverse Participation Ratio
- Expectation values of <img src="https://render.githubusercontent.com/render/math?math=(S^z)^2"> on all eigenstates
- Level spacing statistics
- Entanglement entropy of eigenstates (still missing)

## Installation

Packages needed to run the scripts are included in the 'requirements.txt' file. 
I suggest to use an [Anaconda environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands) to make everything more smooth. 
Then just type from the terminal:  

```
conda install --file requirements.txt
```

Type `conda list` to check that the various packages are installed.  
Then install [QuSpin](https://weinbe58.github.io/QuSpin/) by typing:  

```
conda install -c weinbe58 quspin
```

## How to use the scripts

In order to use this code to make your own simulations of dipole conserving chains there are a few steps to follow:

1. Install the requested packages

2. Clone the code:  
   
   ```
   git clone https://github.com/DaRos97/Dipole_conserving_chain.git
   ```

3. Now you can start running the code in the following **order**: 
   
   1. Construct the states of the sectors you want to analyze with the `States.py` file, where you need to add the total charge and dipole moment. 
   2. Construct the basis and compute the eigenvalues and eigenfunctions with `Basis.py` (carefull that the computation time and memory needed to store the eigenvectors grows exponentially). 
   3. Now you can compute the `IPR`, the `autocorrelator`, the expectation value of `Sz^2` and the `level spacing statistics` running the corresponding scripts. Plots are created automatically to show the relevant differences between the two considered models.

## Notes on single scripts

- `States.py` constructs all the possible basis states with a determined charge and dipole moment. These are the eigenvalues of the operators

  <img src="https://render.githubusercontent.com/render/math?math=\Large Q=\sum_iS_i^z">
  <img src="https://render.githubusercontent.com/render/math?math=\Large P_{i_0}=\sum_i(i-i_0)S_i^z">
  
  where <img src="https://render.githubusercontent.com/render/math?math=i_0"> is a reference point here taken to be the leftmost of the chain. As an example, the legth-7 chain in the state
  
  <img src="https://render.githubusercontent.com/render/math?math=\Large %2B00%2B0-0">
  
  will have charge <img src="https://render.githubusercontent.com/render/math?math=q=1"> and dipole <img src="https://render.githubusercontent.com/render/math?math=p=-1">.
  
  States are saved in the `data/States` folder.

- `Basis.py` creates the basis using the QuSpin library and the two Hamiltonians

  <img src="https://render.githubusercontent.com/render/math?math=\Large H_3=\sum_i\left(S_i^-(S_{i%2B1}^%2B)^2S_{i%2B2}^- %2B h.c.\right)">  
  <img src="https://render.githubusercontent.com/render/math?math=\Large H_{34}=H_3%2B\sum_i\left(S_i^-S_{i%2B1}^%2B S_{i%2B2}^%2BS_{i%2B3}^- %2B h.c\right)">
  
  Then it diagonalizes the Hamiltonians in order to recover the eigenvalues and the eigenvectors. Since the computation of the eigenvectors is pretty long, an option is given to compute all the spectrum (option [1]) or only the eigenvalues (option [0]).

- `Autocorrelator.py` evaluates the autocorrelation function

  <img src="https://render.githubusercontent.com/render/math?math=\Large C_0(t)=\langle S_0^z(t)S_0^z(0)\rangle">
  
  where <img src="https://render.githubusercontent.com/render/math?math=S^z"> is the usual spin-1 Pauli operator and the pedix 0 refers to the chain's central site it acts upon.
  
  In order to run this script are only necessary the states created by `States.py`. Will be displayed a graph showing the time evolution of the autocorrelator of the 2 different Hamiltonians.

- `IPR.py` calculates the Inverse Participation Ratio of the systems' eigenstates. This is defined as

  <img src="https://render.githubusercontent.com/render/math?math=\Large IPR(\epsilon_k)=\sum_lp_l(\epsilon_k)^2">
  
  where each eigenstate is expressed in terms of a basis <img src="https://render.githubusercontent.com/render/math?math=\left\{\mid l\rangle\right\}"> with some coefficients <img src="https://render.githubusercontent.com/render/math?math=c_{l,i}"> as
  
  <img src="https://render.githubusercontent.com/render/math?math=\Large\mid\epsilon_i\rangle=\sum_lc_{l,i}\mid l\rangle">
  
  In order for this script to run are needed all the eigenvalues and eigenstates of the system. Thus it must be first run `States.py` and then `Basis.py` chosing the option [1] (compute the entire spectrum). The value of the IPR is plotted in function of the energy of the eigenstates.

- `Sz2.py`evaluates the expectation value of the operator

  <img src="https://render.githubusercontent.com/render/math?math=\Large(S_0^z)^2">
  
  on all eigenstates of the systems.
  
  In order for this script to run are needed all the eigenvalues and eigenstates of the system. Thus it must be first run `States.py` and then `Basis.py` chosing the option [1] (compute the entire spectrum). The expectation values are plotted in terms of the energy eigenstates.

- `Level_statistic.py`evaluates the energy level's spacings distribution.
  
  This only needs the eigenvalues of the system so one can chose the option [0] (only eigenvalues) when running `Basis.py`.

## Additional notes

The basis dimension grows exponentially in system size as well as the computation time.  The largest subsector of the Hilbert space corresponds to the `(q=charge,p=dipole moment)=(0,0)` quantum numbers.

To see some simulations and comments on the results read [here](https://amslaurea.unibo.it/23919/).
