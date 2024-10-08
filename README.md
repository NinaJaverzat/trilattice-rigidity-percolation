# trilattice-rigidity-percolation

The code implements basic central-force bond rigidity percolation on the triangular lattice, as described in the original works [D. J. Jacobs and M. F. Thorpe Phys. Rev. E 53, 3682 (1996)](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.53.3682) and [D. J. Jacobs, B. Hendrickson Journal of Computational Physics 137, 2 (1997)](https://www.sciencedirect.com/science/article/abs/pii/S0021999197958095).

COMPILATION:

```
make
```
creates an executable named rp.exe

EXECUTION:

```
./rp.exe <LINEAR_LATTICE_SIZE> <NUMBER_OF_TRIALS> <RUN_NUMBER>
```

OUTPUTS:

Snapshots of the network anf of its decomposition in rigid clusters cen be saved for each trial in ./cfgs (the folder must exist) and can be visualised using plotting.py: ```python3 plotting.py``` 

The order parameter (strength of percolating rigid cluster), the probabilities of rigidity percolating, and the susceptibility as function of the filling probability are recorded in ./res (the folder must exist). They can be collapsed using standard scaling analysis as in Scaling_funcs.ipynb (open with ```jupyter notebook```)
