
Minimal Fatal Shocks in the Great Britain Power Grid
==========

This code calculates the Minimal Fatal Shock (MiFaS) for the 
Great Britain power grid as described in *Halekotte & Feudel (2020)*,
"Minimal Fatal Shocks in Multistable Complex Networks".

The main component of the code is a search algorithm which is based on 
an optimization using Lagrange multipliers. It is written in C++. 
Numerical integrations are carried out using the GNU Scientific Library (GSL).

This program is licensed under the GPL v3 or later (see file *COPYING*).


Files and Directories:
--------

The topology and parameter set for generating the dynamical system corresponding
to the Great Britain power grid are contained in the directory **netz_GB**. The 
topology is extracted from *Simonsen et al., Phys. review letters 100, 218701 (2008)*.
Results of the calculations are saved to the directory **results**.

The C++ code can be found in the directory **mifas_pg** and includes the following 
files:

1. **main_pg.cpp**: The main file whose function is to call routines to generate a
plant-pollinator network and the search algorithm which calculates the MiFaS for this
system.

2. **mifas_new.cpp**: This file includes the search algorithm which is used to determine
the MiFaS for a given power grid. Contains functions which are needed for the random 
initial search and the non-random search algorithm.

3. **mifas_smooth.cpp**: This file includes an alternative search algorithm for 
calculating the MiFaS for a given power grid. This algorithm is not described in 
*Halekotte & Feudel (2020)*! 

4. **myDifferentials.cpp**: This file holds the differential equations for the forward-
and backward-integration for the power grid. Both differential equations are required 
to run the optimization process.

5. **powergrid.cpp**: This file holds the class Powergrid which holds information on a 
specific power grid.

6. **network.cpp**: This file holds the class Netz which is the base for the class 
Powergrid.

In accordance with the corresponding cpp-files, the directory 'mifas_pp' also contains 
the header files **mifas_new.h**, **mifas_smooth.h**, **myDifferentials.h**, **powergrid.h** 
and **network.h**.


How to run the program.
-----------------------

Run make and the execute the executable file 'run_mifas_pp'.

```sh
$ make
$ ./run_mifas_pg
```



