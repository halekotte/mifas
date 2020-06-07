
Minimal Fatal Shocks in Plant-Pollinator Networks
==========

This code calculates the Minimal Fatal Shocks (MiFaS) for the 
plant-pollinator networks as described in *Halekotte & Feudel (2020)*,
"Minimal Fatal Shocks in Multistable Complex Networks".

The main component of the code is a search algorithm which is based on 
an optimization using Lagrange multipliers. It is written in C++. 
Numerical integration is carried out using the GNU Scientific Library (GSL).

This program is licensed under the GPL v3 or later (see file *COPYING*).


Files and Directories:
--------

The topologies of the plant-pollinator networks which are required to generate the
corresponding dynamical sytems can be found in the directory **pol_topos**. The 
topologies come from the *Web of Life dataset (www.web-of-life.es)*. The 
original sources of the topologies can be found in in the file **original_sources.txt**.
The parametrization of the plant-pollinator networks is contained in the file 
**bees_para.txt**.
Results of the calculations are saved to the directory **results**.

The C++ code can be found in the directory **mifas_pp** and includes the following 
files:

1. **main_pp.cpp**: The main file whose function is to call routines to generate a
plant-pollinator network and the search algorithm which calculates the MiFaS for this
system.

2. **mifas.cpp**: This file includes the search algorithm which is used to determine
the MiFaS for a given plant pollinator system. Contains functions which are needed 
for the random initial search and the non-random search algorithm.

3. **myDifferentials.cpp**: This file holds the differential equations for the forward-
and backward-integration for the plant-pollinators networks. Both differential equations 
are required for the optimization process.

4. **beegrid.cpp**: This file holds the class Beegrid which holds information on a 
specific plant-pollinator network.

5. **network.cpp**: This file holds the class Netz which is the base for the class 
Beegrid.

In accordance with the corresponding cpp-files, the directory 'mifas_pp' also contains 
the header files **mifas.h**, **myDifferentials.h**, **beegrid.h** and **network.h**.


How to run the program.
-----------------------

Run make and the execute the executable file 'run_mifas_pp' with one
parameter giving the index of the network topology should be use 
(in this, the code is run for topology 1).

```sh
$ make
$ ./run_mifas_pp 1 
```



