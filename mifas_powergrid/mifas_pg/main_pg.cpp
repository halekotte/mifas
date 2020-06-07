// Copyright (C) 2020  Lukas Halekotte <lukas.halekotte@uol.de>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.



#include <iostream>
// #include <string>
// #include <cstdlib>
// #include <sstream>
using namespace std;

#include "myDifferentials.h"
#include "mifas_new.h"
#include "mifas_smooth.h"
#include "network.h"
#include "powergrid.h"


// #include <ctime>


int main(int argc, char *argv[]) {

	int run_id = 0;

	if (argc >= 2) {
		int eingabe = stoi(argv[1]);
		run_id = eingabe;
	}


	// --- declare parameterset:
	struct paraBullo * params;
	params = new paraBullo;

	// --- define functions for the forward- and backward-integration:
	int (*func)(double, const double *, double *, void *) = netVier;
	int (*func_back)(double, const double *, double *, void *) = netVierBack;


	// --- create power grid:
	// Powergrid definieren und Ausgabe des zugehoerigen Parametersatzes:
	Powergrid meinNetz("netz_GB/GB_topo.txt", "netz_GB/GB_para.txt");
	meinNetz.setPeriod(10.0);
	params = meinNetz.createParams(params, true);
	// meinNetz.showNetInfo();


	// --- number of state variables (phases and frequencies):
	long unsigned int numVar = params->nNodes*2;


	// --- if the desired state is not stable, stop the simulation:
	for (int jj=0; jj<numVar/2; jj++) {
		if (params->omega[jj] > 10.0) {
			std::cout << "this network is not locally stable" << std::endl;
			exit(23);
		}
	}


	// --- determine which state variables should be taken into account
	// for the perturbation and the objective function
	bool * perturb_null, * var_norm;
	perturb_null = new bool[numVar];
	var_norm = new bool[numVar];
	for (size_t ii=0; ii<numVar; ii++) {
		if (ii<numVar/2) {
			perturb_null[ii] = false;
			var_norm[ii] = false;
		} else {
			perturb_null[ii] = false;
			var_norm[ii] = true;
		}
	}
   	

	// --- which state variables should be taken into account in the perturbation
	// (here we consider all frequencies or phase velocities)
	for (size_t ii=0; ii<numVar/2; ii++) {
		perturb_null[numVar/2+ii] = true;
	}


	// --- number of runs (local MiFaS):
	unsigned long complete_runs = 1;       

	// --- storage location for results:
	string name_results = "results/mifas_GB_" + std::to_string(run_id) + ".dat";


	// determine_mifas_smooth(params, func, func_back, numVar, perturb_null, var_norm, name_results, complete_runs);
	determine_mifas_new(params, func, func_back, numVar, perturb_null, var_norm, name_results, complete_runs);


	delete [] perturb_null;
	delete [] var_norm;

	

    cout << "\nEnde." << endl;


	return 0;
}

