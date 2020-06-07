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
#include "myDifferentials.h"
#include "network.h"
#include "beegrid.h"
#include "mifas.h"


// =============== MAIN ================
int main(int argc, char *argv[]) {

	std::cout << "Here we go." << std::endl;

	
	if (argc == 2) {

		
		// --- generate plant-pollinator network
		string name_topo = "pol_topos/topo_" + string(argv[1]) + ".txt";
		Beegrid myNet(name_topo, "bees_para.txt");
		
		long unsigned int dimVar = myNet.getNumberNodes();
		std::cout << "Variablen " << dimVar << std::endl;

		// --- choose which state variables to take into account (in this case, all)
		bool * perturb_null;
		perturb_null = new bool[dimVar];
		for (int ii=0; ii<dimVar; ii++) {
			perturb_null[ii] = true;
		}	 

		// --- storage location for results
		string name_results = "results/MiFaS_PP" + string(argv[1]) + ".txt";


		// --- number of runs or number of local MiFaS to calculate ---
		unsigned long complete_runs = 5;

		// --- Choose initial intgeration time for the optimization ---
		double T_ini = 5.0;

		// ===  Start claculation of the MiFaS  ===
		determine_mifas(&myNet, T_ini, perturb_null, name_results, complete_runs); 


	} else {

		std::cout << "Which topology?" << std::endl;

	}

	std::cout << "FINITE" << std::endl;
}



