// Copyright (C) 2020  Lukas Halekotte <lukas.halekotte@uol.de>
//
// This file is part of the program 'mifas_pp': you can redistribute it 
// and/or modify it under the terms of the GNU General Public License as 
// published by the Free Software Foundation, either version 3 of the 
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.



#ifndef MYDIFFERENTIALS_H_
#define MYDIFFERENTIALS_H_

/* ===== Parmeterset 1:
 * Parameter fuer das Netzwerk
 */
struct paraSimple {
	int nb; 		// number of bees 
	int nf;			// number of flowers
	int nEdges; 	// number of connections
	double alpha;	// intrinsic growth/death rate
	double beta_ii;	// interspecific competition
	double beta_ij;	// intraspecific competition 
	double delta; 	// Allee effect	
	double h;
	double * gamma; // mutualistic strength
	int * adList;
	double * nap;
	double * f_help;
};

struct paraBack {
	int nb; 		// number of bees 
	int nf;			// number of flowers
	int nEdges; 	// number of connections
	double alpha;	// intrinsic growth/death rate
	double beta_ii;	// interspecific competition
	double beta_ij;	// intraspecific competition 
	double delta; 	// Allee effect	
	double h;
	double * gamma; // mutualistic strength
	int * adList;
	double * nap;
	double * f_help;
};


int netBee (double t, const double y[], double f[], void *params);
int netBee_B (double t, const double y[], double f[], void *params);


#endif /* MYDIFFERENTIALS_H_ */
