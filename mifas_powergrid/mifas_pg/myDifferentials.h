// Copyright (C) 2020  Lukas Halekotte <lukas.halekotte@uol.de>
//
// This file is part of the program 'mifas_pg': you can redistribute it 
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
struct paraBullo {
	int nNodes;
	int nEdges;
	double alpha;
	double K;
	double T;
	double * P;
	double * omega;
	int * adList;
};

/* ===== Nummer 1:
 * Differentialgleichung fuer beliebiges Netzwerk
 * mit letztem Knoten als Referenz
 */
int netVier (double t, const double y[], double f[], void *params);
int netVierBack (double t, const double y[], double f[], void *params);




#endif /* MYDIFFERENTIALS_H_ */
