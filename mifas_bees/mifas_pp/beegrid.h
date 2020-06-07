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



#ifndef BEEGRID_H_
#define BEEGRID_H_

#include <iostream>
#include <vector>
#include "network.h"
#include "myDifferentials.h"

using namespace std;

class Beegrid : public Netz {
	int nBees; 
	int nFlowers;
	double alpha;
	double beta_ii;
	double beta_ij;
	double delta;
	double h;
	vector<double> gamma;

public:
	// Default-Konstruktor:
	Beegrid();

	// Konstruktor mit Datei fuer Topologie, weieter Parameter auf Def.:
	Beegrid(string);

	// Konstruktor mit Datei fuer Topologie und Datei fuer weitere Parameter:
	Beegrid(string, string);

	// Weitere Konstruktoren (werden erstmal nicht gebraucht):
	Beegrid(int, int);
	Beegrid(const Netz&, int, int);
	Beegrid(const Netz &, vector<double>, int, int);


	// add singel node:
	void addNode(double);
    // add single node mit Angabe von P(i) und omega(i):
	void addNode(double, double);

	// Ausgabe Informationen uerb das Netz (Konsole):
	void showNetInfo();

	// Set single parameter:
	void setAlpha(double);    // alpha
	
	// Set node parameter:
	void setParaFromFile(string);
	void setGamma(int, double);
	void setGammas(vector<double> &);

    
    // Ausgabe Parameterset:
    paraSimple * createParams();

};

#endif /* BEEGRID_H_ */


