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



#ifndef POWERGRID_H_
#define POWERGRID_H_

#include <iostream>
#include <vector>
#include "network.h"
#include "myDifferentials.h"
using namespace std;

class Powergrid : public Netz {
	double alpha;
	double K;
	double T;
	vector<double> P;
	vector<double> omega;
public:
	// Default-Konstruktor:
	Powergrid();

	// Konstruktor mit Datei fuer Topologie, weieter Parameter auf Def.:
	Powergrid(string);

	// Konstruktor mit Datei fuer Topologie und Datei fuer weitere Parameter:
	Powergrid(string, string);

	// Weitere Konstruktoren (werden erstmal nicht gebraucht):
	Powergrid(int);
	Powergrid(const Netz&);
	Powergrid(const Netz &, vector<double>);
	Powergrid(vector<int>, vector<double>, double, double, double);
	Powergrid(vector<int>, vector<double>, vector<double>, double, double, double);


	// add singel node:
	void addNode();
    // add single node mit Angabe von P(i) und omega(i):
	void addNode(double, double);

	// Ausgabe Informationen uerb das Netz (Konsole):
	void showNetInfo();

	// Set single parameter:
	void setAlpha(double);    // alpha
	void setCapacity(double); // K
	void setPeriod(double);   // T

	// Set node parameter:
	void setParaFromFile(string);
	void setOmega(int, double);
	void setP(int, double);
	void setPs(vector<double> &);
	void setOmegas(vector<double> &);

	// Check if P-values sum up to 0:
    bool isBalancedP();

    // Bestimmen der Werte von omega fuer bestehendes Netz:
    void determineOmega(paraBullo *);

    // Ausgabe Parameterset:
    paraBullo * createParams(paraBullo *, bool);

    // Um Referenzknoten zu aendern
    void changePositionOfNode(int a, int b);

    double getP(int);
};

#endif /* POWERGRID_H_ */
