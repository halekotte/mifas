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
#include <vector>
#include <fstream>
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include "network.h"
#include "powergrid.h"
using namespace std;

// Konstruktor
Powergrid::Powergrid() : Netz() {
	alpha = 0.1;
	K = 1;
	T = 10;
}

// Konstruktor mit Datei (nur Adjazenz-Liste und Anzahl Knoten)
Powergrid::Powergrid(string name) : Netz(name) {
	alpha = 0.1;
	K = 1;
	T = 10;

	for (int ii=0; ii<getNumberNodes(); ii++) {
		P.push_back(0);
		omega.push_back(0);
	}
}

// Konstruktor mit zwei Dateien (Topologie und Parameter)
Powergrid::Powergrid(string name1, string name2) : Netz(name1) {
	const char * fname2 = name2.c_str();
	std::ifstream datei;
	datei.open(fname2,ios::in);
	int nN;
	double ph;

	datei >> nN;
	datei >> alpha;
	datei >> K;
	datei >> T;

	for (int ii = 0; ii<nN; ii++) {
		datei >> ph;
		P.push_back(ph);
		omega.push_back(0.0);
	}
	datei.close();
}

// Konstruktor
Powergrid::Powergrid(int n) : Netz(n) {
	alpha = 0.1;
	K = 1;
	T = 10;
	for (int ii=0; ii<n; ii++) {
		P.push_back(0.0);
		omega.push_back(0.0);
	}
}

// Konstruktor
Powergrid::Powergrid(const Netz &N) : Netz(N) {
	alpha = 0.1;
	K = 1;
	T = 10;
	for (int ii=0; ii<getNumberNodes(); ii++) {
		P.push_back(0.0);
		omega.push_back(0.0);
	}
}

// Konstruktor
Powergrid::Powergrid(const Netz &N, vector<double> p) : Netz(N) {
	alpha = 0.1;
	K = 1;
	T = 10;
	P = p;
	for (size_t ii=0; ii<p.size(); ii++) {
		omega.push_back(0.0);
	}
}

// Konstruktor
Powergrid::Powergrid(vector<int> al, vector<double> p,
    double a, double k, double t) : Netz(al), P(p) {

	for (size_t ii=0; ii<p.size(); ii++) {
		omega.push_back(0.0);
	}
	alpha = a;
	K = k;
	T = t;
}

// Konstruktor
Powergrid::Powergrid(vector<int> al, vector<double> p, vector<double> o,
    double a, double k, double t) : Netz(al), P(p), omega(o) {
	alpha = a;
	K = k;
	T = t;
}


// Add single node
void Powergrid::addNode(double p = 0.0, double o = 0.0) {
	Netz::addNode();
	P.push_back(p);
	omega.push_back(o);
}

// Show info on console
void Powergrid::showNetInfo() {
	cout << "-------------- INFO ---------------" << endl;
	cout << "Powergrid with " << getNumberNodes() << " nodes and " <<
	getNumberEdges() << " edges:" << endl;
	cout << "alpha = " << alpha << endl;
	cout << "K     = " << K << endl;
	cout << "T     = " << T << endl;
	cout << "Node info: " << endl;
	for (int ii=0; ii<getNumberNodes(); ii++) {
		cout << "N" << ii << " : P = " << P[ii] <<
		",  omega = " << omega[ii] << endl;
	}
	cout << "Adjacency list: " << endl;
	showAdList();
	cout << "-----------------------------------" << endl;
}

// Set Parameter alpha
void Powergrid::setAlpha(double a) {
	alpha = a;
}

// Set parameter K (capacity):
void Powergrid::setCapacity(double k) {
	K = k;
}

// Set Parameter T (period):
void Powergrid::setPeriod(double t) {
	T = t;
}

// Set Parameter from File (if topology is already initialized)
void Powergrid::setParaFromFile(string datei_name) {
	const char * fname = datei_name.c_str();
	std::ifstream datei;
	datei.open(fname,ios::in);
	int nN;
	double ph;

	datei >> nN;
	datei >> alpha;
	datei >> K;
	datei >> T;

	for (int ii = 0; ii<nN; ii++) {
		datei >> ph;
		setP(ii,ph);
	}
	datei.close();
}

// Set omega of single node
void Powergrid::setOmega(int n, double o){
	if (n>=getNumberNodes()) {
		cout << "Node does not exist." << endl;
	} else {
		omega[n] = o;
	}
}

// Set power of single node
void Powergrid::setP(int n, double p) {
	if (n>=getNumberNodes()) {
		cout << "Node does not exist." << endl;
	} else {
		P[n] = p;
	}
}

// set omegas via vector
void Powergrid::setOmegas(vector<double> &o) {
	if (o.size() == getNumberNodes()) {
		for (size_t ii=0; ii<omega.size(); ii++) {
			omega[ii] = o[ii];
		}
	} else {
		cout << "Vector does not correpond to number of nodes." << endl;
	}
}

// Set powers via vector
void Powergrid::setPs(vector<double> &p) {
	if (p.size() == getNumberNodes()) {
		for (size_t ii=0; ii<P.size(); ii++) {
			P[ii] = p[ii];
		}
	} else {
		cout << "Vector does not correpond to number of nodes." << endl;
	}
}

// Check if consumed and produced power match
bool Powergrid::isBalancedP() {
	double check_P;
	check_P = 0.0;
	for (int ii=0; ii<getNumberNodes(); ii++) {
		check_P += P[ii];
	}
	return (check_P==0.0);
}

// Omega bestimmen (Simulation mit GSL-Paket)
void Powergrid::determineOmega(paraBullo * params) {
	double t, T;
	double * y;
	unsigned long dimVar;
	double K_temp, T_temp;
	bool repeat_procedure = false;
	dimVar = getNumberNodes()*2;
	y = new double[dimVar];
	for (size_t ii=0; ii<dimVar; ii++) {
		y[ii] = 0.0;
	}

	t = 0.0;
	T = 2500.0;

    // Set up system
    gsl_odeiv2_system sys = {netVier, 0, dimVar, params};

    // Choose parameters for integration
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-9, 1e-12, 1e-12);

    // Integrate
    int status = gsl_odeiv2_driver_apply (d, &t, T, y);

    if (status != GSL_SUCCESS) {
        printf ("error, return value=%d\n", status);
    }


	// ----- NEU: Um auch weniger stabile Systeme untersuchen zu koennen.
	for (int jj=0; jj<getNumberNodes(); jj++) {
		if (abs(y[jj+getNumberNodes()])>0.01) {
			repeat_procedure = true;
			break;
		} 
	}

	if (repeat_procedure) {
		gsl_odeiv2_driver_reset(d);
			
		for (int jj=0; jj<dimVar; jj++) {
			y[jj] = 0.0;
		} 

		K_temp = 10.0;
		T_temp = 100.0;

		params->K = K_temp;
		t = 0.0;

		status = gsl_odeiv2_driver_apply (d, &t, T, y);
		gsl_odeiv2_driver_reset(d);

		while (K_temp > this->K+0.2) {
			t = 0.0;
			K_temp -= 0.1;
			params->K = K_temp;

			status = gsl_odeiv2_driver_apply (d, &t, T_temp, y);
			gsl_odeiv2_driver_reset(d);
		} 
		
		t = 0.0;
		params->K = this->K; 
		status = gsl_odeiv2_driver_apply (d, &t, T, y);
		gsl_odeiv2_driver_reset(d);	
		
	}
	// ------------------------------------------------------------------


    for (int jj = 0; jj < getNumberNodes(); jj++) {
    	omega[jj] = y[jj];
    }

    gsl_odeiv2_driver_free (d);
    delete[] y;
    y = 0;

}


// Ausgabe des Parametersatzes fuer die eigentliche Simulation
paraBullo * Powergrid::createParams(paraBullo * para, bool detOmega = true) {
	int * ad;
	double * pm, * om;
	ad = new int[getNumberEdges()*2];
	pm = new double[getNumberNodes()];
	om = new double[getNumberNodes()];

	para->nNodes = getNumberNodes();
	para->nEdges = getNumberEdges();
	para->alpha  = alpha;
	para->K      = K;
	para->T      = T;

	for (int ii=0; ii<getNumberNodes(); ii++) {
		pm[ii] = P[ii];
		om[ii] = omega[ii];
	}

	for (int ii=0; ii<getNumberEdges()*2; ii++) {
		ad[ii] = gimmeThatAdList().at(ii);
	}

	delete [] para->P;
	delete [] para->omega;
	delete [] para->adList;

	para->P = new double[getNumberNodes()];
	para->omega = new double[getNumberNodes()];
	para->adList = new int[getNumberEdges()*2];

	for (int ii=0; ii<getNumberNodes(); ii++) {
		para->P[ii] = P[ii];
		para->omega[ii] = om[ii];
	}

	for (int ii=0; ii<getNumberEdges()*2; ii++) {
		para->adList[ii] = ad[ii];
	}


	if (detOmega) {
		determineOmega(para);
	}

	for (int ii=0; ii<getNumberNodes(); ii++) {
		para->omega[ii] = omega[ii];
	}

	delete [] ad;
	delete [] om;
	delete [] pm;

	return para;
}


void Powergrid::changePositionOfNode(int a, int b)  {
    double P_help;
    int index;
    P_help = P[a];
    P[a] = P[b];
    P[b] = P_help;
    for (size_t ii=0; ii<getNumberNodes(); ii++) {
    	omega[ii] = 0.0;
    }

    for(size_t ii=0; ii<getNumberEdges()*2; ii++) {
    	index = getAdListEntry(ii);
    	if (index == a) {
    		setAdListEntry(ii, b);
    	} else if (index == b) {
    		setAdListEntry(ii,a);
    	}
    }

}

double Powergrid::getP(int ind) {
	return P[ind];
}



