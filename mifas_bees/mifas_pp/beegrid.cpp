// ...
//
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



#include <fstream>
#include "network.h"
#include "beegrid.h"

using namespace std;

// Konstruktor
Beegrid::Beegrid() : Netz() {
	nBees = 0;
	nFlowers = 0;
	alpha = -0.1;
	beta_ii = 1.0;
	beta_ij = 1.0;
	delta = 0.1;
	h = 0.1;
}

// Konstruktor mit Datei (nur Adjazenz-Liste und Anzahl Knoten)
Beegrid::Beegrid(string name) : Netz(name) {
	nBees = 0;
	nFlowers = 0;
	alpha = -0.1;
	beta_ii = 1.0;
	beta_ij = 1.0;
	delta = 0.1;
	h = 0.1;

	for (int ii=0; ii<getNumberNodes(); ii++) {
		gamma.push_back(0.0);
	}
}

// Konstruktor mit zwei Dateien (Topologie und Parameter)
Beegrid::Beegrid(string name1, string name2) : Netz(name1) {
	const char * fname2 = name2.c_str();
	std::ifstream datei;
	datei.open(fname2,ios::in);
	int nN, test_nb;
	double ph, gamma_null;

	nN = getNumberNodes();
	datei >> alpha;
	datei >> beta_ii;
	datei >> beta_ij;
	datei >> delta;
	datei >> h;
	datei >> gamma_null;

	nBees = 0;
	for (int ii = 0; ii<2*getNumberEdges(); ii+=2) {
		test_nb = getAdListEntry(ii);		
		if (test_nb>nBees) {
			nBees = test_nb;
		}
	}
	nBees += 1;
	nFlowers = nN-nBees;
	

	for (int ii = 0; ii<nN; ii++) {
		gamma.push_back(gamma_null/getDegree(ii));
	}
	datei.close();
}

// Konstruktor
Beegrid::Beegrid(int nB, int nF) : Netz(nB) {
	nBees = nB; 
	nFlowers = nF;
	alpha = 0.1;
	beta_ii = 1.0;
	beta_ij = 1.0;
	delta = 0.1;
	h = 0.1;
	for (int ii=0; ii<(nB+nF); ii++) {
		gamma.push_back(0.0);
	}
}

// Konstruktor
Beegrid::Beegrid(const Netz &N, int nB, int nF) : Netz(N) {
	nBees = nB;
	nFlowers = nF;
	alpha = 0.1;
	beta_ii = 1.0;
	beta_ij = 1.0;
	delta = 0.1;
	h = 0.1;
	for (int ii=0; ii<getNumberNodes(); ii++) {
		gamma.push_back(0.0);
	}
}

// Konstruktor
Beegrid::Beegrid(const Netz &N, vector<double> p, int nB, int nF) : Netz(N) {
	nBees = nB;
	nFlowers = nF;
	alpha = 0.1;
	beta_ii = 1.0;
	beta_ij = 1.0;
	delta = 0.1;
	h = 0.1;
	gamma = p;
}


// Add single node
void Beegrid::addNode(double p = 0.0) {
	Netz::addNode();
	nBees++;
	gamma.push_back(p);
}

// Show info on console
void Beegrid::showNetInfo() {
	cout << "-------------- INFO ---------------" << endl;
	cout << "Beegrid with " << getNumberNodes() << " nodes and " <<
	getNumberEdges() << " edges:" << endl;
	cout << "alpha = " << alpha << endl;
	cout << "beta_ii = " << beta_ii << endl;
	cout << "beta_ij = " << beta_ij << endl;
	cout << "delta = " << delta << endl;
	cout << "h = " << h << endl;
	cout << "Node info: " << endl;
	for (int ii=0; ii<getNumberNodes(); ii++) {
		cout << "N" << ii << " : gamma = " << gamma[ii] << endl;
	}
	cout << "Adjacency list: " << endl;
	showAdList();
	cout << "-----------------------------------" << endl;
}

// Set Parameter alpha
void Beegrid::setAlpha(double a) {
	alpha = a;
}


// Set Parameter from File (if topology is already initialized)
void Beegrid::setParaFromFile(string datei_name) {
	const char * fname = datei_name.c_str();
	std::ifstream datei;
	datei.open(fname,ios::in);
	int nN;
	double ph;

	datei >> nN;
	datei >> nBees;
	datei >> nFlowers;
	datei >> alpha;
	datei >> beta_ii;
	datei >> beta_ij;
	datei >> delta;
	datei >> h;

	for (int ii = 0; ii<nN; ii++) {
		datei >> ph;
		setGamma(ii,ph);
	}
	datei.close();
}



// Set power of single node
void Beegrid::setGamma(int n, double p) {
	if (n>=getNumberNodes()) {
		cout << "Node does not exist." << endl;
	} else {
		gamma[n] = p;
	}
}


// Set powers via vector
void Beegrid::setGammas(vector<double> &p) {
	if (p.size() == getNumberNodes()) {
		for (size_t ii=0; ii<gamma.size(); ii++) {
			gamma[ii] = p[ii];
		}
	} else {
		cout << "Vector does not correpond to number of nodes." << endl;
	}
}


// Ausgabe des Parametersatzes fuer die eigentliche Simulation
paraSimple* Beegrid::createParams() {
	
	struct paraSimple * para;
	para = new paraSimple;
	// BITTE HIER NOCHMAL REINSCHAUEN, DAS GEHT NOCH SCHOENER.
	int * ad;
	// double * pm;
	ad = new int[getNumberEdges()*2];
	// pm = new double[getNumberNodes()];


	para->nb 		= nBees;
	para->nf 		= nFlowers;
	para->nEdges 	= getNumberEdges();
	para->alpha  	= alpha;
	para->beta_ii	= beta_ii;
	para->beta_ij	= beta_ij;
	para->delta 	= delta;
	para->h			= h;

	//for (int ii=0; ii<getNumberNodes(); ii++) {
	//	pm[ii] = gamma[ii];
	//}

	for (int ii=0; ii<getNumberEdges()*2; ii++) {
		ad[ii] = gimmeThatAdList().at(ii);
	}

	// delete [] para->gamma;
	// delete [] para->adList;

	para->gamma = new double[getNumberNodes()];
	para->adList = new int[getNumberEdges()*2];
	para->nap = new double[getNumberNodes()];
	para->f_help = new double[getNumberNodes()];

	for (int ii=0; ii<getNumberNodes(); ii++) {
		para->gamma[ii] = gamma[ii];
		para->nap[ii] = 0.0;
		para->f_help[ii] = 0.0;
	}

	for (int ii=0; ii<getNumberEdges()*2; ii++) {
		para->adList[ii] = ad[ii];
	}

	

	delete [] ad;
	// delete [] pm;

	return para;
}



