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
#include <fstream>
#include <vector>
#include "network.h"
using namespace std;

// Default-Konstruktor (Netzwerk ohne Knoten)
Netz::Netz() {
	nNodes = 0;
	nEdges = 0;
}

// Konstruktur mit Angabe von Anzahl der Knoten
Netz::Netz(int n) {
	nNodes = n;
	nEdges = 0;
}

// Copy-Konstruktor
Netz::Netz(const Netz &N) : adList(N.adList) {
	nNodes = N.nNodes;
	nEdges = N.nEdges;
}

// Konstruktor mit Datei
Netz::Netz(string name) {
	const char * fname = name.c_str();
    std::ifstream datei;
	datei.open(fname,ios::in);
	int nE, hn1, hn2;

	datei >> nNodes;
	nEdges = 0;
	datei >> nE;
	// wichtig: immer direkt initialisieren

	for (int ii = 0; ii<nE*2; ii+=2) {
		datei >> hn1;
		datei >> hn2;
		addEdge(hn1, hn2);
	}
	datei.close();
}

// Konstruktur mit Adjazenzlist in Vektor
Netz::Netz(vector<int> connections) : adList(connections) {
	int g_node = 0;
	nEdges = (connections.size())/2;
	if (connections.size()==0) {
		nNodes = 0;
	} else {
		for (size_t ii=0; ii<connections.size(); ii++) {
			g_node = connections[ii]>g_node ? connections[ii] : g_node;
		}
		nNodes = g_node+1;
	}

}


// Rueckgabe Anzahl der Knoten
int Netz::getNumberNodes() {
	return nNodes;
}

// Rueckgabe Anzahl der Kanten
int Netz::getNumberEdges() {
	return nEdges;
}

// Einen Knoten hinzufuegen
void Netz::addNode() {
	nNodes++;
	sayNumberNodes();
}

// Mehrere (n) Knoten hinzufuegen
void Netz::addNodes(int n) {
	nNodes += n;
	sayNumberNodes();
}

// 1. Variante: Gleich mehrere Kanten einfuegen
void Netz::addEdges(vector<int> connections) {
	if (connections.size()%2 != 0) {
		cout << "Vector need's to have an even number of entries." << endl;
		return;
	}
	for (size_t ii=0; ii<connections.size(); ii+=2) {
		addEdge(connections[ii],connections[ii+1]);
	}
}

// 2. Variante (besser die erste verwenden)
void Netz::addEdges(int * connections, int n) {

	for (int ii=0; ii < 2*n; ii+=2) {
		addEdge(*(connections+ii),*(connections+ii+1));
	}

}

// Eine Kante hinzufuegen
void Netz::addEdge(int a, int b) {
	int h;

	if (a==b) {
		cout << "No, same node." << endl;
		return;
	}

	h = b;
	b = b>a ? b : a;
	a = a<=h ? a : h;
	if (b>=nNodes) {
		cout << "Not enough nodes in the network" << endl;
		sayNumberNodes();
		return;
	}

	if (!existsEdge(a,b)) {
	   	adList.push_back(a);
	   	adList.push_back(b);
	   	nEdges++;
	} else {
		// cout << "Connection already exists." << endl;
	}
}

// Abfrage ob Kante schon existiert
bool Netz::existsEdge(int a, int b) {
	int h;
	if (a==b) {
		cout << "Actually, that's the same node." << endl;
		return false;
	} else if (adList.size()==0) {
		return false;
	}
	h = b;
	b = b>a ? b : a;
	a = a<h ? a : h;

	for (int ii=0; ii<2*nEdges; ii+=2) {
		if (adList[ii]==a && adList[ii+1]==b) {
			return true;
		}
	}
	return false;

}

// Ausgabe der Adjazenzliste auf Konsole
void Netz::showAdList() {
	for (int ii=0; ii<2*nEdges; ii+=2) {
		cout << adList[ii] << " - " << adList[ii+1] << endl;
	}
}

// Rueckgabe der Adjazenzliste
vector<int> Netz::gimmeThatAdList() {
	return adList;
}

// Ausgabe der Anzahl der Knoten
void Netz::sayNumberNodes() {
	cout << "Network now consists of " << nNodes << " nodes." << endl;
}

// Ausgabe der Anzahl der Kanten
void Netz::sayNumberEdges() {
	cout << "Network now has " << nEdges << " edges." << endl;
}

int Netz::getAdListEntry(int ind) {
	return adList[ind];
}

void Netz::setAdListEntry(int ind, int eintrag) {
	adList[ind] = eintrag;
}


