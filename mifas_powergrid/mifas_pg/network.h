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



#ifndef NETWORK_H_
#define NETWORK_H_


#include <iostream>
#include <vector>
using namespace std;

// Klasse fuer ungerichtete Netzwerke
class Netz {
	int nNodes;
	int nEdges;
	vector<int> adList;

public:
	Netz();
	Netz(int);
	Netz(vector<int>);
	Netz(string);
	Netz(const Netz &);

	int getNumberNodes();
	int getNumberEdges();
	void addNode();
	void addNodes(int);
	void addEdge(int, int);
	void addEdges(vector<int>);
	void addEdges(int*, int);
	bool existsEdge(int, int);
	void sayNumberNodes();
	void sayNumberEdges();
	void showAdList();
	int getAdListEntry(int);
	void setAdListEntry(int, int);
	vector<int> gimmeThatAdList();

	// fuer die Initialisierung aus einer Datei waere vielleicht eine
	// friend-Funktion ganz gut geeignet
};


#endif /* NETWORK_H_ */
