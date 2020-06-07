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



#include "myDifferentials.h"

#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>
#include <math.h>



int netVier (double t, const double y[], double f[], void *params) {
	(void) (t);
    paraBullo pa = *(paraBullo *) params;

    int ind_from, ind_to;
    double conflow;

    for (int ii=0; ii<pa.nNodes; ii++) {
    	f[ii] = y[ii+pa.nNodes];
    	f[ii+pa.nNodes] = -pa.alpha*y[ii+pa.nNodes] + pa.P[ii];
    }

    for (int jj = 0; jj<pa.nEdges; jj++) {
    	ind_from = pa.adList[jj*2];
    	ind_to = pa.adList[jj*2+1];

    	conflow = pa.K * sin( y[ind_to]-y[ind_from] );
        f[ind_from+pa.nNodes] += conflow;
        f[ind_to+pa.nNodes]   -= conflow;
    }

	return GSL_SUCCESS;
}


int netVierBack(double t, const double y[], double f[], void *params) {
	(void) (t);
    paraBullo pa = *(paraBullo *) params;

    int ind_from, ind_to;
    double conflow;

    f[4*pa.nNodes-1] = 0.0;
    for (int ii=0; ii<pa.nNodes; ii++) {
    	f[ii] = -y[ii+pa.nNodes];
    	f[ii+pa.nNodes] = pa.alpha*y[ii+pa.nNodes] - pa.P[ii];
    	f[ii+2*pa.nNodes] = 0.0;
    	f[ii+3*pa.nNodes] = y[ii+2*pa.nNodes] - pa.alpha*y[ii+3*pa.nNodes] - (2.0/pa.T * y[ii+pa.nNodes]);

    }


    for (int jj=0; jj<pa.nEdges; jj++) {
    	ind_from = pa.adList[jj*2];
    	ind_to = pa.adList[jj*2+1];

    	conflow = pa.K * sin( y[ind_to]-y[ind_from] );
        f[ind_from+pa.nNodes] -= conflow;
        f[ind_to+pa.nNodes]   += conflow;

        conflow = (y[ind_to+3*pa.nNodes]-y[ind_from+3*pa.nNodes])*pa.K*cos( y[ind_to]-y[ind_from] );
        f[ind_from+2*pa.nNodes] += conflow;
        f[ind_to+2*pa.nNodes]   -= conflow;

    }

	return GSL_SUCCESS;
}



