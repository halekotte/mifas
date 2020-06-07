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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>


#define NO_ALLEE
#ifdef NO_ALLEE
int netBee (double t, const double y[], double f[], void *params) {
	(void) (t);
    paraSimple pa = *(paraSimple *) params;
  
	int ind_a, ind_b;
	double sumBees, sumFlowers, mutualism; 
	sumBees = 0.0;
	sumFlowers = 0.0;
	mutualism = 0.0;

	for (int ii=0; ii<pa.nb; ii++) {
		sumBees += y[ii];
		f[ii] = 0.0;
	}

	for (int ii=pa.nb; ii<(pa.nf+pa.nb); ii++) {
		sumFlowers += y[ii];
		f[ii] = 0.0;
	}

	for (int ii=0; ii<pa.nEdges; ii++) {
		ind_a = pa.adList[2*ii];
		ind_b = pa.adList[2*ii+1];

		f[ind_a] += y[ind_b];
		f[ind_b] += y[ind_a];
	}

	for (int ii=0; ii<pa.nb; ii++) {
		mutualism = pa.gamma[ii]*f[ii]/(1.0+pa.h*pa.gamma[ii]*f[ii]);
		f[ii] = ( pa.alpha  -  pa.beta_ii*y[ii]  -  (pa.beta_ij/(pa.nb-1.0))*(sumBees-y[ii])  + 
			mutualism ) * y[ii];
	}

	for (int ii=pa.nb; ii<(pa.nb+pa.nf); ii++) {
		mutualism = pa.gamma[ii]*f[ii]/(1.0+pa.h*pa.gamma[ii]*f[ii]);
		f[ii] = ( pa.alpha  -  pa.beta_ii*y[ii]  -  (pa.beta_ij/(pa.nf-1.0))*(sumFlowers-y[ii])  + 
			mutualism ) * y[ii];
	}


	return GSL_SUCCESS;
}




int netBee_B (double t, const double y[], double f[], void *params) {
	// zumindest die Formulierung sollte jetzt stimmen
	(void) (t);
    // paraSimple pa = *(paraSimple *) params; 
	paraSimple pa = *(paraSimple *) params;

	int ind_a, ind_b;
	double sumBees, sumBees_extra, sumFlowers, sumFlowers_extra;
	
	sumBees = 0.0;
	sumBees_extra = 0.0;
	sumFlowers = 0.0;
	sumFlowers_extra = 0.0;


	// (1) Schleifen ueber Knoten
	for (int ii=0; ii<pa.nb; ii++) {
		sumBees += pa.nap[ii];
		sumBees_extra += pa.nap[ii]*y[ii];
		f[ii] = 0.0;
		pa.f_help[ii] = 0.0;
	}

	for (int ii=pa.nb; ii<(pa.nf+pa.nb); ii++) {
		sumFlowers += pa.nap[ii];
		sumFlowers_extra += pa.nap[ii]*y[ii];
		f[ii] = 0.0;
		pa.f_help[ii] = 0.0;
	}


	
	// (2) Schleife ueber Adjazenzliste
	for (int ii=0; ii<pa.nEdges; ii++) {
		ind_a = pa.adList[2*ii];
		ind_b = pa.adList[2*ii+1];

		pa.f_help[ind_a] += pa.nap[ind_b];
		pa.f_help[ind_b] += pa.nap[ind_a];
	}



	// (3) Schleifen ueber Knoten
	for (int ii=0; ii<pa.nb; ii++) {
		f[ii] = pa.alpha*y[ii]  -  (2.0*pa.beta_ii*pa.nap[ii] + (pa.beta_ij/(pa.nb-1.0))*(sumBees-pa.nap[ii])) * y[ii] -
			(pa.beta_ij/(pa.nb-1.0))*(sumBees_extra-pa.nap[ii]*y[ii]) + 
			y[ii] * 
			pa.gamma[ii]*pa.f_help[ii]/(1.0+pa.h*pa.gamma[ii]*pa.f_help[ii]);
	}
	
	for (int ii=pa.nb; ii<(pa.nf+pa.nb); ii++) {
		f[ii] = pa.alpha*y[ii] - (2.0*pa.beta_ii*pa.nap[ii] + (pa.beta_ij/(pa.nf-1.0))*(sumFlowers-pa.nap[ii])) * y[ii] -
			(pa.beta_ij/(pa.nf-1.0))*(sumFlowers_extra-pa.nap[ii]*y[ii]) + 
			y[ii] *
			pa.gamma[ii]*pa.f_help[ii]/(1.0+pa.h*pa.gamma[ii]*pa.f_help[ii]);		
	}

	

	// (4) Schleife ueber Adjazenzliste
	for (int ii=0; ii<pa.nEdges; ii++) {
		ind_a = pa.adList[2*ii];
		ind_b = pa.adList[2*ii+1];

		f[ind_a] += pa.nap[ind_b]*y[ind_b] * pa.gamma[ind_b]/ 
			( (1.0 + pa.h * pa.gamma[ind_b] * pa.f_help[ind_b]) * (1.0 + pa.h * pa.gamma[ind_b] * pa.f_help[ind_b]) );
	
		f[ind_b] += pa.nap[ind_a]*y[ind_a] * pa.gamma[ind_a]/
			( (1.0 + pa.h * pa.gamma[ind_a] * pa.f_help[ind_a]) * (1.0 + pa.h * pa.gamma[ind_a] * pa.f_help[ind_a]) );
	}
	

	return GSL_SUCCESS;
}


#elif
int netBee (double t, const double y[], double f[], void *params) {
	(void) (t);
    paraSimple pa = *(paraSimple *) params;
  
	int ind_a, ind_b;
	double sumBees, sumFlowers, mutualism; 
	sumBees = 0.0;
	sumFlowers = 0.0;
	mutualism = 0.0;

	for (int ii=0; ii<pa.nb; ii++) {
		sumBees += y[ii];
		f[ii] = 0.0;
	}

	for (int ii=pa.nb; ii<(pa.nf+pa.nb); ii++) {
		sumFlowers += y[ii];
		f[ii] = 0.0;
	}

	for (int ii=0; ii<pa.nEdges; ii++) {
		ind_a = pa.adList[2*ii];
		ind_b = pa.adList[2*ii+1];

		f[ind_a] += y[ind_b];
		f[ind_b] += y[ind_a];
	}

	for (int ii=0; ii<pa.nb; ii++) {
		mutualism = pa.gamma[ii]*f[ii]/(1.0+pa.h*pa.gamma[ii]*f[ii]);
		f[ii] = ( pa.alpha  -  pa.beta_ii*y[ii]  -  (pa.beta_ij/(pa.nb-1.0))*(sumBees-y[ii])  + 
			y[ii]/(y[ii]+pa.delta)  *  mutualism ) * y[ii];
	}

	for (int ii=pa.nb; ii<(pa.nb+pa.nf); ii++) {
		mutualism = pa.gamma[ii]*f[ii]/(1.0+pa.h*pa.gamma[ii]*f[ii]);
		f[ii] = ( pa.alpha  -  pa.beta_ii*y[ii]  -  (pa.beta_ij/(pa.nf-1.0))*(sumFlowers-y[ii])  + 
			y[ii]/(y[ii]+pa.delta)  *  mutualism ) * y[ii];
	}


	return GSL_SUCCESS;
}



int netBee_B (double t, const double y[], double f[], void *params) {
	// zumindest die Formulierung sollte jetzt stimmen
	(void) (t);
    // paraSimple pa = *(paraSimple *) params; 
	paraSimple pa = *(paraSimple *) params;

	int ind_a, ind_b;
	double sumBees, sumBees_extra, sumFlowers, sumFlowers_extra;
	
	sumBees = 0.0;
	sumBees_extra = 0.0;
	sumFlowers = 0.0;
	sumFlowers_extra = 0.0;


	// (1) Schleifen ueber Knoten
	for (int ii=0; ii<pa.nb; ii++) {
		sumBees += pa.nap[ii];
		sumBees_extra += pa.nap[ii]*y[ii];
		f[ii] = 0.0;
		pa.f_help[ii] = 0.0;
	}

	for (int ii=pa.nb; ii<(pa.nf+pa.nb); ii++) {
		sumFlowers += pa.nap[ii];
		sumFlowers_extra += pa.nap[ii]*y[ii];
		f[ii] = 0.0;
		pa.f_help[ii] = 0.0;
	}


	
	// (2) Schleife ueber Adjazenzliste
	for (int ii=0; ii<pa.nEdges; ii++) {
		ind_a = pa.adList[2*ii];
		ind_b = pa.adList[2*ii+1];

		pa.f_help[ind_a] += pa.nap[ind_b];
		pa.f_help[ind_b] += pa.nap[ind_a];
	}



	// (3) Schleifen ueber Knoten
	for (int ii=0; ii<pa.nb; ii++) {
		f[ii] = pa.alpha*y[ii]  -  (2.0*pa.beta_ii*pa.nap[ii] + (pa.beta_ij/(pa.nb-1.0))*(sumBees-pa.nap[ii])) * y[ii] -
			(pa.beta_ij/(pa.nb-1.0))*(sumBees_extra-pa.nap[ii]*y[ii]) + 
			y[ii] * ( 1.0 - pa.delta*pa.delta/((pa.nap[ii]+pa.delta)*(pa.nap[ii]+pa.delta)) ) * 
			pa.gamma[ii]*pa.f_help[ii]/(1.0+pa.h*pa.gamma[ii]*pa.f_help[ii]);
	}
	
	for (int ii=pa.nb; ii<(pa.nf+pa.nb); ii++) {
		f[ii] = pa.alpha*y[ii] - (2.0*pa.beta_ii*pa.nap[ii] + (pa.beta_ij/(pa.nf-1.0))*(sumFlowers-pa.nap[ii])) * y[ii] -
			(pa.beta_ij/(pa.nf-1.0))*(sumFlowers_extra-pa.nap[ii]*y[ii]) + 
			y[ii] * ( 1.0 - pa.delta*pa.delta/((pa.nap[ii]+pa.delta)*(pa.nap[ii]+pa.delta)) ) * 
			pa.gamma[ii]*pa.f_help[ii]/(1.0+pa.h*pa.gamma[ii]*pa.f_help[ii]);		
	}

	

	// (4) Schleife ueber Adjazenzliste
	for (int ii=0; ii<pa.nEdges; ii++) {
		ind_a = pa.adList[2*ii];
		ind_b = pa.adList[2*ii+1];

		f[ind_a] += (1.0-pa.delta/(pa.nap[ind_b]+pa.delta))*pa.nap[ind_b]*y[ind_b] * pa.gamma[ind_b]/ 
			( (1.0 + pa.h * pa.gamma[ind_b] * pa.f_help[ind_b]) * (1.0 + pa.h * pa.gamma[ind_b] * pa.f_help[ind_b]) );
	
		f[ind_b] += (1.0-pa.delta/(pa.nap[ind_a]+pa.delta))*pa.nap[ind_a]*y[ind_a] * pa.gamma[ind_a]/
			( (1.0 + pa.h * pa.gamma[ind_a] * pa.f_help[ind_a]) * (1.0 + pa.h * pa.gamma[ind_a] * pa.f_help[ind_a]) );
	}
	

	return GSL_SUCCESS;
}

#endif

