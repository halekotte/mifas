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
#include "myDifferentials.h"
#include "network.h"
#include "beegrid.h"
#include "mifas.h"

#include <fstream>
#include <string>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>


int determine_mifas(Beegrid * mouiNet, double T_ini, const bool * perturb_null, 
	std::string name_results, unsigned long complete_runs) {
	

	int t_steps, status, count, n_ini;
	double t, h;
	bool coexist_all, finIni;


	double epsilon;                // step size of adaptation
    double eps_min;                // * Default- und minimal value for eps_stern
    double max_eps;                // maximal possible absolute value forr epsilon
    double max_eps_h1, max_eps_h2; // helpers for determining max_eps
    double lambda;                 // for ratio of x_null to nu
    double wurz, nicht_wurz;       // helpers for calulating lambda

    double eps_a, eps_b, eps_c;
    double dist_b, dist_c, dist_current;
    double eps_b_ini, eps_c_ini, eps_dist_ac;
    double max_eps_fak;
    bool stay_tuned;
    double eps_help, dist_help;
    double T;

    double d;        		// * current radius (Distanz der Anfangsbedingungen zum FP)
    double d_min;         	// current minimal value of d
    double dStep;         	// step size for reducing d
    double test_laenge;   	// helper for scaling of x_null

    unsigned long jj;     	// counter

	
	double var_phase2, ini_phase;


	std::ofstream datei;
	// ----- Dateinamen in char-Arrays umwandeln:
    const char * fname = name_results.c_str();

	size_t nod = 7;

	double dStep_vec [7] = {0.1, 0.05, 0.01, 0.001, 0.05, 0.05, 0.001};
	double crit_enhance = 1.0e-6;
	double T_plus [30] =  {0.0, 0.0, 0.0, 0.0, 2.0, 3.0, 0.0, 5.0, 5.0, 5.0, 
							1.0, 1.0, 3.0, 3.0, 3.0, 1.0, 1.0, 1.0, 1.0, 1.0,
							3.0, 3.0, 3.0, 3.0, 3.0, 1.0, 1.0, 1.0, 1.0, 1.0};


	// ----- NEW:
	int count_in_basin = 0;
	const int count_in_crit [4] = {1, 1, 1, 1};



	// ===== PARAMETERS TO CHOOSE:

	h = 0.01;
	T = T_ini;

	ini_phase = 1;
	var_phase2 = 0;


	// =============== CONSTANT ================
	t_steps = int(T_ini/h);
	t = 0.0;


	// ----
    eps_b_ini = 0.01;
    eps_c_ini = 0.02;


    max_eps_fak = 0.1;
    eps_dist_ac = 1.0e-4;
    eps_min = 1.0e-5; // es ist wichtig, dass diese kleiner ist als eps_dist_ac


    // Anfangswerte fuer Parameter, die nicht extern vorgegeben werden muessen:
    d     = 1.0;
    d_min = 100.0;


	// a couple of arrays
	double * N_forward;	
	double * x_T;
	double * nu;
	double * x_null;
	double * x_null_s;
	double * x_star;
	
	
	// get parameters
	struct paraSimple * pams = mouiNet->createParams();
	unsigned long int dimVar = mouiNet->getNumberNodes();


	// WICHTIG: Sollte ich die Integrationszeit erhoehen, muss ich auch die 
	// Groesse des Arrays N_forward veraendern.
	N_forward = new double[dimVar*t_steps];
	x_T = new double[dimVar];
	nu = new double[dimVar];
	x_null = new double[dimVar];
	x_null_s = new double[dimVar];
	x_star = new double[dimVar];

 
	// Arrays initialisieren (das brauch man im Prinzip nicht zwingend)
	for (int ii_x=0; ii_x<dimVar; ii_x++) {
		x_star[ii_x] = 0.0;
		x_T[ii_x] = 0.0;
		nu[ii_x] = 0.0;
		x_null[ii_x] = 0.0;
		x_null_s[ii_x] = 0.0;
	}

	for (int ii_x=0; ii_x<(dimVar*t_steps); ii_x++) {
		N_forward[ii_x] = 0.0;
	}

	// determine desired state
	coexist_all = determineStar(pams, dimVar, x_star);	

	//for (int ii=0; ii<dimVar; ii++)
	//	std::cout << ii << ' '<< x_star[ii] << '\n';
	//std::cout << '\n';


	if (!coexist_all) {
		std::cout << "No fixed point with all species coexisting." << std::endl;
		complete_runs = 0;
	}


	// initialize integrators
	int (*func_for)(double, const double *, double *, void *) = netBee;
	gsl_odeiv2_system sys_for = {func_for, 0, dimVar, pams};

	int (*func_back)(double, const double *, double *, void *) = netBee_B;
	gsl_odeiv2_system sys_back = {func_back, 0, dimVar, pams};

	const gsl_odeiv2_step_type * st_for = gsl_odeiv2_step_rkf45;
	gsl_odeiv2_step * s_for = gsl_odeiv2_step_alloc(st_for, dimVar);
	gsl_odeiv2_control * c_for = gsl_odeiv2_control_y_new(1e-8, 1e-8);
	gsl_odeiv2_evolve * e_for = gsl_odeiv2_evolve_alloc(dimVar);

	// Initialisierung der Zufallszahlengeneratoren (uniform und normal distribution)
	std::random_device rd;
	std::mt19937 engine(rd());
	std::uniform_real_distribution<double> distribution_uniform(0.0,1.0);
	std::normal_distribution<double> distribution_normal(0.0,1.0);

	// ==============================================
	
	
	// Hilfsvariable initialisieren (ist auch nichz zwingend notwendig)
	for (int ii=0; ii<dimVar; ii++) {
		pams->nap[ii] = 0.0;
	}	


	// =============== CENTRAL PART ============== 

	datei.open(fname,std::ios::out);

	n_ini = 0;

    // ===== (3) Schleife ueber mehrere Durchlaeufe:
    for (size_t ii = 0; ii < complete_runs; ii++) {

		// Integrationszeit zuruecksetzen
    	T = T_ini;


        // ===== (4) Stochastic Initialization (Prerun):
		if (ini_phase == 1) {
			d = stochastic_ini(&sys_for, s_for, c_for, e_for, dimVar, h, engine, distribution_uniform,
				distribution_normal, x_null, perturb_null, x_T, x_star, n_ini);

			// std::cout << n_ini << std::endl;

			if (ii==0 && n_ini>25000) {
				datei << 100;
				break;
			}

			if (ii > 10000) {
				ini_phase = 2;
			}
		
		} else if (ini_phase == 2) {
			std::cout << "Phase 2" << std::endl;
			d = 2*d_min;
			while (d > 1.01*d_min && var_phase2<dimVar) {
				d = unstochastic_ini(&sys_for, s_for, c_for, e_for, dimVar, h, x_null, x_T, x_star, var_phase2);
				var_phase2++;
			}
			if (var_phase2 >= dimVar) {
				ini_phase = 3;
			}

		} else {
			d = stochastic_ini(&sys_for, s_for, c_for, e_for, dimVar, h, engine, distribution_uniform,
				distribution_normal, x_null, perturb_null, x_T, x_star, n_ini);

		}
        std::cout << d << std::endl;

		finIni = false;

        // ===== (4) Ende (siehe Funktion unten)



        // ===== (5) Loop for descreasing d
        for (size_t ii_s=0; ii_s<nod; ii_s++) {   

        	dStep = dStep_vec[std::min(int(ii_s), 6)];
            T += T_plus[ii_s];

			// Wichtig Extra-Speicherplatz zur Verfuegung zu stellen
			if (T_plus[ii_s]>0) {
				t_steps = int(T/h);
				delete [] N_forward;
				N_forward = new double[dimVar*t_steps];
			}

			// Variablen fuer Abbruchkriterium zuruecksetzen
        	stay_tuned = true;
			count_in_basin = 0;


			// ===== OPTIMIZATION =====
			while (stay_tuned and d > 0.0) {

				dist_b = 1.0e5;
				dist_current = 1.0e-5;


				// ===== (6) method according to Kerswell
				// criterion for terminating the process
			    while ( (dist_b/dist_current) > (1.0+crit_enhance) ) {


				    // !!!!!!!!!!!!!!!!!!!!!!
				    test_laenge = 0.0;
				    for (jj=0; jj<dimVar; jj++) {
					    if (perturb_null[jj]) {
						    test_laenge += x_null[jj]*x_null[jj];
					    }
				    }
				    test_laenge = sqrt(test_laenge);


				    // Initialisirung von x_T und x_null auf die korrekte neue Laenge bringen:
				    for (jj = 0; jj<dimVar; jj++) {
					    if (perturb_null[jj]) {
						    x_null[jj] *= (d/test_laenge);
						    x_T[jj] = x_star[jj] + x_null[jj];
					    } else {
                            x_null[jj] = 0.0;
						    x_T[jj] = x_star[jj];
					    }
				    }
				    // !!!!!!!!!!!!!!!!!!!!!!


					// Sole backward integration
					dist_b = integrate_forward_backward(pams, &sys_for, &sys_back, s_for, c_for, e_for, 
						x_T, nu, x_star, N_forward, dimVar, T, h);	
					dist_current = dist_b;

					// -> zu benutzen sind: nu (entspricht nu(0)), x_null (entspricht der Stoerung)
					// x_T hingegen wird weiterhin nur als Hilfsvariable fuer Simulationen verwendet
					// x_star sollte ohnehin in keinem Fall veraendert werden




					// ===== (9) Adapt step size:
					// ----- nu anpassen und maximalen Wert fuer epsilon bestimmen:
					
					max_eps_h1 = 0.0;
					max_eps_h2 = 0.0;

					for (jj = 0; jj<dimVar; jj++) {
						if (!perturb_null[jj]) {
							nu[jj] = 0.0;
						} else {
							max_eps_h1 += nu[jj]*nu[jj];
							max_eps_h2 += x_null[jj]*nu[jj];
						}
					}
					max_eps = max_eps_fak * std::sqrt( 1.0 / ( max_eps_h1/(d*d) - (max_eps_h2*max_eps_h2)/(d*d*d*d) ) );

					// ========================


					// einige Hilfsvariablen initiieren
					eps_a = 0.0;
					dist_b = 0.0;
					eps_b = eps_b_ini;
					eps_c = eps_c_ini;
					dist_c = -1.0;


					// ===== (10) Loop for adapting step size:
					// simply serach for the step size leading to the highest distance

                    // DINGDINGDINGDINGDINGDING --->>> (A)
					// right border
					while (dist_b<=dist_current and eps_b > eps_min) {

						dist_b =  det_distance(pams, &sys_for, s_for, c_for, e_for,
							dimVar, T, h, eps_b, max_eps, d, x_null, nu, max_eps_h1, max_eps_h2, 
							perturb_null, x_T, x_star);

						if (dist_b <= dist_current) {
							eps_c = eps_b;
							dist_c = dist_b;
							eps_b *= 0.5;
						}
					}

					if (dist_c<0.0 and eps_c < 1.0) {
						dist_c = dist_b+1.0;
					}


					// Sollte das Maximum in dieser Richtung bereits erreicht sein
					if (eps_b<=eps_min) {
						break;
					}


					// DINGDINGDINGDINGDINGDING --->>> (B)
					// left border
                    while (dist_c >= dist_b and eps_c < 1.0) {

						dist_c =  det_distance(pams, &sys_for, s_for, c_for, e_for,
							dimVar, T, h, eps_c, max_eps, d, x_null, nu, max_eps_h1, max_eps_h2, 
							perturb_null, x_T, x_star);

                    	if (dist_c >= dist_b) {
                    		eps_a = eps_b;
                    		dist_b = dist_c;
                    		eps_b = eps_c;
                    		eps_c += 0.01;
                    	}

                    }


                    // DINGDINGDINGDINGDINGDING --->>> (C)
					// narrow it down
                    while ( (eps_c-eps_a) > eps_dist_ac ) {

                    	// --- Check left side:
                        eps_help = (eps_a+eps_b)/2.0;

						dist_help =  det_distance(pams, &sys_for, s_for, c_for, e_for,
							dimVar, T, h, eps_help, max_eps, d, x_null, nu, max_eps_h1, max_eps_h2, 
							perturb_null, x_T, x_star);

                        if (dist_help > dist_b) {
                        	eps_c = eps_b;
                        	eps_b = eps_help;
                        	dist_b = dist_help;
                        } else {
                        	eps_a = eps_help;
                        }


                        // --- Check right side:
                        eps_help = (eps_b+eps_c)/2.0;

						dist_help =  det_distance(pams, &sys_for, s_for, c_for, e_for,
							dimVar, T, h, eps_help, max_eps, d, x_null, nu, max_eps_h1, max_eps_h2, 
							perturb_null, x_T, x_star);

                        if (dist_help > dist_b) {
                        	eps_a = eps_b;
                        	eps_b = eps_help;
                        	dist_b = dist_help;
                        } else {
                        	eps_c = eps_help;
                        }

                    }


                    // --->>> set new x_null					
                    epsilon = eps_b*max_eps;

                    nicht_wurz = max_eps_h2/(d*d) - 1.0/epsilon;
                    wurz = 1.0/(epsilon*epsilon) + max_eps_h2*max_eps_h2/(d*d*d*d) - max_eps_h1/(d*d);
                    lambda = nicht_wurz + std::sqrt(wurz);

                    for (jj=0; jj<dimVar; jj++) {
                    	if (perturb_null[jj]) {
                    		x_null[jj] = x_null[jj] + epsilon*(lambda*x_null[jj]-nu[jj]);
                    	} else {
                    		x_null[jj] = 0.0;
                    	}
                    }


				} // === (6) Ende



				// ----- Abbruch Schleife (5)
				for (jj=0; jj<dimVar; jj++) {
					x_T[jj] = x_star[jj] + x_null[jj];				
				}


				// NEW: Mit Erreichen des Basins des FP wird die Suche nicht zwingend abgebrochen.
				if (isIn(x_T, dimVar, x_star, &sys_for, s_for, c_for, e_for, h)) {

					count_in_basin++;

				} else {
					d_min = d<d_min ? d:d_min;
	
					for (jj=0; jj<dimVar; jj++) {
						x_null_s[jj] = x_null[jj];
					}
					count_in_basin = 0;

					if (!finIni) {
						finIni = true;
					}
				}

				if (count_in_basin>=count_in_crit[std::min(int(ii_s), 3)]) {
					stay_tuned = false;
					
					for (jj=0; jj<dimVar; jj++) {
					    x_null[jj] = x_null_s[jj];
				    }
				}

				d -= dStep;

			} // === (5) Ende


			if (finIni) {
				d += (1.0+count_in_crit[std::min(int(ii_s), 3)])*dStep;
            	std::cout << "Lauf " << ii+1 << " - " << ii_s+1 << "/" << nod << ": " << d << std::endl;
			} else {
				ii_s = nod;
				// std::cout << "Still initializing." << std::endl;
			}
        }


        // Hier muss der Abschnitt zur Speicherung der Daten hin:
		if (finIni) {
        	datei << ii+1 << ' ' << d-dStep << ' ' << d << ' ';
        	for (jj=0; jj<dimVar; jj++) {
        		if (perturb_null[jj]) {
        		    datei << x_null_s[jj] << ' ';
        		} else {
        			datei << 0.0 << ' ';
        		}
        	}
        	datei << '\n';
			n_ini = 0;
		} else {
			ii--;
			n_ini++;
		}


    } // === (3) Ende

    datei.close();
	
	// ============================================

	delete [] x_T;
	delete [] nu;
	delete [] x_null;
	delete [] x_null_s;
	delete [] x_star;
	delete [] N_forward;

}


// ===========================================================================================================
// -----> function for the stochastic initialization
double stochastic_ini(gsl_odeiv2_system * sys_for,gsl_odeiv2_step * s_for, 
	gsl_odeiv2_control * c_for, gsl_odeiv2_evolve * e_for, const unsigned long dimVar, const double h, 
	std::mt19937 &engine, std::uniform_real_distribution<double> &distribution_uniform, 
	std::normal_distribution<double> &distribution_normal,
	double * x_null, const bool * perturb_null, double * x_run, const double * xStar, const int n_ini) {
	
	// in x_null soll die Stoerung geschrieben werden
	// x_run ist fuer die Simulation gedacht
	// perturb_null zeigt an welche Knoten ueberhaupt gestoert werden sollen

	// In dieser Routine sollen nun mit zunehmender Dauer der Initialisierung groessere
	// Radien wahrscheinlicher werden.

	double x_help, d, r_min_squared, r_max_squared, r_min, r_max, test_laenge, n_dim;

	n_dim = min(1.0 + n_ini/10.0, double(dimVar));
	n_dim = min(n_dim,100.0);
	// std::cout << n_dim << std::endl;

	r_min = 0.1*xStar[0];
	r_max = 0.999*xStar[0];
	r_min_squared = pow(r_min,n_dim);
	r_max_squared = pow(r_max,n_dim);
	
	d = pow( distribution_uniform(engine) * (r_max_squared - r_min_squared) + r_min_squared, 1.0/n_dim ); 
	d = floor(d*1000.0)/1000.0;


	// eigentlich unnoetig, aber safety first und so
	while (d>=xStar[0]) {
		d *= 0.99;
	}

	test_laenge = 0.0;
	for (int ii=0; ii<dimVar; ii++) {
		if (perturb_null[ii]) {
			x_help = (distribution_normal(engine));
			x_null[ii] = x_help;
			test_laenge += x_help*x_help;
		} else {
			x_null[ii] = 0.0;
		}
	}


	test_laenge = sqrt(test_laenge);
	for (int ii=0; ii<dimVar; ii++) {
		x_null[ii] = x_null[ii]*d/test_laenge;
		x_run[ii] = xStar[ii] + x_null[ii];
	}

	return d;

}



double unstochastic_ini(gsl_odeiv2_system * sys_for,gsl_odeiv2_step * s_for, 
	gsl_odeiv2_control * c_for, gsl_odeiv2_evolve * e_for, const unsigned long dimVar, const double h, 
	double * x_null, double * x_run, const double * xStar, const int n_var) {
	
	double shock_out, shock_in, shock_test;
	bool in;


	shock_in = 0.025;
	shock_out = xStar[0] - 0.001;
		
	while (shock_out-shock_in > 0.1) {
		shock_test = (shock_in+shock_out)/2.0;

		for (int ii=0; ii<dimVar; ii++) {
			x_run[ii] = xStar[ii];
		}
		x_run[n_var] -= shock_test;
		

		in = isIn(x_run, dimVar, xStar, sys_for, s_for, c_for, e_for, h);

		if (in) {
			shock_in = shock_test;
		} else {
			shock_out = shock_test;
		}
	
	}



	for (int ii=0; ii<dimVar; ii++) {
		x_null[ii] = 0.0;
	}
	x_null[n_var] = shock_out;


	return shock_out;

}


// ===========================================================================================================
// -----> function to determine distance (objective function)
double det_distance(void * params, gsl_odeiv2_system * sys_for,
	gsl_odeiv2_step * s_for, gsl_odeiv2_control * c_for, gsl_odeiv2_evolve * e_for,
	unsigned long dimVar, double T, double h,
	double eps_new, double max_eps, double d, double * x_null, double * nu,
	double max_eps_h1, double max_eps_h2, const bool * perturb_null, double * x_run, double * xStar) {

	double dist_new, epsilon, nicht_wurz, wurz, lambda;

	epsilon = eps_new*max_eps;

	nicht_wurz = max_eps_h2/(d*d) - 1.0/epsilon;
	wurz = 1.0/(epsilon*epsilon) + max_eps_h2*max_eps_h2/(d*d*d*d) - max_eps_h1/(d*d);
	lambda = nicht_wurz + sqrt(wurz);

	// std::cout << "nicht " << nicht_wurz << "   wurz " << sqrt(wurz) << std::endl;  

	for (size_t ii=0; ii<dimVar; ii++) {
		if (perturb_null[ii]) {
			x_run[ii] = xStar[ii] + x_null[ii] + epsilon*(lambda*x_null[ii]-nu[ii]);
		} else {
		    x_run[ii] = xStar[ii] + x_null[ii];
		}
	}

	dist_new = integrate_forward(params, sys_for, s_for, c_for, e_for, x_run, xStar, dimVar, T, h);


	return dist_new;
}



// ===========================================================================================================
// -----> function for forwrad- and backward-integration
double integrate_forward_backward(void * params, 
	gsl_odeiv2_system * sys_for, gsl_odeiv2_system * sys_back, gsl_odeiv2_step * s_for, gsl_odeiv2_control * c_for, 
	gsl_odeiv2_evolve * e_for, double * y, double * nu, const double * xStar, double * N_forward, 
	const unsigned long dimVar, const double T, const double h) {

	int count, t_steps, status; 
	double dist, t, y_help;
	y_help = 0.0;
	t_steps = int(T/h);

	dist = 0.0;
	t = 0.0;
	
	paraSimple pa = *(paraSimple *) params;


	// ----- forward-integration:
	count = 0;
	for (int ii=0; ii<t_steps; ii++) {	
		status = gsl_odeiv2_evolve_apply_fixed_step(e_for, c_for, s_for, sys_for, &t, h, y);
		for (int jj=dimVar-1; jj>=0; jj--) {
			N_forward[count] = y[jj];
			count++;
		}
	}


	// ----- Anfangsbedingungen fuer Rueckwaertsintegration:
	for (int jj=0; jj<dimVar; jj++) {
		y_help = y[jj]-xStar[jj];
		dist += y_help*y_help;
		nu[jj] = -2.0*y_help; 
	}


	// ----- backward-integration:
	t = 0.0;
	for (int ii=0; ii<t_steps; ii++) {
		for (int jj=0; jj<dimVar; jj++) {
			count--;
			pa.nap[jj] = N_forward[count];
		}
		status = gsl_odeiv2_evolve_apply_fixed_step(e_for, c_for, s_for, sys_back, &t, h, nu);
	}


	return dist;
}



// =============================================================================================================
// -----> function for forward-integration
double integrate_forward(void * params, 
	gsl_odeiv2_system * sys_for, gsl_odeiv2_step * s_for, gsl_odeiv2_control * c_for, gsl_odeiv2_evolve * e_for,
	double * y, const double * xStar, const unsigned long dimVar, const double T, const double h) {

	int t_steps, status; 
	double dist, t;
	t_steps = T/h;

	dist = 0.0;
	t = 0.0;
	
	// ----- Vorwaertsintegration:
	for (int ii=0; ii<t_steps; ii++) {
		status = gsl_odeiv2_evolve_apply_fixed_step(e_for, c_for, s_for, sys_for, &t, h, y);
	}

	// ----- Anfangsbedingungen fuer Rueckwaertsintegration:
	for (int jj=0; jj<dimVar; jj++) {
		dist += (y[jj]-xStar[jj])*(y[jj]-xStar[jj]); 
	}


	return dist;
}



// ===============================
// ===============================
// ===============================
bool determineStar(struct paraSimple * params, long unsigned int dimVar, double * xStar) {

	double t = 0.0;
	double T = 10000.0;
	double min_xStar = 100.0;
	double max_xStar = 0.0;


	int (*func)(double, const double *, double *, void *) = netBee;

	gsl_odeiv2_system sys = {func, 0, dimVar, params};
	gsl_odeiv2_driver * ddd = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, 1e-8, 1e-8, 1e-8);

	for (unsigned long int ii=0; ii<dimVar; ii++) {
		xStar[ii] = 5.0 + 0.01*ii;
	}

	int status = gsl_odeiv2_driver_apply (ddd, &t, T, xStar);


	for (int ii=0; ii<dimVar; ii++) {
		if (xStar[ii]<min_xStar) {
			min_xStar = xStar[ii];
		}
		if (xStar[ii]>max_xStar) {
			max_xStar = xStar[ii];
		}
	}

	return (min_xStar>0.1) && (max_xStar-min_xStar<0.001); 
}


// ===============================
bool isIn(double * y, const long unsigned int dimVar, const double * xStar, gsl_odeiv2_driver * ddd) {

	double t = 0.0;
	double Tinter = 10.0;
	double T = Tinter;
	double bee_min = 10.0; 
	double dist_all = 0.0;
	double test_dist = 0.0;
	int count_in = 0;


	while (abs(count_in) < 3 && T < 100000) {
		int status = gsl_odeiv2_driver_apply (ddd, &t, T, y);

		dist_all = 0.0;
		bee_min = xStar[0];

		for (int ii=0; ii<dimVar; ii++) {
			test_dist = y[ii]-xStar[ii];
			test_dist *= test_dist;
			dist_all += test_dist;
			if (y[ii] < bee_min) {
				bee_min = y[ii]; 			
			}
		}

		dist_all = sqrt(dist_all);

		if (bee_min<0.05) {
			count_in--;
		} else if (dist_all<0.01) {
			count_in++;
		} else {
			count_in = 0;
		}

		T += Tinter;
	}

	return count_in>0;
}



// ===============================
// this function aims at determining whether a point is in the basin for the actual 
// mifas-algorithm n which a constant step width is used, one might think about using 
// the simpler approach of varying steps sizes in order to reduce time consumption
bool isIn(double * y, const long unsigned int dimVar, const double * xStar, gsl_odeiv2_system * sys_for,
	gsl_odeiv2_step * s_for, gsl_odeiv2_control * c_for, gsl_odeiv2_evolve * e_for, const double h) {

	double t = 0.0;
	double Tinter = 10.0;
	double bee_min = 10.0; 
	double dist_all = 0.0;
	double test_dist = 0.0;
	int count_in = 0;
	double t_steps = Tinter/h;
	int status;
	double dist_crit, dist_ini;

	dist_ini = 0.0;

	for (int ii=0; ii<dimVar; ii++) {
		test_dist = y[ii]-xStar[ii];
		test_dist *= test_dist;
		dist_ini += test_dist;
	}
	
	dist_ini = sqrt(dist_ini);

	// Bestimme Anfangsabstand und waehle entsprechende Abbruchbedingung:
	dist_crit = 0.1*dist_ini<0.05 ? 0.1*dist_ini : 0.05; 



	while (abs(count_in) < 3 && t < 100000) {
		// ----- Vorwaertsintegration:
		for (int ii=0; ii<t_steps; ii++) {
			status = gsl_odeiv2_evolve_apply_fixed_step(e_for, c_for, s_for, sys_for, &t, h, y);
		}

		// std::cout << status << ' ' << t << ' ' << y[0] << std::endl;

		if (status!=GSL_SUCCESS) {
			std::cerr << "REASONABLE ERROR: Schrittweite der Integration ist zu grob." << std::endl;
			exit(23);
		}

		dist_all = 0.0;
		bee_min = xStar[0];

		for (int ii=0; ii<dimVar; ii++) {
			test_dist = y[ii]-xStar[ii];
			test_dist *= test_dist;
			dist_all += test_dist;
			if (y[ii] < bee_min) {
				bee_min = y[ii]; 			
			}
		}

		dist_all = sqrt(dist_all);

		if (bee_min<0.05) {
			count_in--;
		} else if (dist_all<dist_crit) {
			count_in++;
		} else {
			count_in = 0;
		}

	}

	// std::cout << "Still: Drin oder nicht drin ..." << std::endl;

	return count_in>0;
}
