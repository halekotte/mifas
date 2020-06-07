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



#include "mifas_new.h"
#include "myDifferentials.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include <random>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>



// =========== HAUPTROUTINE ~~~~~~~~~~~~~~~~~~~~~
void determine_mifas_new(void * params, int (*func)(double, const double *, double *, void *),
		                         int (*func_back)(double, const double *, double *, void *), unsigned long dimVar,
		                         const bool * perturb_null, const bool * var_norm,
								 std::string name_results, unsigned long complete_runs) {


	// ----- ========== DE-KLARA-TION ========== -----

	// ===== (1) Parameter:
    double epsilon;                // step size of adaptation
    double eps_min;                // default- und minimal value of eps_stern
    double max_eps;                // maximal absolute value of epsilon
    double max_eps_h1, max_eps_h2; // helpers for determining max_eps
    double lambda;                 // factor for relation of x_null to nu
    double wurz, nicht_wurz;       // helpers for dtermining lambda

    double eps_a, eps_b, eps_c;
    double dist_b, dist_c, dist_current;
    double eps_b_ini, eps_c_ini, eps_dist_ac;
    double max_eps_fak;
    bool stay_tuned;
    double eps_help, dist_help;
    double T, T_ini;

    double d;             // current radius
    double d_min;         // minimal value of d, so far
    double dStep;         // step size for descressing d
    double test_laenge;   // helper for scaling x_null

    unsigned long jj;     // counter


	// === NEW INI === 
	// choose which initialization to use
	bool OLD_INI = true;
	bool finIni = false;	
	int n_ini;
	


	std::ofstream datei;

	// we need to increase the integration times to suffieciently large values
	size_t nod = 10;

	double dStep_vec [7] = {0.1, 0.05, 0.01, 0.001, 0.001, 0.001, 0.0001};
	double crit_enhance = 1.0e-5;
	double T_plus [30] =  {0.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 
							1.0, 1.0, 3.0, 3.0, 3.0, 1.0, 1.0, 1.0, 1.0, 1.0,
							3.0, 3.0, 3.0, 3.0, 3.0, 1.0, 1.0, 1.0, 1.0, 1.0};



	// ----- NEU:
	int count_in_basin = 0;
	const int count_in_crit [4] = {3, 3, 10, 15}; // {3, 3, 10, 15};

	bool phase_1 = true;
	int start_phase_2 = 100; //std::max(10, int(complete_runs)-10);
	


    // ----- Arrays:
    double *x_null, *x_T, *x_null_s, *nu, *x_star, *y, *x_null_bestguess;

    x_null   = new double[dimVar];        // dynamical initial conditions
    x_T      = new double[dimVar];        // storage address for numerical simulations

    x_null_s = new double[dimVar];        // temporary storage
	x_null_bestguess = new double[dimVar]; // best guess for each run 

    nu       = new double[dimVar];        // derived state variables

    x_star = new double[dimVar];		  // fixed point

	y = new double[dimVar*2];			  // backwards


    // ----- Dateinamen in char-Arrays umwandeln:
    const char * fname = name_results.c_str();


    // ----- ========== DEFINITION ========== -----
    eps_b_ini = 0.01;
    eps_c_ini = 0.02;


    max_eps_fak = 0.1;
    eps_dist_ac = 1.0e-3;
    eps_min = 1.0e-3; 


    // starting values
    d     = 1.0;
    d_min = 100.0;
	T_ini = 10.0;


    // load fixed point
    paraBullo pa = *(paraBullo *) params;

	for (jj=0; jj<dimVar; jj++) {
		y[jj] = 0.0;
	}

    for (jj=0; jj<dimVar; jj++) {
    	if (jj < dimVar/2) {
            x_star[jj] = pa.omega[jj];
    	} else {
    	    x_star[jj] = 0.0;
    	}
    }

    std::cout << " --- Looking for MiFaS --- " << std::endl;

    // ----- ================================= -----



    // ----- ========== INITIALIZATION ========== -----

	// ===== (2a) Initialization of the diff.-equations:
    // Forwards (a):
	gsl_odeiv2_system sys_for = {func, 0, dimVar, params};
    gsl_odeiv2_driver * d_for = gsl_odeiv2_driver_alloc_y_new(&sys_for, gsl_odeiv2_step_rkf45, 1e-6, 1e-8, 1e-8);
    // gsl_odeiv2_driver_set_hmax (d_for, 0.01);

    // Forwards (b):
    const gsl_odeiv2_step_type * st_for = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step * s_for = gsl_odeiv2_step_alloc(st_for, dimVar);
    gsl_odeiv2_control * c_for = gsl_odeiv2_control_y_new(1e-3, 1e-3);
    gsl_odeiv2_evolve * e_for = gsl_odeiv2_evolve_alloc(dimVar);

    // Backwards:
	gsl_odeiv2_system sys_back = {func_back, 0, 2*dimVar, params};
    gsl_odeiv2_driver * d_back = gsl_odeiv2_driver_alloc_y_new(&sys_back, gsl_odeiv2_step_rkf45, 1e-6, 1e-8, 1e-8);
    // gsl_odeiv2_driver_set_hmax (d_back, 0.01);

    // Backwards (alternative):
    const gsl_odeiv2_step_type * st_back = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step * s_back = gsl_odeiv2_step_alloc(st_back, 2*dimVar);
    gsl_odeiv2_control * c_back = gsl_odeiv2_control_y_new(1e-3, 1e-3);
    gsl_odeiv2_evolve * e_back = gsl_odeiv2_evolve_alloc(2*dimVar);


    // ===== (2b) Initialisierung des Zufallszahlengenerators:
	std::random_device rd;
	std::mt19937 engine(rd());
	std::uniform_real_distribution<double> distribution_uniform(0.0,1.0);
	std::normal_distribution<double> distribution(0.0,1.0);

	// ----- ===================================== -----


	datei.open(fname,std::ios::out);


	// === NEW INI === 
	n_ini = 0;

    // ===== (3) Loop over several runs:
    for (size_t ii = 0; ii < complete_runs; ii++) {

		T = 5.0 + 15.0*distribution_uniform(engine);
		T = floor(T*100.0)/100.0;	
		T_ini = T;

		std::cout << "T_ini = " << T << std::endl;

		d_min = 100.0;


		// === NEW INI ===
		if (ii==0 && n_ini>complete_runs*1000) {
			datei << 100;
			break;
		}



        // ===== (4) Stochastische Initialisierung (Vorlauf):
    	for (jj=0; jj<dimVar; jj++) {
    		x_null[jj] = x_star[jj];
            x_T[jj] = x_null[jj];
			x_null_bestguess[jj] = 0.0;
    	}


		if (OLD_INI) { 

			
        	d = stochastic_ini_shrinkingsphere(&params, func, d_for, engine, distribution_uniform, distribution, x_null, x_T, 
				perturb_null, dimVar);
				

			finIni = true;

		} else {

			d = stochastic_cool_and_easy(dimVar, engine, distribution_uniform, distribution, 
				x_null, perturb_null, n_ini);

			// === NEW INI ===
			finIni = false;

		}

		std::cout << d << std::endl;


        // ===== (4) Ende (siehe Funktion unten)



        // ===== (5) Loop for decreasing d
		size_t ii_s = 0;
		while (T < 25.0 and ii_s<nod) {

			dStep = dStep_vec[std::min(int(ii_s), 6)];
			// T += T_plus[ii_s];
			if (ii_s > 3) {
				T += 5.0;
			}

        	stay_tuned = true;
			count_in_basin = 0;

			while (stay_tuned and d > 0.0) {

				dist_b = 1.0e5;
				dist_current = 1.0e-5;

				// ===== (6) optimization according to Kerswell

				int count_to_five = 0;
			    while ( (dist_b/dist_current) > (1.0+crit_enhance) && count_to_five < 15) {

					count_to_five++;

				    // !!!!!!!!!!!!!!!!!!!!!!

				    test_laenge = 0.0;
				    for (jj=0; jj<dimVar; jj++) {
					    if (perturb_null[jj]) {
						    test_laenge += x_null[jj]*x_null[jj];
					    }
				    }
				    test_laenge = std::sqrt(test_laenge);


				    // Initialization of x_T und x_null and scaling them:
				    for (jj = 0; jj<dimVar; jj++) {
					    if (perturb_null[jj]) {
						    x_null[jj] *= (d/test_laenge);
						    x_T[jj] = x_null[jj];
					    } else {
                            x_null[jj] = x_star[jj];
						    x_T[jj] = x_star[jj];
					    }
				    }

				    // !!!!!!!!!!!!!!!!!!!!!!



				    // ===== (7) First Forward-Integration:
				    dist_b = integrate_forward_second_new(params, func, &sys_for, s_for, c_for, e_for, x_T, dimVar, T, var_norm);
	                dist_current = dist_b;


					for (jj = 0; jj<dimVar; jj++) {
						nu[jj] = 0.0;
					}


					// ===== (8) Backwards-Integration:
					integrate_backward_second_new(params, func_back, &sys_back, s_back, c_back, e_back, x_T, nu, dimVar, T, y);


					// ===== (9) Adaptation of step size:
					// ----- adapt nu and determine maximal value for epsilon:
					max_eps_h1 = 0.0;
					max_eps_h2 = 0.0;

					for (jj = 0; jj<dimVar; jj++) {
						if (!perturb_null[jj]) {
							nu[jj] = 0;
						} else {
							max_eps_h1 += nu[jj]*nu[jj];
							max_eps_h2 += x_null[jj]*nu[jj];
						}
					}
					max_eps = max_eps_fak * std::sqrt( 1.0 / ( max_eps_h1/(d*d) - (max_eps_h2*max_eps_h2)/(d*d*d*d) ) );




					// ========================

					eps_a = 0.0;
					dist_b = 0.0;
					eps_b = eps_b_ini;
					eps_c = eps_c_ini;
					dist_c = -1.0;


					// ===== (10) Loop for step size adaptation:


                    // DINGDINGDINGDINGDINGDING
					while (dist_b<=dist_current and eps_b > eps_min) {

						dist_b = det_distance(params, func, &sys_for,
								s_for, c_for, e_for, dimVar, T, var_norm,
								eps_b,max_eps,d,x_null,nu,
								max_eps_h1,max_eps_h2,perturb_null,x_T);

						if (dist_b <= dist_current) {
							eps_c = eps_b;
							dist_c = dist_b;
							eps_b *= 0.5;
						}
					}

					if (dist_c<0.0 and eps_c < 1.0) {
						dist_c = dist_b+1.0;
					}


					// maximum is already reached
					if (eps_b<=eps_min) {
						break;
					}

					// DINGDINGDINGDINGDINGDING
                    while (dist_c >= dist_b and eps_c < 1.0) {

                    	dist_c = det_distance(params, func, &sys_for,
                    	    s_for, c_for, e_for, dimVar, T, var_norm,
                    	    eps_c, max_eps, d, x_null, nu,
                    	    max_eps_h1, max_eps_h2, perturb_null, x_T);

                    	if (dist_c >= dist_b) {
                    		eps_a = eps_b;
                    		dist_b = dist_c;
                    		eps_b = eps_c;
                    		eps_c += 0.01;
                    	}

                    }


                    // DINGDINGDINGDINGDINGDING
                    while ( (eps_c-eps_a) > eps_dist_ac ) {

                    	// --- Check left side:
                        eps_help = (eps_a+eps_b)/2.0;

                        dist_help = det_distance(params, func, &sys_for,
                            s_for, c_for, e_for, dimVar, T, var_norm,
                            eps_help, max_eps, d, x_null, nu,
                            max_eps_h1, max_eps_h2, perturb_null, x_T);

                        if (dist_help > dist_b) {
                        	eps_c = eps_b;
                        	eps_b = eps_help;
                        	dist_b = dist_help;
                        } else {
                        	eps_a = eps_help;
                        }


                        // --- Check right side:
                        eps_help = (eps_b+eps_c)/2.0;

                        dist_help = det_distance(params, func, &sys_for,
                            s_for, c_for, e_for, dimVar, T, var_norm,
                            eps_help, max_eps, d, x_null, nu,
                            max_eps_h1, max_eps_h2, perturb_null, x_T);

                        if (dist_help > dist_b) {
                        	eps_a = eps_b;
                        	eps_b = eps_help;
                        	dist_b = dist_help;
                        } else {
                        	eps_c = eps_help;
                        }

                    }

                    // set new x_null
                    epsilon = eps_b*max_eps;

                    nicht_wurz = max_eps_h2/(d*d) - 1.0/epsilon;
                    wurz = 1.0/(epsilon*epsilon) + max_eps_h2*max_eps_h2/(d*d*d*d) - max_eps_h1/(d*d);
                    lambda = nicht_wurz + std::sqrt(wurz);

                    for (jj=0; jj<dimVar; jj++) {
                    	if (perturb_null[jj]) {
                    		x_null[jj] = x_null[jj] + epsilon*(lambda*x_null[jj]-nu[jj]);
                    	} else {
                    		x_null[jj] = x_star[jj];
                    	}
                    }


				} // === (6) Ende



				// ----- terminate loop (5)
				for (jj=0; jj<dimVar; jj++) {
					x_T[jj] = x_null[jj];				
				}
                


				// The arrival at the desired basin does not necessary mean that the search is over
				if (inside_basin_second_new(&params, func, d_for, x_T, dimVar, 0)) {

					// std::cout << "got in" << std::endl;

					count_in_basin++;
				} else {
	
					if (d < d_min) {
						d_min = d;
						for (jj=0; jj<dimVar; jj++) {
							x_null_bestguess[jj] = x_null[jj];
						}
					}
	
					for (jj=0; jj<dimVar; jj++) {
						x_null_s[jj] = x_null[jj];
					}
					count_in_basin = 0;


					// === NEW INI ===
					if (!finIni) {
						finIni = true;
						std::cout << "Ini is finished" << std::endl;
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



			// === NEW INI ===
			if (finIni) {
				d += (1.0+count_in_crit[std::min(int(ii_s), 3)])*dStep;
            	std::cout << "Lauf " << ii+1 << " - (T=" << T << "): " << d << std::endl;

				ii_s++;
			} else {
				ii_s = nod;
				// std::cout << "Still initializing." << std::endl;
			}

        }

		// === NEW INI ===
		// Save results.
		if (finIni) {
        	datei << T_ini << ' ' << d_min-dStep << ' ' << d_min << ' ';
        	for (jj=0; jj<dimVar; jj++) {
        		if (perturb_null[jj]) {
        		    datei << x_null_bestguess[jj] << ' ';
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

		// gsl_odeiv2_driver_reset(d_for);
    	// gsl_odeiv2_driver_reset(d_back);

    } // === (3) Ende

    datei.close();

    // Differentialgleichungen befreien:
    gsl_odeiv2_driver_free (d_for);
    gsl_odeiv2_driver_free (d_back);
    gsl_odeiv2_evolve_free (e_for);
    gsl_odeiv2_control_free (c_for);
    gsl_odeiv2_step_free (s_for);

    gsl_odeiv2_evolve_free (e_back);
    gsl_odeiv2_control_free (c_back);
    gsl_odeiv2_step_free (s_back);

    delete[] x_null;
    // delete[] x_check;
    delete[] x_null_s;
    delete[] x_T;
    delete[] nu;

	x_null = 0;
    // x_check = 0;
    x_null_s = 0;
    x_T = 0;
    nu = 0;

    fcloseall();
}










// ===============================================
// ===============================================
// ========== HILFSROUTINEN ~~~~~~~~~~~~~~~~~~~~~~


// ========== zu (4) ---------------
double stochastic_cool_and_easy(const unsigned long dimVar, std::mt19937 &engine, 
	std::uniform_real_distribution<double> &distribution_uniform, std::normal_distribution<double> &distribution_normal,
	double * x_null, const bool * perturb_null, const int n_ini) {
	

	double x_help, d, r_min_squared, r_max_squared, r_min, r_max, test_laenge, n_dim;

	n_dim = std::min(1.0 + n_ini/10.0, double(dimVar));
	n_dim = std::min(n_dim,100.0);
	// std::cout << n_dim << std::endl;

	r_min = 2.0;
	r_max = 12.0;
	r_min_squared = pow(r_min,n_dim);
	r_max_squared = pow(r_max,n_dim);
	
	d = pow( distribution_uniform(engine) * (r_max_squared - r_min_squared) + r_min_squared, 1.0/n_dim ); 
	d = floor(d*1000.0)/1000.0;

	d = 4.0 + n_ini/100.0;


	test_laenge = 0.0;
	for (int ii=0; ii<dimVar; ii++) {
		if (perturb_null[ii]) {
			x_help = (distribution_normal(engine));
			x_help = pow(x_help,5.0);
			x_null[ii] = x_help;
			test_laenge += x_help*x_help;
		} else {
			x_null[ii] = 0.0;
		}
	}


	test_laenge = sqrt(test_laenge);
	for (int ii=0; ii<dimVar; ii++) {
		x_null[ii] = x_null[ii]*d/test_laenge;
	}

	return d;
}



double stochastic_ini_shrinkingsphere(void * params, int (*func)(double, const double *, double *, void *), gsl_odeiv2_driver * d_for,
		            std::mt19937 &engine, std::uniform_real_distribution<double> &distribution_uniform, 
					std::normal_distribution<double> &distribution,
					double * x_null, double * x_check, const bool * perturb_null, unsigned long dimVar) {


	// Hilfsparameter
	double d, laenge, dStep, d_best;
	bool is_in;
	double x_help;
    laenge = 0.0;

	double r_min, r_max;

	double * x_null_best;

	x_null_best = new double[dimVar];
	for (size_t ii=0; ii<dimVar; ii++) {
		x_null_best[ii] = 0.0;
	}

	d_best = 20.0;

	r_max = 20.0;
	r_min = 0.75*r_max;
	d = r_min + (r_max-r_min)*distribution_uniform(engine);
	d = floor(d*1000.0)/1000.0;
	dStep = 0.1;

	std::cout << "Start ini at d = " << d << std::endl;

	
    is_in = true;
    // ===== (A) Suche nach Wert ausserhalb des Einzuggebietes:
    for (int jj_ini=0; jj_ini<50; jj_ini++) {

    	// ===== (D) d erhoehen?
		d = r_min + (r_max-r_min)*distribution_uniform(engine);
		d = floor(d*1000.0)/1000.0;		

        // ===== (D) Ende

    	// ===== (B) Zufaellige Anfangsbedingungen:
        // ----- zufaellige Richtung bestimmen
        laenge = 0.0;
	    for (size_t ii = 0; ii<dimVar; ii++) {
		    if (perturb_null[ii]) {
				x_help = distribution(engine);
		        x_null[ii] = x_help;
		        laenge += x_help*x_help;
		    }
	    } // ---
	    laenge = sqrt(laenge);

	    // ----- Richtungsvektor (Anfangsbedingungen) auf d normieren
	    for (size_t ii = 0; ii<dimVar; ii++) {
	    	if (perturb_null[ii]) {
		        x_null[ii] *= (d/laenge);
	    	}
			x_check[ii] = x_null[ii];
	    } // ---
	    // ===== (B) Ende



	    // ===== (C) Noch im Einzugsgebiet:
	    // ----- Integration:
	    is_in = inside_basin_second_new(&params, func, d_for, x_check, dimVar, 0);

        // ===== (C) Ende

		if (!is_in) {
			d_best = d;
			r_max = d_best;
			r_min = 0.75*r_max;
			for (size_t ii=0; ii<dimVar; ii++) {
				x_null_best[ii] = x_null[ii];
			}

			// std::cout << jj_ini << " ,  d_best = " << d_best << std::endl;
		}


		while (false) { // (!is_in)
			laenge = d;
			d -= dStep;
			for (size_t ii = 0; ii<dimVar; ii++) {
	    		if (perturb_null[ii]) {
		        	x_null[ii] *= (d/laenge);
	    		}
				x_check[ii] = x_null[ii];
			}

			is_in = inside_basin_second_new(&params, func, d_for, x_check, dimVar, 0);

			if (!is_in) {
				d_best = d;
				r_max = d_best;
				r_min = 0.75*r_max;
				for (size_t ii=0; ii<dimVar; ii++) {
					x_null_best[ii] = x_null[ii];
				}

				// std::cout << "it does ,  d_best = " << d_best << std::endl;
			}

		}

    } // === (A) Ende


    for (size_t ii = 0; ii<dimVar; ii++) {
		x_null[ii] = x_null_best[ii];	
    }

	delete [] x_null_best;


	if (d_best>=20.0) {
		std::cout << "Initialization did not work out." << std::endl;
		exit(23);
	}

	return d_best;
}

// ========== (4) Ende -------------



// ========== zu (5) ---------------

bool inside_basin_second_new(void * params, int (*func)(double, const double *, double *, void *), gsl_odeiv2_driver * d_for,
		          double * x_start ,unsigned long dimVar, signed long count) {


	// ----- Hilfsparameter:
	double t, max_xT, test_xT;
	t = 0.0;
	max_xT = 0.0;
	test_xT = 0.0;


	// ----- Integration:
	int status = gsl_odeiv2_driver_apply (d_for, &t, 50.0, x_start);

	if (status != GSL_SUCCESS) {
	    printf ("error, return value=%d\n", status);
		exit(23);
	}

	// ----- Sind die Anfangsbedingungen schon im Einzugsgebiet:
	max_xT = 0.0;
	for (size_t jj = dimVar/2; jj<dimVar; jj++) {
	    test_xT = x_start[jj]*x_start[jj];
	    if (test_xT > max_xT) {
	    	max_xT = test_xT;
	    }
	}

	max_xT = sqrt(max_xT);


	if (max_xT < 0.1) {
		count+=7;
	} else if (max_xT > 5.0) {
		count-=6;
	} else {
		count--;
	}


	if (count >= 40) {
		gsl_odeiv2_driver_reset(d_for);
		return true;
	} else if (count <= -40) {
		gsl_odeiv2_driver_reset(d_for);
		return false;
	} else {
		return inside_basin_second_new(&params, func, d_for, x_start, dimVar, count);
	}



}
// ========== (5) Ende -------------




double det_distance(void * params, int (*func)(double, const double *, double *, void *), gsl_odeiv2_system * sys_for,
        gsl_odeiv2_step * s_for, gsl_odeiv2_control * c_for, gsl_odeiv2_evolve * e_for,
        unsigned long dimVar, double T, const bool * var_norm,
		double eps_new, double max_eps, double d, double * x_null, double * nu,
		double max_eps_h1, double max_eps_h2, const bool * perturb_null, double * x_run) {

	double dist_new, epsilon, nicht_wurz, wurz, lambda;

	epsilon = eps_new*max_eps;

	nicht_wurz = max_eps_h2/(d*d) - 1.0/epsilon;
	wurz = 1.0/(epsilon*epsilon) + max_eps_h2*max_eps_h2/(d*d*d*d) - max_eps_h1/(d*d);
	lambda = nicht_wurz + std::sqrt(wurz);

	for (size_t ii=0; ii<dimVar; ii++) {
		if (perturb_null[ii]) {
			x_run[ii] = x_null[ii] + epsilon*(lambda*x_null[ii]-nu[ii]);
		} else {
		    x_run[ii] = x_null[ii];
		}
	}

	dist_new = integrate_forward_second_new(params, func, sys_for,
			s_for, c_for, e_for, x_run, dimVar, T, var_norm);


	return dist_new;
}




// ========== zu (7) ---------------
double integrate_forward_second_new(void * params, int (*func)(double, const double *, double *, void *), gsl_odeiv2_system * sys_for,
		                 gsl_odeiv2_step * s_for, gsl_odeiv2_control * c_for, gsl_odeiv2_evolve * e_for,
		                 double * x_T, unsigned long dimVar, double T, const bool * var_norm) {

	double t = 0.0;
	// double h = 1e-3;
	double dist_current = 0.0;
	double dist_help = 0.0;
	double dt;

	int status;
	double h = 0.01;


	// =========== ??? ==========
	while (t < T-1e-10) {
		dt = -t;

		status = GSL_FAILURE;
		while (status != GSL_SUCCESS) {

		    status = gsl_odeiv2_evolve_apply_fixed_step(e_for, c_for, s_for, sys_for, &t, h, x_T);

		    if (status != GSL_SUCCESS) {
		        h /= 2.0;
		    }

		}

        dt += t;
        dist_help = 0.0;

        for(size_t ii = 0; ii < dimVar; ii++) {
        	if (var_norm[ii]) {
        		dist_help += x_T[ii]*x_T[ii];
        	}
        }

        dist_current += 1.0*dt*dist_help;

	}

    gsl_odeiv2_evolve_reset(e_for);
	gsl_odeiv2_step_reset(s_for);

	return dist_current/T;
}
// ========== (7) Ende -------------



// ========== zu (8) ---------------
void integrate_backward_second_new(void * params, int (*func_back)(double, const double *, double *, void *),
		                gsl_odeiv2_system * sys_back, gsl_odeiv2_step * s_back, gsl_odeiv2_control * c_back,
						gsl_odeiv2_evolve * e_back, double * x_T, double * nu, unsigned long dimVar, double T, 
						double * y) {

	double t = 0.0;
	// double * y;

	int status;

	double h = 0.01;

	// ----- Vorbereitung Zustandsvariablen:
	// y = new double[dimVar*2];

	for (size_t ii = 0; ii<dimVar; ii++) {
		y[ii] = x_T[ii];
	}
	for (size_t ii=dimVar; ii<dimVar*2; ii++) {
		y[ii] = nu[ii-dimVar];
	}

    // ----- Eigentliche Integration:
	while (t < T-1e-10) {

		status = GSL_FAILURE;

		while (status != GSL_SUCCESS) {
		    status = gsl_odeiv2_evolve_apply_fixed_step(e_back, c_back, s_back, sys_back, &t, h, y);

            if (status != GSL_SUCCESS) {
	            h /= 2.0;
	        }
		}
	}

    // ----- Auslesen von nu(0):
	for (size_t ii = dimVar; ii<dimVar*2; ii++) {
	    nu[ii-dimVar] = y[ii];
	    x_T[ii-dimVar] = y[ii-dimVar];
	}

	gsl_odeiv2_evolve_reset(e_back);
	gsl_odeiv2_step_reset(s_back);

	// delete[] y;
	// y = 0;

}
// ========== (8) Ende -------------


