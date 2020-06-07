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



#include "mifas_smooth.h"
#include "mifas_new.h"
#include "myDifferentials.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
//#include <gsl/gsl_math.h>
#include <random>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>


// =========== HAUPTROUTINE ~~~~~~~~~~~~~~~~~~~~~
void determine_mifas_smooth(void * params, int (*func)(double, const double *, double *, void *),
		                         int (*func_back)(double, const double *, double *, void *), unsigned long dimVar,
		                         const bool * perturb_null, const bool * var_norm,
								 std::string name_results, unsigned long complete_runs) {


	// ----- ========== DE-KLARA-TION ========== -----

	// ===== (1) Parameter:
    // mit * gekennzeichnete Variablen sollen langfristig von aussen vorgegeben werden
    double epsilon;                // Schrittweite der Anpassung
    double eps_min;                // * Default- und minimaler Wert fuer eps_stern
    double max_eps;                // maximal moeglicher absoluter Wert fuer epsilon
    double max_eps_h1, max_eps_h2; // Hilfsgroessen zur Bestimmung von max_eps
    double lambda;                 // Groesse fuer Verhaeltnis von x_null zu nu
    double wurz, nicht_wurz;       // Hilfsgroessen zur Berechnung von lambda

    double eps_a, eps_b, eps_c;
    double dist_b, dist_c, dist_current;
    double eps_b_ini, eps_c_ini, eps_dist_ac;
    double max_eps_fak;
    bool stay_tuned;
    double eps_help, dist_help;
    double T, T_ini;

    double d;             // * Aktueller Radius (Distanz der Anfangsbedingungen zum FP)
    double d_min;         // bislang minimaler Wert von d
    double dStep;         // Schrittweite zur Herabsetzung von d
    double test_laenge;   // Hilfsgroesse fuer die Normierung von x_null

    unsigned long jj;     // Laufvariable


	// === NEW INI === 
	bool OLD_INI = true;
	bool finIni = false;	
	int n_ini;
	

	int count_smooth, count_smooth_crit, count_smooth_crit_a, count_smooth_crit_b;
	double d_s, dStep_ini, T_plus;
	bool inside_smooth;
	int small_buffer, small_buffer_default;

	inside_smooth = false;
	d_s = 10.0;

	count_smooth_crit_a = 30;
	count_smooth_crit_b = 5;

	std::ofstream datei;


	

	// !!! ???  !!! TRY THESE !!! ??? !!!
	T_plus = 0.0;
	dStep_ini = 0.01;
	small_buffer_default = 3;





	// ----- NEU:
	int count_in_basin = 0;
	int count_in_basin_crit = 5;

	


    // ----- Arrays:
    double *x_null, *x_T, *x_null_s, *nu, *x_star, *y;

    x_null   = new double[dimVar];        // dynamische Anfangsbedingungen
    x_T      = new double[dimVar];        // Speicheradresse fuer die Simulation

    x_null_s = new double[dimVar];        // Zwischenspeicher

    nu       = new double[dimVar];        // Abgeleitete Zustandsgroessen

    x_star = new double[dimVar];		  // Fixpunkt
    // x_check = new double[dimVar];		  // zum Testen

	y = new double[dimVar*2];			  // Rueckwaerts


    // ----- Dateinamen in char-Arrays umwandeln:
    const char * fname = name_results.c_str();


    // ----- ========== DEFINITION ========== -----
    eps_b_ini = 0.01;
    eps_c_ini = 0.02;


    max_eps_fak = 0.5;
    eps_dist_ac = 0.5*1.0e-2;
    eps_min = 1.0e-3; 

    d     = 1.0;
    d_min = 100.0;
	T_ini = 10.0;


    // Neu: Fixpunkt auslesen:
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



    // ----- ========== INITIALISIERUNG ========== -----

	// ===== (2a) Initialiserung der Diff.-Gleichungen:
    // Forwards (a):
	gsl_odeiv2_system sys_for = {func, 0, dimVar, params};
    gsl_odeiv2_driver * d_for = gsl_odeiv2_driver_alloc_y_new(&sys_for, gsl_odeiv2_step_rkf45, 1e-6, 1e-8, 1e-8);
    gsl_odeiv2_driver_set_hmax (d_for, 0.01);

    // Forwards (b):
    const gsl_odeiv2_step_type * st_for = gsl_odeiv2_step_rkf45;
    gsl_odeiv2_step * s_for = gsl_odeiv2_step_alloc(st_for, dimVar);
    gsl_odeiv2_control * c_for = gsl_odeiv2_control_y_new(1e-3, 1e-3);
    gsl_odeiv2_evolve * e_for = gsl_odeiv2_evolve_alloc(dimVar);

    // Backwards:
	gsl_odeiv2_system sys_back = {func_back, 0, 2*dimVar, params};
    gsl_odeiv2_driver * d_back = gsl_odeiv2_driver_alloc_y_new(&sys_back, gsl_odeiv2_step_rkf45, 1e-6, 1e-8, 1e-8);
    gsl_odeiv2_driver_set_hmax (d_back, 0.01);

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

    // ===== (3) Schleife ueber mehrere Durchlaeufe:
    for (size_t ii = 0; ii < complete_runs; ii++) {

		dStep = dStep_ini;

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
    	}

		if (OLD_INI) { 

			
        	d = stochastic_ini_shrinkingsphere(&params, func, d_for, engine, distribution_uniform, distribution, x_null, x_T, 
				perturb_null, dimVar);
				
			for (jj=0; jj<dimVar; jj++) {
				x_null_s[jj] = x_null[jj];
			}
			d_s = d;


			finIni = true;


		} else {

			d = stochastic_cool_and_easy(dimVar, engine, distribution_uniform, distribution, 
				x_null, perturb_null, n_ini);

			// === NEW INI ===
			finIni = false;

		}

		std::cout << d << std::endl;


        // ===== (4) Ende (siehe Funktion unten)



		dist_b = 1.0e5;
		dist_current = 1.0e-5;


		// ===== (6) Eigentliche Methode nach Kerswell

		count_smooth = 0;
		count_smooth_crit = count_smooth_crit_a;
		count_in_basin_crit = 5;
		count_in_basin = 0;
		small_buffer = 100;

		while ( dStep > 1e-5) {


			// !!!!!!!!!!!!!!!!!!!!!!


			test_laenge = 0.0;
			for (jj=0; jj<dimVar; jj++) {
				if (perturb_null[jj]) {
					test_laenge += x_null[jj]*x_null[jj];
				}
			}
			test_laenge = std::sqrt(test_laenge);


			// Initialisirung von x_T und x_null auf die korrekte neue Laenge bringen:
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



			// ===== (7) Erste Vorwaerts-Integration:
			dist_b = integrate_forward_second_new(params, func, &sys_for, s_for, c_for, e_for, x_T, dimVar, T, var_norm);
	        dist_current = dist_b;


			// VARIABEL: nu muss nicht zwingend auf 0 gesetzt werden (je nach Norm).
			for (jj = 0; jj<dimVar; jj++) {
				nu[jj] = 0.0;
			}


			// ===== (8) Rueckwaerts-Integration:
			integrate_backward_second_new(params, func_back, &sys_back, s_back, c_back, e_back, x_T, nu, dimVar, T, y);


			// ===== (9) Anpassung Schrittweite:
			// ----- nu anpassen und maximalen Wert fuer epsilon bestimmen:
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


			// ===== (10) Schleife Schrittweitenanpassung:


            // DINGDINGDINGDINGDINGDING
			while (dist_b<=dist_current and eps_b > eps_min) {

				dist_b = det_distance(params, func, &sys_for,
					s_for, c_for, e_for, dimVar, T, var_norm,
					eps_b,max_eps,d,x_null,nu,
					max_eps_h1,max_eps_h2,perturb_null,x_T);

				// std::cout << "See: " << dist_b << std::endl;

				if (dist_b <= dist_current) {
					eps_c = eps_b;
					dist_c = dist_b;
					eps_b *= 0.5;
				}
			}

			if (dist_c<0.0 and eps_c < 1.0) {
				dist_c = dist_b+1.0;
			}


			// maximum is already reached ??
			if (eps_b > eps_min) {

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


	            // Setze neues x_null nur einmal am Ende
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

			}




			// =================================================================================
			// =================================================================================
			// =================================================================================

			if (small_buffer <= 0) {

				count_smooth++;

				small_buffer = small_buffer_default;

				if (count_smooth > count_smooth_crit || eps_b > 0.9) {
					for (jj=0; jj<dimVar; jj++) {
						x_T[jj] = x_null[jj];				
					}

					inside_smooth = inside_basin_second_new(&params, func, d_for, x_T, dimVar, 0);
					// std::cout << "inside -> " << inside_smooth << std::endl;

					if (!inside_smooth) {
						for (jj=0; jj<dimVar; jj++) {
							x_null_s[jj] = x_null[jj];
						}
						d_s = d;
						count_in_basin = 0;
					} else {
						count_in_basin++;
						if (count_in_basin > count_in_basin_crit) {
							for (jj=0; jj<dimVar; jj++) {
								x_null[jj] = x_null_s[jj];
							}

							T += T_plus;

							if (count_smooth_crit == count_smooth_crit_a) {
								count_smooth_crit = count_smooth_crit_b;
								count_in_basin_crit = 5;
							} else if (count_smooth_crit == count_smooth_crit_b) { 
								count_smooth_crit = 0;
								if (dStep < 0.0025) {
									count_in_basin_crit = 15;
								} else {
									count_in_basin_crit = 5;
								}
							} else { 
								count_smooth_crit = count_smooth_crit_b;
								count_in_basin_crit = 5;
								dStep *= 0.1;
							}

							d = d_s + dStep;
							count_in_basin = 0;
						}

						std::cout << "went in at  d = " << d << " ;   cibc = " << count_in_basin_crit << "  " <<
							count_in_basin << std::endl;
					}

					small_buffer = 10 + small_buffer_default;
					count_smooth = 0;
				} 

				d -= dStep;
					
			} else {

				small_buffer--;

			}

			// std::cout << "d = " << d << " ;   dist_b = " << dist_b << " ;   eps_b = " << eps_b << std::endl;

			// =================================================================================
			// =================================================================================
			// =================================================================================




		} // === (6) Ende





		std::cout << "local MiFaS = " << d_s << std::endl;






		// === NEW INI ===
		// Hier muss der Abschnitt zur Speicherung der Daten hin:
		if (finIni) {
        	datei << T_ini << ' ' << d_s-dStep*10.0 << ' ' << d_s << ' ';
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

