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



#ifndef MIFAS_NEW_H_
#define MIFAS_NEW_H_

#include <random>
#include <string>
#include <gsl/gsl_odeiv2.h>

void determine_mifas_new(void *, int (*)(double, const double *, double *, void *),
		          int (*)(double, const double *, double *, void *), unsigned long, const bool *, const bool *,
				  std::string, unsigned long);

double stochastic_cool_and_easy(const unsigned long, std::mt19937 &, 
	std::uniform_real_distribution<double> &, std::normal_distribution<double> &,
	double *, const bool *, const int);

double stochastic_ini_shrinkingsphere(void *, int (*)(double, const double *, double *, void *), gsl_odeiv2_driver *,
		            std::mt19937 &,std::uniform_real_distribution<double> &, std::normal_distribution<double> &, 
					double*, double*, const bool*, unsigned long);

bool inside_basin_second_new(void *, int (*)(double, const double *, double *, void *),
		          gsl_odeiv2_driver *, double *, unsigned long, signed long);

double integrate_forward_second_new(void *, int (*)(double, const double *, double *, void *), gsl_odeiv2_system *,
		                 gsl_odeiv2_step *, gsl_odeiv2_control *, gsl_odeiv2_evolve *,
		                 double *, unsigned long, double, const bool *);

void integrate_backward_second_new(void *, int (*)(double, const double *, double *, void *),
		                  gsl_odeiv2_system *, gsl_odeiv2_step *, gsl_odeiv2_control *,
						  gsl_odeiv2_evolve *, double *, double *, unsigned long, double, double *);

double det_distance(void *, int (*func)(double, const double *, double *, void *), gsl_odeiv2_system *,
        gsl_odeiv2_step *, gsl_odeiv2_control *, gsl_odeiv2_evolve *,
        unsigned long, double, const bool *,
		double, double, double, double *, double *,
		double, double, const bool *, double *);

#endif /* MIFAS_NEW_H_ */
