// Copyright (C) 2020  Lukas Halekotte <lukas.halekotte@uol.de>
//
// This file is part of the program 'mifas_pp': you can redistribute it 
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


#ifndef MIFAS_H_
#define MIFAS_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

#include <random>


int determine_mifas(Beegrid *, double, const bool *, std::string, unsigned long);

double stochastic_ini(gsl_odeiv2_system *, gsl_odeiv2_step *, gsl_odeiv2_control *, gsl_odeiv2_evolve *, 
	const unsigned long, const double, std::mt19937 &, std::uniform_real_distribution<double> &, 
	std::normal_distribution<double> &,	double *, const bool *, double *, const double *, const int);

double stochastic_ini_B(gsl_odeiv2_system *, gsl_odeiv2_step *, gsl_odeiv2_control *, gsl_odeiv2_evolve *, 
	const unsigned long, const double, std::mt19937 &, std::uniform_real_distribution<double> &, 
	std::normal_distribution<double> &,	double *, const bool *, double *, const double *, const int);


double unstochastic_ini(gsl_odeiv2_system *,gsl_odeiv2_step *, gsl_odeiv2_control *, gsl_odeiv2_evolve *, 
	const unsigned long, const double, double *, double *, const double *, const int);
	

double det_distance(void *, gsl_odeiv2_system *, gsl_odeiv2_step *, gsl_odeiv2_control *, gsl_odeiv2_evolve *,
	unsigned long, double, double, double, double, double, double *, double *, double, double, const bool *, 
	double *, double *);

double integrate_forward_backward(void *, 
	gsl_odeiv2_system *, gsl_odeiv2_system *, gsl_odeiv2_step *, gsl_odeiv2_control *, 
	gsl_odeiv2_evolve *, double *, double *, const double *, double *, 
	const unsigned long, const double, const double);

double integrate_forward(void *, 
	gsl_odeiv2_system *, gsl_odeiv2_step *, gsl_odeiv2_control *, gsl_odeiv2_evolve *,
	double *, const double *, const unsigned long, const double, const double);


bool determineStar(struct paraSimple *, long unsigned int, double *);
bool isIn(double *, const long unsigned int, const double *, gsl_odeiv2_driver *);
bool isIn(double *, const long unsigned int, const double *, gsl_odeiv2_system *,
	gsl_odeiv2_step *, gsl_odeiv2_control *, gsl_odeiv2_evolve *, const double);

#endif
