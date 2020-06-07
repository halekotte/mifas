/// Copyright (C) 2020  Lukas Halekotte <lukas.halekotte@uol.de>
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



#ifndef MIFAS_SMOOTH_H_
#define MIFAS_SMOOTH_H_

#include <random>
#include <string>
#include <gsl/gsl_odeiv2.h>

void determine_mifas_smooth(void *, int (*)(double, const double *, double *, void *),
		          int (*)(double, const double *, double *, void *), unsigned long, const bool *, const bool *,
				  std::string, unsigned long);



#endif /* MIFAS_SMOOTH_H_ */
