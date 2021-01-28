/*--
    Fresnel.h  

    This file is part of the Cornucopia curve sketching library.
    Copyright (C) 2010 Ilya Baran (baran37@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CORNUCOPIA_FRESNEL_H_INCLUDED
#define CORNUCOPIA_FRESNEL_H_INCLUDED

#include "defs.h"
#include <Eigen/Core>

NAMESPACE_Cornu

//almost full double-precision accuracy, using rational approximations
void fresnel(double xxa, double *ssa, double *cca);
void fresnel(const Eigen::VectorXd &t, Eigen::VectorXd *s, Eigen::VectorXd *c);

//roughly single-precision accuracy, using polynomial approximations
void fresnelApprox(double xxa, double *ssa, double *cca);
void fresnelApprox(const Eigen::VectorXd &t, Eigen::VectorXd *s, Eigen::VectorXd *c); //sse vectorized

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_FRESNEL_H_INCLUDED
