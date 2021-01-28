/*--
    TwoCurveCombine.h  

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

#ifndef CORNUCOPIA_TWOCURVECOMBINE_H_INCLUDED
#define CORNUCOPIA_TWOCURVECOMBINE_H_INCLUDED

#include "defs.h"
#include "smart_ptr.h"

NAMESPACE_Cornu

class Fitter;
CORNU_SMART_FORW_DECL(CurvePrimitive);

struct Combination
{
    CurvePrimitivePtr c1;
    CurvePrimitivePtr c2;
    double err1;
    double err2;
};

Combination twoCurveCombine(int p1, int p2, int continuity, const Fitter &fitter);

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_TWOCURVECOMBINE_H_INCLUDED
