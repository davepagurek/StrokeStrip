/*--
    defs.h  

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

#ifndef CORNUCOPIA_DEFS_H_INCLUDED
#define CORNUCOPIA_DEFS_H_INCLUDED

#ifdef _WIN32
#pragma warning(disable:4251)
#pragma warning(disable:4275)
#pragma warning(disable:4996)
#endif

//The macros are here so there's no base indentation in IDE's
#define NAMESPACE_Cornu namespace Cornu {
#define END_NAMESPACE_Cornu }

namespace std {}
namespace Eigen {}

#include "Debugging.h"

namespace Cornu
{
    static const double PI =         3.1415926535897932385;
    static const double TWOPI =      6.2831853071795864769;
    static const double HALFPI =     1.5707963267948966192;

    template<typename T> T SQR(const T &in) { return in * in; }
    template<typename T> T CUBE(const T &in) { return in * in * in; }
}

#endif //CORNUCOPIA_DEFS_H_INCLUDED
