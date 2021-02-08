/*--
    Bezier.cpp  

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

#include "Bezier.h"

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

CubicBezier CubicBezier::hermite(const Vec &start, const Vec &end, const Vec &startDer, const Vec &endDer)
{
    CubicBezier out;
    out._pts[0] = start;
    out._pts[3] = end;
    out._pts[1] = start + startDer / 3.;
    out._pts[2] = end - endDer / 3.;

    return out;
}

void CubicBezier::eval(double t, Vec *pos, Vec *der, Vec *der2) const
{
    if(pos)
        *pos = CUBE(1 - t) * _pts[0] + CUBE(t) * _pts[3] + 3. * t * (1 - t) * ((1 - t) * _pts[1] + t * _pts[2]);
    if(der)
        *der = 3 * (-SQR(1 - t) * _pts[0] + SQR(t) * _pts[3] + (3 * SQR(t) - 4 * t + 1) * _pts[1] + (2 * t - 3 * SQR(t)) * _pts[2]);
    if(der2)
        *der2 = 6 * ((1 - t) * _pts[0] + t * _pts[3] + (3 * t - 2) * _pts[1] + (1 - 3 * t) * _pts[2]);
}

void BezierSpline::eval(double t, Vec *pos, Vec *der, Vec *der2) const
{
    int idx = int(t);
    idx = max(0, min((int)_primitives.size() - 1, idx));
    _primitives[idx].eval(t - idx, pos, der, der2);
}

END_NAMESPACE_Cornu


