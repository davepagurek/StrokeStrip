/*--
    Line.cpp  

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

#include "Line.h"

using namespace std;
using namespace Eigen;

NAMESPACE_Cornu

Line::Line(const Vec &p1, const Vec &p2)
{
    _params.resize(numParams());

    _params.head<2>() = p1;
    Vec dir = p2 - p1;
    double len = dir.norm();

    if(fabs(len) < NumTraits<double>::dummy_precision()) //fail gracefully on zero-length line
    {
        _params[LENGTH] = 0;
        _params[ANGLE] = 0;
        _der = Vec(1., 0.);
    }
    else
    {
        _params[LENGTH] = len;
        _params[ANGLE] = AngleUtils::angle(dir);
        _der = dir * (1. / len);
    }
}

bool Line::isValidImpl() const
{
    if(_params[LENGTH] < 0.)
        return false;
    return true;
}

double Line::project(const Vec &point) const
{
    return min(_length(), max(0., _der.dot(point - _startPos())));
}

void Line::eval(double s, Vec *pos, Vec *der, Vec *der2) const
{
    if(pos)
        *pos = _startPos() + s * _der;
    if(der)
        *der = _der;
    if(der2)
        *der2 = Vec::Zero();
}

void Line::trim(double sFrom, double sTo)
{
    Vec newStart = _startPos() + sFrom * _der;
    _params[LENGTH] = sTo - sFrom;
    _params[X] = newStart[0];
    _params[Y] = newStart[1];
}

void Line::flip()
{
    Vec newStart = endPos();
    _params[ANGLE] = PI + _params[ANGLE];
    _der = -_der;
    _params[X] = newStart[0];
    _params[Y] = newStart[1];
}

void Line::derivativeAt(double s, ParamDer &out, ParamDer &outTan) const
{
    outTan = out = ParamDer::Zero(2, 4);
    out(0, X) = 1;
    out(1, Y) = 1;
    out(0, ANGLE) = -s * _der(1);
    out(1, ANGLE) = s * _der(0);
    outTan(0, ANGLE) = -_der(1);
    outTan(1, ANGLE) = _der(0);
}

void Line::derivativeAtEnd(int continuity, EndDer &out) const
{
    out = EndDer::Zero(2 + continuity, 4);
    out(0, X) = 1;
    out(1, Y) = 1;
    out(0, LENGTH) = _der(0);
    out(1, LENGTH) = _der(1);
    out(0, ANGLE) = -_length() * _der(1);
    out(1, ANGLE) = _length() * _der(0);
    if(continuity >= 1)
        out(2, ANGLE) = 1.;
}

END_NAMESPACE_Cornu


