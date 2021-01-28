/*--
    Bezier.h  

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

#ifndef CORNUCOPIA_BEZIER_H_INCLUDED
#define CORNUCOPIA_BEZIER_H_INCLUDED

#include "defs.h"
#include "smart_ptr.h"
#include "VectorC.h"

NAMESPACE_Cornu

class CubicBezier
{
public:
    typedef Eigen::Vector2d Vec;

    static CubicBezier hermite(const Vec &start, const Vec &end, const Vec &startDer, const Vec &endDer);

    //t is between zero and one
    void eval(double t, Vec *pos, Vec *der = NULL, Vec *der2 = NULL) const;

    const Vec &controlPoint(int i) const { return _pts[i]; }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    Vec _pts[4];
};

CORNU_SMART_FORW_DECL(BezierSpline);

class BezierSpline : public smart_base
{
public:
    typedef VectorC<CubicBezier, Eigen::aligned_allocator<CubicBezier> > PrimitiveVector;
    typedef Eigen::Vector2d Vec;

    BezierSpline(const PrimitiveVector &primitives) : _primitives(primitives) {}

    //t should be between 0 and primitives().size()
    void eval(double t, Vec *pos, Vec *der = NULL, Vec *der2 = NULL) const;

    const PrimitiveVector &primitives() const { return _primitives; }

private:
    PrimitiveVector _primitives;
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_BEZIER_H_INCLUDED
