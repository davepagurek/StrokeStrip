/*--
    PrimitiveSequence.h  

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

#ifndef CORNUCOPIA_PRIMITIVESEQUENCE_H_INCLUDED
#define CORNUCOPIA_PRIMITIVESEQUENCE_H_INCLUDED

#include "defs.h"
#include "CurvePrimitive.h"
#include "VectorC.h"

NAMESPACE_Cornu

CORNU_SMART_FORW_DECL(PrimitiveSequence);
CORNU_SMART_FORW_DECL(BezierSpline);

class PrimitiveSequence : public Curve
{
public:
    PrimitiveSequence(const VectorC<CurvePrimitiveConstPtr> &primitives);

    //overrides
    double length() const { return _lengths.back(); }
    bool isClosed() const { return _primitives.circular() == CIRCULAR; }

    void eval(double s, Vec *pos, Vec *der = NULL, Vec *der2 = NULL) const;

    double project(const Vec &point) const;

    //utility functions
    int paramToIdx(double param, double *outParam = NULL) const;
    bool isParamValid(double param) const { return isClosed() || (param >= 0 && param <= _lengths.back()); }

    //Returns the curve trimmed between the input arguments, if possible, but at most
    //the length of the original curve.  Arguments can be negative and to can be smaller
    //than from for a closed curve.
    PrimitiveSequencePtr trimmed(double from, double to) const;
    PrimitiveSequencePtr flipped() const;

    const VectorC<CurvePrimitiveConstPtr> &primitives() const { return _primitives; }

    BezierSplinePtr toBezierSpline(double tolerance) const;

private:
    VectorC<CurvePrimitiveConstPtr> _primitives;
    //lengths[x] = \sum_{i=1}^{i=x} _primitives[i-1]->length(), i.e., length up to the start of the primitive at x
    std::vector<double> _lengths; 
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_PRIMITIVESEQUENCE_H_INCLUDED
