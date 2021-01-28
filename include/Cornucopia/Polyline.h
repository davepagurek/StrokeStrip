/*--
    Polyline.h  

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

#ifndef CORNUCOPIA_POLYLINE_H_INCLUDED
#define CORNUCOPIA_POLYLINE_H_INCLUDED

#include "defs.h"
#include "Curve.h"
#include "VectorC.h"

NAMESPACE_Cornu

CORNU_SMART_FORW_DECL(Polyline);

//A polyline must have at least two points
class Polyline : public Curve
{
public:
    Polyline(const VectorC<Eigen::Vector2d> &pts);

    //overrides
    double length() const { return _lengths.back(); }
    bool isClosed() const { return _pts.circular() == CIRCULAR; }

    void eval(double s, Vec *pos, Vec *der = NULL, Vec *der2 = NULL) const;

    double project(const Vec &point) const;

    //utility functions
    int paramToIdx(double param, double *outParam = NULL) const;
    double idxToParam(int idx) const { return _lengths[idx]; }
    bool isParamValid(double param) const { return isClosed() || (param >= 0 && param <= _lengths.back()); }
    double lengthFromTo(int fromIdx, int toIdx) const; //for closed polylines, toIdx can be smaller than fromIdx

    //Returns the curve trimmed between the input arguments, if possible, but at most
    //the length of the original curve.  Arguments can be negative and to can be smaller
    //than from for a closed polyline.
    PolylinePtr trimmed(double from, double to) const;

    const VectorC<Eigen::Vector2d> &pts() const { return _pts; }

private:
    VectorC<Eigen::Vector2d> _pts;
    //lengths[x] = \sum_{i=1}^{i=x} ||pts[i]-pts[i-1]||, i.e., length up to point x
    //For closed curves, _lengths has a last member which is the total curve length.
    std::vector<double> _lengths; 
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_POLYLINE_H_INCLUDED
