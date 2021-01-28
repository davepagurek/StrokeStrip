/*--
    Arc.h  

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

#ifndef CORNUCOPIA_ARC_H_INCLUDED
#define CORNUCOPIA_ARC_H_INCLUDED

#include "defs.h"
#include "CurvePrimitive.h"

NAMESPACE_Cornu

CORNU_SMART_FORW_DECL(Arc);

class Arc : public CurvePrimitive
{
public:
    Arc() {} //uninitialized
    Arc(const Vec &start, double startAngle, double length, double curvature);
    Arc(const Vec &start, const Vec &mid, const Vec &end);

    //overrides
    void eval(double s, Vec *pos, Vec *der = NULL, Vec *der2 = NULL) const;

    double project(const Vec &point) const;

    double angle(double s) const { return _startAngle() + s * _params[CURVATURE]; }
    double curvature(double s) const { return _params[CURVATURE]; }

    double endAngle() const { return _startAngle() + _angleDiff; }
    double startCurvature() const { return _params[CURVATURE]; }
    double endCurvature() const { return _params[CURVATURE]; }

    PrimitiveType getType() const { return ARC; }

    void trim(double sFrom, double sTo);
    void flip();
    CurvePrimitivePtr clone() const { ArcPtr out = new Arc(); out->setParams(_params); return out; }
    void derivativeAt(double s, ParamDer &out, ParamDer &outTan) const;
    void derivativeAtEnd(int continuity, EndDer &out) const;

    //arc specific--UNDEFINED if arc is flat
    Vec center() const { return _center; }
    double radius() const { return _radius; }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
protected:
    //override
    void _paramsChanged();
    bool isValidImpl() const;

private:
    Vec _tangent; //at start
    Vec _center; //if arc is not flat
    double _radius; //if arc is not flat, 1 / curvature
    double _angleDiff; //length * curvature
    bool _flat; //if true, arc is almost flat and we should use an approximation
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_ARC_H_INCLUDED
