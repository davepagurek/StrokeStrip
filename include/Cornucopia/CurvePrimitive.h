/*--
    CurvePrimitive.h  

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

#ifndef CORNUCOPIA_CURVEPRIMITIVE_H_INCLUDED
#define CORNUCOPIA_CURVEPRIMITIVE_H_INCLUDED

#include "defs.h"
#include "Curve.h"

NAMESPACE_Cornu

CORNU_SMART_FORW_DECL(CurvePrimitive);

/*
    A CurvePrimitive is a line, arc, or a clothoid, defined by 4-6 parameters.
*/
class CurvePrimitive : public Curve
{
public:
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1, Eigen::AutoAlign, 6, 1> ParamVec;
    typedef Eigen::Matrix<double, 2, Eigen::Dynamic, Eigen::AutoAlign, 2, 6> ParamDer; //derivative of x and y w.r.t. parameters
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign, 4, 6> EndDer; //derivative of x, y, angle, curvature w.r.t. parameters

    enum PrimitiveType
    {
        LINE = 0,
        ARC,
        CLOTHOID
    };

    enum Param
    {
        X = 0,
        Y,
        ANGLE,
        LENGTH,
        CURVATURE,
        DCURVATURE
    };

    int numParams() const { return 4 + getType(); } //line has 4 params, etc.

    //interface functions
    virtual PrimitiveType getType() const = 0;

    virtual void trim(double sFrom, double sTo) = 0; //this should work even for sFrom < 0 and sTo > length()
    virtual void flip() = 0;
    virtual CurvePrimitivePtr clone() const = 0;

    //derivative (of curve and its tangent vector) with respect to paramters
    virtual void derivativeAt(double s, ParamDer &out, ParamDer &outTan) const = 0;
    virtual void derivativeAtEnd(int continuity, EndDer &out) const = 0; //continuity: 0 = position, 1 = +angle, 2 = +curvature

    //for clothoids, converts derivative w.r.t. dcurvature into der w.r.t. end curvature
    virtual void toEndCurvatureDerivative(Eigen::MatrixXd &) const {} 

    void setParams(const ParamVec &params) { _params = params; _paramsChanged(); }
    const ParamVec &params() const { return _params; }

    //overrides
    double length() const { return _length(); }
    Vec startPos() const { return _startPos(); }
    double startAngle() const { return _startAngle(); }
    bool isValid() const; //calls isValidImpl internally

    //utility functions
    CurvePrimitivePtr flipped() const { CurvePrimitivePtr out = clone(); out->flip(); return out; }
    CurvePrimitivePtr trimmed(double sFrom, double sTo) const { CurvePrimitivePtr out = clone(); out->trim(sFrom, sTo); return out; }

protected:
    //non-virtual inline functions -- use them in derived classes
    double _length() const { return _params[LENGTH]; }
    Vec _startPos() const { return _params.head<2>(); }
    double _startAngle() const { return _params[ANGLE]; }

    virtual void _paramsChanged() = 0;
    virtual bool isValidImpl() const = 0;

    ParamVec _params;
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_CURVEPRIMITIVE_H_INCLUDED
