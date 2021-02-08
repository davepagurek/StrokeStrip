/*--
    Clothoid.h  

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

#ifndef CORNUCOPIA_CLOTHOID_H_INCLUDED
#define CORNUCOPIA_CLOTHOID_H_INCLUDED

#include "defs.h"
#include "CurvePrimitive.h"

NAMESPACE_Cornu

CORNU_SMART_FORW_DECL(Clothoid);

class Clothoid : public CurvePrimitive
{
public:
    Clothoid() {} //uninitialized
    Clothoid(const Vec &start, double startAngle, double length, double curvature, double endCurvature);

    //overrides
    void eval(double s, Vec *pos, Vec *der = NULL, Vec *der2 = NULL) const;

    double project(const Vec &point) const;

    double angle(double s) const;
    double curvature(double s) const;

    double startCurvature() const { return _params[CURVATURE]; }
    double endCurvature() const { return _params[CURVATURE] + _params[LENGTH] * _params[DCURVATURE]; }

    PrimitiveType getType() const { return CLOTHOID; }

    void trim(double sFrom, double sTo);
    void flip();
    CurvePrimitivePtr clone() const { ClothoidPtr out = new Clothoid(); out->setParams(_params); return out; }
    void derivativeAt(double s, ParamDer &out, ParamDer &outTan) const;
    void derivativeAtEnd(int continuity, EndDer &out) const;

    void toEndCurvatureDerivative(Eigen::MatrixXd &der) const;

    class _ClothoidProjector //internal singleton class
    {
    public:
        virtual double project(const Vec &pt, double from, double to) const = 0;
    };

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
protected:
    //override
    void _paramsChanged();
    bool isValidImpl() const;

private:
    Vec _startShift; //translation component of transformation from canonical clothoid
    Eigen::Matrix2d _mat; //rotation and scale component of transformation from canonical clothoid
    double _t1; //start parameter on the canonical clothoid
    double _tdiff;
    bool _arc;
    bool _flat;

    class _ClothoidProjectorImpl;
    static _ClothoidProjector *_clothoidProjector(); //projects onto a generic clothoid
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_CLOTHOID_H_INCLUDED
