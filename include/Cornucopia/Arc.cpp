/*--
    Arc.cpp  

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

#include "Arc.h"

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

Arc::Arc(const Vec &start, double startAngle, double length, double curvature)
{
    _params.resize(numParams());

    _params.head<2>() = start;
    _params[ANGLE] = AngleUtils::toRange(startAngle);
    _params[LENGTH] = length;
    _params[CURVATURE] = curvature;

    _paramsChanged();
}

void Arc::_paramsChanged()
{
    _tangent = Vec(cos(_startAngle()), sin(_startAngle()));

    _angleDiff = _length() * _params[CURVATURE];
    _flat = fabs(_params[CURVATURE]) < 1e-6;

    if(!_flat)
    {
        Vec toCenter(-_tangent[1], _tangent[0]);
        _radius = 1. / _params[CURVATURE];
        _center = _startPos() + _radius * toCenter;
    }
}

bool Arc::isValidImpl() const
{
    if(_params[LENGTH] < 0.)
        return false;

    if(fabs(_angleDiff) > TWOPI)
        return false;

    return true;
}

void Arc::eval(double s, Vec *pos, Vec *der, Vec *der2) const
{
    double angle = _startAngle() + s * _params[CURVATURE];
    double cosa, sina;
    if(!_flat || der || der2)
    {
        cosa = cos(angle);
        sina = sin(angle);
    }

    if(_flat && pos)
        (*pos) = _startPos() + _tangent * s;
    else if(pos)
        (*pos) = _center + _radius * Vec(sina, -cosa);

    if(der)
        (*der) = Vec(cosa, sina);
    if(der2)
        (*der2) = Vec(-sina, cosa) * _params[CURVATURE];
}

double Arc::project(const Vec &point) const
{
    double t;
    if(_flat)
    {
        t = (point - _startPos()).dot(_tangent);
    }
    else
    {
        double projAngle = atan2(point[1] - _center[1], point[0] - _center[0]);
        //To compute the projection, get the angle difference into the range centered on the midpoint of
        //the arc.
        if(_params[CURVATURE] > 0.)
        {
            double projAngleDiff = HALFPI + projAngle - _startAngle();
            t = AngleUtils::toRange(projAngleDiff, 0.5 * _angleDiff - PI) * _radius;
        }
        else
        {
            double projAngleDiff = HALFPI - projAngle + _startAngle();
            t = -AngleUtils::toRange(projAngleDiff, -0.5 * _angleDiff - PI) * _radius;
        }
    }

    return max(0., min(_length(), t));
}

void Arc::trim(double sFrom, double sTo)
{
    Vec newStart = pos(sFrom);
    _params[ANGLE] = angle(sFrom);
    _params[LENGTH] = sTo - sFrom;
    _params[X] = newStart[0];
    _params[Y] = newStart[1];

    _paramsChanged();
}

void Arc::flip()
{
    Vec newStart = endPos();
    _params[ANGLE] = PI + endAngle();
    _params[CURVATURE] = -_params[CURVATURE];
    _params[X] = newStart[0];
    _params[Y] = newStart[1];

    _paramsChanged();
}

void Arc::derivativeAt(double s, ParamDer &out, ParamDer &outTan) const
{
    outTan = out = ParamDer::Zero(2, 5);
    out(0, X) = 1;
    out(1, Y) = 1;

    Vec diff = pos(s) - _startPos();
    out(0, ANGLE) = -diff[1];
    out(1, ANGLE) = diff[0];

    if(_flat)
    {
        outTan.col(ANGLE) = Vec(-_tangent[1], _tangent[0]);
        out.col(CURVATURE) = (0.5 * s * s) * Vec(-_tangent[1], _tangent[0]);
        outTan.col(CURVATURE) = s * Vec(-_tangent[1], _tangent[0]);
    }
    else
    {
        double angle = _startAngle() + s * _params[CURVATURE];

        double cosa, sina;
        cosa = cos(angle);
        sina = sin(angle);

        outTan.col(ANGLE) = Vec(-sina, cosa);
        out(0, CURVATURE) = (s * cosa + (_tangent[1] - sina) * _radius) * _radius;
        out(1, CURVATURE) = (s * sina - (_tangent[0] - cosa) * _radius) * _radius;
        outTan(0, CURVATURE) = -s * sina;
        outTan(1, CURVATURE) = s * cosa;
    }
}

void Arc::derivativeAtEnd(int continuity, EndDer &out) const
{
    out = EndDer::Zero(2 + continuity, 5);

    ParamDer pDer, dummy;
    derivativeAt(_length(), pDer, dummy);
    Vector2d endTan;
    eval(_length(), NULL, &endTan);

    out.block(0, 0, 2, 5) = pDer;
    out.col(LENGTH).head<2>() = endTan;

    if(continuity >= 1)
    {
        out(2, ANGLE) = 1.;
        out(2, LENGTH) = _params[CURVATURE];
        out(2, CURVATURE) = _params[LENGTH];
    }

    if(continuity == 2)
        out(3, CURVATURE) = 1.;
}

Arc::Arc(const Vec &start, const Vec &mid, const Vec &end)
{
    _params.resize(numParams());

    for(int i = 0; i < numParams(); ++i)
        _params[i] = 0;

    Vec mid1 = mid - start;
    Vec end1 = end - start;

    double twiceSignedArea = (mid1[0] * end1[1] - mid1[1] * end1[0]);
    double abc = sqrt(mid1.squaredNorm() * end1.squaredNorm() * (mid - end).squaredNorm());

    if(fabs(twiceSignedArea) < 1e-16 || abc < 1e-16)
        return; //degenerate arc

    _params.head<2>() = start;
    _params[CURVATURE] = 2. * twiceSignedArea / abc;
    double halfArcAngle = asin(fabs(0.5 * end1.norm() * _params[CURVATURE]));

    if(mid1.dot(end - mid) < 0.)
        halfArcAngle = PI - halfArcAngle;

    _params[LENGTH] = fabs(2. * halfArcAngle / _params[CURVATURE]);

    if(twiceSignedArea < 0.)
        _params[ANGLE] = AngleUtils::angle(end1) + halfArcAngle;
    else
        _params[ANGLE] = AngleUtils::angle(end1) - halfArcAngle;

    _paramsChanged();    
}

END_NAMESPACE_Cornu
