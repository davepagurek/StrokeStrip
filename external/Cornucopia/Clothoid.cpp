/*--
    Clothoid.cpp  

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

#include "Clothoid.h"
#include "Fresnel.h"
#include "Eigen/LU"

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

Clothoid::Clothoid(const Vec &start, double startAngle, double length, double curvature, double endCurvature)
{
    _params.resize(numParams());

    _params.head<2>() = start;
    _params[ANGLE] = AngleUtils::toRange(startAngle);
    _params[LENGTH] = length;
    _params[CURVATURE] = curvature;
    _params[DCURVATURE] = (endCurvature - curvature) / length;

    _paramsChanged();
}

void Clothoid::_paramsChanged()
{
    Vector2d startcs;

    _arc = fabs(_params[DCURVATURE]) < 1e-12;
    _flat = false;

    if(_arc)
    {
        _flat = fabs(_params[CURVATURE]) < 1e-6;

        if(_flat)
        {
            _t1 = 0;
            _tdiff = 1;

            double angle = _params[ANGLE];
            double cosA = cos(angle), sinA = sin(angle);

            _mat << cosA, -sinA,
                sinA, cosA;

            startcs = Vector2d(0, 0);
        }
        else //non-flat arc
        {
            _t1 = 0;
            _tdiff = _params[CURVATURE];

            double angle = _params[ANGLE] - HALFPI;
            double cosAS = cos(angle) / _params[CURVATURE], sinAS = sin(angle) / _params[CURVATURE];

            _mat << cosAS, -sinAS,
                sinAS, cosAS;

            startcs = Vector2d(1, 0);
        }
    }
    else //clothoid
    {
        double scale = sqrt(fabs(1. / (PI * _params[DCURVATURE])));

        _t1 = _params[CURVATURE] * scale;
        _tdiff = _params[DCURVATURE] * scale;

        if(_tdiff > 0)
        {
            double angleShift = _params[ANGLE] - _t1 * _t1 * HALFPI;
            double cosAS = PI * scale * cos(angleShift), sinAS = PI * scale * sin(angleShift);
            _mat << cosAS, -sinAS,
                sinAS, cosAS;
        }
        else //we need a reflection here
        {
            double angleShift = _params[ANGLE] + _t1 * _t1 * HALFPI;
            double cosAS = PI * scale * cos(angleShift), sinAS = PI * scale * sin(angleShift);
            _mat << -cosAS, -sinAS,
                -sinAS, cosAS;
        }

        fresnel(_t1, &(startcs[1]), &(startcs[0]));
    }

    _startShift = _startPos() - _mat * startcs;    
}

bool Clothoid::isValidImpl() const
{
    if(_params[LENGTH] < 0.)
        return false;
    return true;
}

void Clothoid::eval(double s, Vec *pos, Vec *der, Vec *der2) const
{
    if(pos)
    {
        double t = _t1 + s * _tdiff;

        Vector2d cs;
        if(_flat)
            cs = Vector2d(t, 0);
        else if(_arc)
            cs = Vector2d(cos(t), sin(t));
        else
            fresnel(t, &(cs[1]), &(cs[0]));

        (*pos) = _startShift + _mat * cs;
    }
    if(der || der2)
    {
        double angle = _params[ANGLE] + s * (_params[CURVATURE] + 0.5 * s * _params[DCURVATURE]);
        double cosa = cos(angle), sina = sin(angle);

        if(der)
            (*der) = Vector2d(cosa, sina);
        if(der2)
            (*der2) = (_params[CURVATURE] + s * _params[DCURVATURE]) * Vector2d(-sina, cosa);
    }
}

double Clothoid::angle(double s) const
{
    return _params[ANGLE] + s * (_params[CURVATURE] + 0.5 * s * _params[DCURVATURE]);
}

double Clothoid::curvature(double s) const
{
    return _params[CURVATURE] + s * _params[DCURVATURE];
}

double Clothoid::project(const Vec &point) const
{
    if(_flat)
    {
        Vector2d tangent(cos(_params[ANGLE]), sin(_params[ANGLE]));
        return min(_length(), max(0., tangent.dot(point - _startPos())));
    }
    else if(_arc)
    {
        Vector2d tangent(cos(_params[ANGLE]), sin(_params[ANGLE]));
        Vec toCenter(-tangent[1], tangent[0]);
        Vector2d center = _startPos() + toCenter / _params[CURVATURE];
        double angleDiff = _length() * _params[CURVATURE];

        double projAngle = atan2(point[1] - center[1], point[0] - center[0]);
        //To compute the projection, get the angle difference into the range centered on the midpoint of
        //the arc.
        double t;
        if(_params[CURVATURE] > 0.)
        {
            double projAngleDiff = HALFPI + projAngle - _startAngle();
            t = AngleUtils::toRange(projAngleDiff, 0.5 * angleDiff - PI) / _params[CURVATURE];
        }
        else
        {
            double projAngleDiff = HALFPI - projAngle + _startAngle();
            t = AngleUtils::toRange(projAngleDiff, -0.5 * angleDiff - PI) / -_params[CURVATURE];
        }

        return max(0., min(_length(), t));
    }
    //Go to the canonical clothoid
    double endT = _t1 + _tdiff * _length();
    Vec pt = _mat.inverse() * (point - _startShift);

    double bestT = _clothoidProjector()->project(pt, min(_t1, endT), max(_t1, endT));
    return (bestT - _t1) / _tdiff;
}

void Clothoid::trim(double sFrom, double sTo)
{
    Vec newStart = pos(sFrom);
    _params[ANGLE] = angle(sFrom); //this line must come before the curvature change
    _params[CURVATURE] = curvature(sFrom); //the angle change in the previous line doesn't affect curvature
    _params[LENGTH] = sTo - sFrom;
    _params[X] = newStart[0];
    _params[Y] = newStart[1];

    _paramsChanged();
}

void Clothoid::flip()
{
    Vec newStart = endPos();
    _params[ANGLE] = PI + endAngle();
    _params[CURVATURE] = -endCurvature();
    _params[X] = newStart[0];
    _params[Y] = newStart[1];

    _paramsChanged();
}

void Clothoid::derivativeAt(double s, ParamDer &out, ParamDer &outTan) const
{
    outTan = out = ParamDer::Zero(2, 6);
    out(0, X) = 1;
    out(1, Y) = 1;

    Vec pos, tangent;
    eval(s, &pos, &tangent);

    Vec diff = pos - _startPos();
    out(0, ANGLE) = -diff[1];
    out(1, ANGLE) = diff[0];

    outTan.col(ANGLE) = Vec(-tangent[1], tangent[0]);
    outTan.col(CURVATURE) = s * outTan.col(ANGLE);
    outTan.col(DCURVATURE) = 0.5 * s * s * outTan.col(ANGLE);

    //Now compute derivatives with respect to curvature and the curvature
    //derivative.  These were derived with Mathematica and limits.
    if(_flat)
    {
        Vec nsinCos = Vec(-sin(_params[ANGLE]), cos(_params[ANGLE]));
        out.col(CURVATURE) = (0.5 * s * s) * nsinCos;
        out.col(DCURVATURE) = (0.5 / 3. * s * s * s) * nsinCos;
    }
    else if(_arc)
    {
        double curv = _params[CURVATURE];
        double curAngle = angle(s);
        double cosCur = cos(curAngle), sinCur = sin(curAngle);
        double cosStart = cos(_params[ANGLE]), sinStart = sin(_params[ANGLE]);
        double curvs = curv * s;
        Vec curvDer(curvs * cosCur + sinStart - sinCur, curvs * sinCur + cosCur - cosStart);
        out.col(CURVATURE) = curvDer / (curv * curv);
        Vec dcurvDer(cosStart + (curvs * curvs * 0.5 - 1.) * cosCur - curvs * sinCur,
            sinStart + (curvs * curvs * 0.5 - 1.) * sinCur + curvs * cosCur);
        out.col(DCURVATURE) = dcurvDer / (curv * curv * curv);
    }
    else //non-degenerate
    {
        //This block is half-derived from Mathematica, and half from the considerations below.
        //Unfortunately it's a mess, but it works (see CurveDerivativesTest).
        //Any change to what we need to compute will require rederivation.
        //There's a bunch of optimizations possible, but they'll make the code even worse.
        //
        //dt/dx = dt1/dx + s * dtdiff/dx
        //dcs/dx = cossin(pi t^2 / 2) * dt/dx
        //dstartcs/dx = cossin(pi t1^2 / 2) * dt1/dx
        //dp/dx = dmat/dx * (cs - startcs) + mat * (dcs/dx - dstartcs/dx)

        double t = _t1 + s * _tdiff;
        double scale = sqrt(fabs(1. / (PI * _params[DCURVATURE])));
        RowVector2d dt1dx(scale, -_params[CURVATURE] * scale / (2. * _params[DCURVATURE]));
        RowVector2d dtdx = dt1dx + RowVector2d(0, scale * s * 0.5);
        Vector2d cs, startcs;
        fresnel(t, &(cs[1]), &(cs[0]));
        fresnel(_t1, &(startcs[1]), &(startcs[0]));
        double basicAngle = HALFPI * _t1 * _t1;
        Vector2d dstartcs(cos(basicAngle), sin(basicAngle));
        Vector2d dcs(cos(HALFPI * t * t), sin(HALFPI * t * t));

        Matrix2d result = (_mat * dcs) * dtdx - (_mat * dstartcs) * dt1dx;
        Matrix2d dmatdc, dmatdd;

        double angleShift;
        if(_tdiff > 0.)
            angleShift = _params[ANGLE] - _t1 * _t1 * HALFPI;
        else
            angleShift = _params[ANGLE] + _t1 * _t1 * HALFPI;

        double cosAS = cos(angleShift), sinAS = sin(angleShift);
        dmatdc << sinAS, cosAS,
            -cosAS, sinAS;
        dmatdc *= PI * scale * _params[CURVATURE] / _params[DCURVATURE];
        double curvSqr = _params[CURVATURE] * _params[CURVATURE];
        dmatdd(0, 0) = -_params[DCURVATURE] * cosAS - curvSqr * sinAS;
        dmatdd(1, 1) = dmatdd(0, 0);
        dmatdd(0, 1) = _params[DCURVATURE] * sinAS - curvSqr * cosAS;
        dmatdd(1, 0) = -dmatdd(0, 1);
        dmatdd *= HALFPI * scale / (_params[DCURVATURE] * _params[DCURVATURE]);

        if(_tdiff < 0.)
        {
            dmatdc.col(0) *= -1;
            dmatdd.col(0) *= -1;
        }

        result.col(0) += dmatdc * (cs - startcs);
        result.col(1) += dmatdd * (cs - startcs);

        out.block<2, 2>(0, CURVATURE) = result;
    }
}

void Clothoid::derivativeAtEnd(int continuity, EndDer &out) const
{
    out = EndDer::Zero(2 + continuity, 6);

    ParamDer pDer, dummy;
    derivativeAt(_length(), pDer, dummy);
    Vector2d endTan;
    eval(_length(), NULL, &endTan);

    out.block(0, 0, 2, 6) = pDer;
    out.col(LENGTH).head<2>() = endTan;

    if(continuity >= 1)
    {
        out(2, ANGLE) = 1.;
        out(2, LENGTH) = _params[CURVATURE] + _params[LENGTH] * _params[DCURVATURE];
        out(2, CURVATURE) = _params[LENGTH];
        out(2, DCURVATURE) = 0.5 * SQR(_params[LENGTH]);
    }

    if(continuity == 2)
    {
        out(3, CURVATURE) = 1.;
        out(3, LENGTH) = _params[DCURVATURE];
        out(3, DCURVATURE) = _params[LENGTH];
    }
}

void Clothoid::toEndCurvatureDerivative(MatrixXd &der) const
{
    double invLength = 1. / _params[LENGTH];
    der.col(DCURVATURE) *= invLength;

    //compute the derivative with respect to curvature and length as well.
    der.col(CURVATURE) -= der.col(DCURVATURE);
    der.col(LENGTH) -= der.col(DCURVATURE) * _params(DCURVATURE);
}

END_NAMESPACE_Cornu


