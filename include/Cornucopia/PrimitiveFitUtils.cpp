/*--
    PrimitiveFitUtils.cpp  

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

#include "PrimitiveFitUtils.h"

#include "Line.h"
#include "Arc.h"
#include "Clothoid.h"
#include "Fresnel.h"

#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

void LineFitter::addPointW(const Vector2d &pt, double weight)
{
    ++_numPts;
    _totWeight += weight;

    _lastPoint = pt;
    if(_numPts == 1)
        _firstPoint = pt;

    _sum += weight * pt;
    _squaredSum[0] += weight * SQR(pt[0]);
    _squaredSum[1] += weight * SQR(pt[1]);
    _crossSum += weight * pt[0] * pt[1];
}

LinePtr LineFitter::getCurve() const
{
    if(_numPts < 2)
        return LinePtr();

    Matrix2d cov = Matrix2d::Zero();
    double factor = 1. / double(_totWeight);
    cov(0, 0) = _squaredSum[0] - SQR(_sum[0]) * factor;
    cov(1, 1) = _squaredSum[1] - SQR(_sum[1]) * factor;
    cov(0, 1) = cov(1, 0) = _crossSum - _sum[0] * _sum[1] * factor;
    cov *= factor;

    SelfAdjointEigenSolver<Matrix2d> eigenSolver(cov);

    Vector2d eigVs = eigenSolver.eigenvalues();
    Vector2d dir = eigenSolver.eigenvectors().col(1).normalized(); //1 is the index of the larger eigenvalue
    Vector2d pt = _sum * factor;
    Vector2d pt0 = pt + ((_firstPoint - pt).dot(dir)) * dir;
    Vector2d pt1 = pt + ((_lastPoint - pt).dot(dir)) * dir;

    return new Line(pt0, pt1);
}

void ArcFitter::addPointW(const Vector2d &pt, double weight)
{
    _pts.push_back(pt);

    Vector3d pt3(pt[0] - _pts[0][0], pt[1] - _pts[0][1], (pt - _pts[0]).squaredNorm());

    _totWeight += weight;
    _sum += weight * pt3;
    _squaredSum += weight * pt3 * pt3.transpose();
}

ArcPtr ArcFitter::getCurve() const
{
    if((int)_pts.size() < 2)
        return ArcPtr();

    double factor = 1. / _totWeight;
    Vector3d pt = _sum * factor;
    Matrix3d cov = factor * _squaredSum - pt * pt.transpose();

    SelfAdjointEigenSolver<Matrix3d> eigenSolver(cov);
    Vector3d eigVs = eigenSolver.eigenvalues();

    Vector3d dir = eigenSolver.eigenvectors().col(0); //0 is the index of the smallest eigenvalue
    dir /= (1e-16 + dir[2]);

    double dot = dir.dot(pt);
    //circle equation is:
    //dir[0] * x + dir[1] * y + (x^2+y^2) = dot
    Vector2d center = -0.5 * Vector2d(dir[0], dir[1]);
    double radius = sqrt(1e-16 + dot + center.squaredNorm());
    center += _pts[0];

    //TODO: convert code to use AngleUtils
    //Now get the arc
    Vector2d c[3] = { _pts[0], _pts[_pts.size() / 2], _pts.back() };
    double angle[3];
    for(int i = 0; i < 3; ++i) {
        c[i] = (c[i] - center).normalized() * radius;
        angle[i] = atan2(c[i][1], c[i][0]);
    }
    if(angle[2] < angle[0])
        angle[2] += PI * 2.;
    if(angle[1] < angle[0])
        angle[1] += PI * 2.;
    if(angle[1] <= angle[2]) { //OK--CCW arc
        double a = angle[2] - angle[0];
        return new Arc(c[0] + center, angle[0] + PI * 0.5, a * radius, 1. / radius);
    }
    else { //Backwards--CW arc
        for(int i = 0; i < 3; ++i)
            angle[i] = atan2(c[i][1], c[i][0]);
        if(angle[0] <= angle[2])
            angle[0] += PI * 2.;
        if(angle[1] < angle[2])
            angle[1] += PI * 2.;
        assert(angle[1] <= angle[0]);
        double a = angle[0] - angle[2];
        return new Arc(c[0] + center, angle[0] - PI * 0.5, a * radius, -1. / radius);
    }
}

void ClothoidFitter::addPoint(const Vector2d &pt)
{
    _pts.push_back(pt);

    if((int)_pts.size() < 2)
        return;

    const Vector2d &prevPt = _pts[_pts.size() - 2];
    double segmentLength = (pt - prevPt).norm();
    _centerOfMass += (pt + prevPt) * (0.5 * segmentLength);

    double angle = atan2(pt[1] - prevPt[1], pt[0] - prevPt[0]);
    if(angle < _prevAngle) //make sure it's not far from the previous angle
        angle += TWOPI * int(0.5 + (_prevAngle - angle) / TWOPI);
    else
        angle -= TWOPI * int(0.5 + (angle - _prevAngle) / TWOPI);
    _prevAngle = angle;

    double x0 = _totalLength;
    double x1 = (_totalLength += segmentLength);

    double y = angle;
    double z = _angleIntegral - y * x0;
    _rhs += _getRhs(x1, y, z) - _getRhs(x0, y, z);

    _angleIntegral += segmentLength * angle;
}

ClothoidPtr ClothoidFitter::getCurve() const
{
    Matrix4d lhs;

    lhs = _getLhs(_totalLength);

    Vector4d abcd = lhs.inverse() * _rhs;
    return getClothoidWithParams(abcd);
}

ClothoidPtr ClothoidFitter::getCurveWithZeroCurvature(double param) const
{
    Matrix<double, 5, 5> lhs;
    Matrix<double, 5, 1> rhs;

    Vector4d constraint;
    constraint << 6 * param, 2, 0, 0; // second derivative of ax^3+bx^2+cx+d is 6ax+2b

    //For constrained least squares,
    //lhs is now [A^T A    C]
    //           [ C^T     0]
    lhs << _getLhs(_totalLength),   constraint,
           constraint.transpose(),  0;
    rhs << _rhs, 0;

    Vector4d abcd = (lhs.inverse() * rhs).head<4>();
    return getClothoidWithParams(abcd);
}

ClothoidPtr ClothoidFitter::getClothoidWithParams(const Eigen::Vector4d &abcd) const
{
    Vector3d abc = Vector3d(abcd[0] * 3, abcd[1] * 2, abcd[2]);

    double startAngle = abc[2];
    double startCurvature = abc[1];
    double endCurvature = 2 * abc[0] * _totalLength + abc[1];

    //now compute the center of mass of the clothoid at 0
    double x, y;

    if(fabs(abc[0]) > 1e-8) //if it's a real clothoid
    {
        //The following comes from the expression that Mathematica generates with:
        //Integrate[Integrate[{Cos[a x^2 + b x + c], Sin[a x^2 + b x + c]}, {x, 0, t}], {t, 0, s}]/s
        bool negative = false;
        if(abc[0] < 0)
        {
            negative = true;
            abc = -abc;
        }
        double a = abc[0], b = abc[1], c = abc[2];
        double s = _totalLength;

        double invRtApi2 = 1. / sqrt(0.5 * PI * a);
        
        double f1s, f2s, f1c, f2c;
        fresnelApprox(b * invRtApi2 * 0.5, &f1s, &f1c);
        fresnelApprox((b + 2 * a * s) * invRtApi2 * 0.5, &f2s, &f2c);

        double disc = b * b / (4 * a) - c;
        double sind = sin(disc), cosd = cos(disc);

        y = (cos(c + s * (b + a * s)) - cos(c)) / (2. * a);
        y -= (b + 2 * a * s) * PI * 0.25 * (cosd * (f1s - f2s) + sind * (f2c - f1c)) * invRtApi2 / a;
        if(negative)
            y = -y;

        x = (sin(c) - sin(c + s * (b + a * s))) / (2. * a);
        x += (b + 2 * a * s) * PI * 0.25 * (-sind * (f1s - f2s) + cosd * (f2c - f1c)) * invRtApi2 / a;

        x /= s;
        y /= s;
    }
    else
    {
        double b = abc[1], c = abc[2];
        double s = _totalLength;

        if(fabs(b) < 1e-8) //line
        {
            x = cos(c) * s * 0.5;
            y = sin(c) * s * 0.5;
        }
        else //arc
        {
            double c1 = cos(c), s1 = sin(c);
            double c2 = cos(c + b * s), s2 = sin(c + b * s);

            x = (c1 - c2) / (b * b * s) - s1 / b;
            y = (s1 - s2) / (b * b * s) + c1 / b;
        }
    }

    return new Clothoid(Vector2d(_centerOfMass[0] / _totalLength - x, _centerOfMass[1] / _totalLength - y),
                        startAngle, _totalLength, startCurvature, endCurvature);
}

Matrix4d ClothoidFitter::_getLhs(double x)
{
    double xp[8] = {1, x, 0, 0, 0, 0, 0, 0 }; //powers of totalLength
    for(int i = 2; i < 8; ++i)
        xp[i] = xp[i - 1] * x;

    Matrix4d out;

    for(int i = 0; i < 4; ++i) for(int j = 0; j < 4; ++j)
    {
        int p = 7 - i - j;
        out(i, j) = (2. / double(p)) * xp[p];
    }

    return out;
}

Vector4d ClothoidFitter::_getRhs(double x, double y, double z)
{
    double xp[6] = {1, x, 0, 0, 0, 0 }; //powers of totalLength
    for(int i = 2; i < 6; ++i)
        xp[i] = xp[i - 1] * x;

    Vector4d out;

    out[0] = (2. / 5.) * y * xp[5] + 0.5 * z * xp[4];
    out[1] = 0.5 * y * xp[4] + (2. / 3.) * z * xp[3];
    out[2] = (2. / 3.) * y * xp[3] + z * xp[2];
    out[3] = y * xp[2] + 2. * z * xp[1];

    return out;
}

END_NAMESPACE_Cornu
