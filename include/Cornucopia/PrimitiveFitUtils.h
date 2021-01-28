/*--
    PrimitiveFitUtils.h  

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

#ifndef CORNUCOPIA_PRIMITIVEFITUTILS_H_INCLUDED
#define CORNUCOPIA_PRIMITIVEFITUTILS_H_INCLUDED

#include "defs.h"
#include "smart_ptr.h"
#include <vector>
#include <Eigen/Core>
#include <Eigen/StdVector>
#include "Line.h"
#include "Arc.h"
#include "Clothoid.h"

NAMESPACE_Cornu

class FitterBase : public smart_base
{
public:
    virtual ~FitterBase() {}

    virtual void addPointW(const Eigen::Vector2d &pt, double weight) = 0;
    virtual void addPoint(const Eigen::Vector2d &pt) { addPointW(pt, 1.); }
    virtual CurvePrimitivePtr getPrimitive() const = 0;
};
CORNU_SMART_TYPEDEFS(FitterBase);

class LineFitter : public FitterBase
{
public:
    typedef Eigen::Vector2d Vec;

    LineFitter() : _numPts(0), _totWeight(0.), _crossSum(0.),
        _firstPoint(Vec::Zero()), _lastPoint(Vec::Zero()), _sum(Vec::Zero()), _squaredSum(Vec::Zero()) {}

    LinePtr getCurve() const;

    //overrides
    void addPointW(const Eigen::Vector2d &pt, double weight);
    CurvePrimitivePtr getPrimitive() const { return getCurve(); }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    Vec _firstPoint, _lastPoint;
    Vec _sum;
    Vec _squaredSum;
    int _numPts;
    double _totWeight;
    double _crossSum;
};

class ArcFitter : public FitterBase
{
public:
    ArcFitter() : _totWeight(0.), _squaredSum(_squaredSum.Zero()), _sum(_sum.Zero()) {}

    ArcPtr getCurve() const;

    //overrides
    void addPointW(const Eigen::Vector2d &pt, double weight);
    CurvePrimitivePtr getPrimitive() const { return getCurve(); }

private:
    Eigen::Matrix3d _squaredSum;
    Eigen::Vector3d _sum;
    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d> > _pts;
    double _totWeight;
};

//This fits a clothoid by fitting a cubic polynomial to the integral of the angle function,
//using its derivative as the angle function of the clothoid,
//and making the clothoid center of mass align with that of the input
class ClothoidFitter : public FitterBase
{
public:
    ClothoidFitter() : _totalLength(0), _prevAngle(0), _angleIntegral(0),
                       _centerOfMass(_centerOfMass.Zero()), _rhs(_rhs.Zero()) {}

    ClothoidPtr getCurve() const;
    ClothoidPtr getCurveWithZeroCurvature(double param) const;

    //overrides
    void addPoint(const Eigen::Vector2d &pt);
    void addPointW(const Eigen::Vector2d &pt, double weight) { addPoint(pt); } //weight is unsupported
    CurvePrimitivePtr getPrimitive() const { return getCurve(); }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    ClothoidPtr getClothoidWithParams(const Eigen::Vector4d &abcd) const;

    static Eigen::Matrix4d _getLhs(double x);
    static Eigen::Vector4d _getRhs(double x, double y, double z);

    Eigen::Vector2d _centerOfMass;
    Eigen::Vector4d _rhs;
    std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d> > _pts;
    double _prevAngle;
    double _angleIntegral;

    double _totalLength;
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_PRIMITIVEFITUTILS_H_INCLUDED
