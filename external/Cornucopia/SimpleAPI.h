/*--
    SimpleAPI.h  

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

#ifndef CORNUCOPIA_SIMPLEAPI_H_INCLUDED
#define CORNUCOPIA_SIMPLEAPI_H_INCLUDED

//The point of this file is to provide a minimalistic API without dependencies on Eigen or anything else.
//To use Cornucopia, you only need to include this file, construct a vector of Points, Parameters, and run the fit(...) function.

#include "Parameters.h"

namespace Cornu
{

struct Point //basic point class
{
    Point() : x(0), y(0) {}
    Point(double inX, double inY) : x(inX), y(inY) {}

    double x;
    double y;
};

struct BasicPrimitive
{
    enum PrimitiveType
    {
        LINE = 0,
        ARC,
        CLOTHOID
    };

    PrimitiveType type;
    Point start;
    double length;
    double startAngle;
    double startCurvature;
    double curvatureDerivative;

    //Evaluates this primitive at (arclength) parameter s (between 0 and length) and returns the position and first and second derivatives.
    //This function is slower than using CurvePrimitive's evaluators
    void eval(double s, Point *outPos, Point *outDer = NULL, Point *outDer2 = NULL) const;
};

//The basic API function: takes a vector of points and a Parameters object (see Parameters.h)
//and returns a vector of primitives and (optionally) whether the curve is closed
std::vector<BasicPrimitive> fit(const std::vector<Point> &points, const Parameters &parameters, bool *outClosed = NULL);

struct BasicBezier
{
    Point controlPoint[4];
};

std::vector<BasicBezier> toBezierSpline(const std::vector<BasicPrimitive> &curve, double tolerance);

} //end of namespace Cornu

#endif //CORNUCOPIA_SIMPLEAPI_H_INCLUDED
