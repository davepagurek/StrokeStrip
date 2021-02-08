/*--
    SimpleAPI.cpp  

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

#include "SimpleAPI.h"
#include "Cornucopia.h"

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

vector<BasicPrimitive> fit(const vector<Point> &points, const Parameters &parameters, bool *outClosed)
{
    Fitter fitter;
    fitter.setParams(parameters);

    VectorC<Vector2d> pts((int)points.size(), NOT_CIRCULAR);
    for(int i = 0; i < pts.size(); ++i)
        pts[i] = Vector2d(points[i].x, points[i].y);

    //pass it to the fitter and process it
    fitter.setOriginalSketch(new Cornu::Polyline(pts));
    fitter.run();

    //process the output -- count the number of primitives of each type
    PrimitiveSequenceConstPtr output = fitter.finalOutput();

    if(outClosed)
        (*outClosed) = output->isClosed();

    vector<BasicPrimitive> out(output->primitives().size());

    for(int i = 0; i < (int)out.size(); ++i)
    {
        CurvePrimitiveConstPtr cur = output->primitives()[i];
        out[i].type = (BasicPrimitive::PrimitiveType)cur->getType();
        out[i].start = Point(cur->startPos()[0], cur->startPos()[1]);
        out[i].length = cur->length();
        out[i].startAngle = cur->startAngle();
        out[i].startCurvature = cur->startCurvature();
        out[i].curvatureDerivative = 0;
        if(cur->getType() == CurvePrimitive::CLOTHOID)
            out[i].curvatureDerivative = cur->params()[CurvePrimitive::DCURVATURE];
    }

    return out;
}

//converts a BasicPrimitive to a CurvePrimitive
CurvePrimitivePtr _toCurvePrimitive(const BasicPrimitive &primitive)
{
    CurvePrimitive::ParamVec params;
    params.resize(4 + primitive.type, 1);

    params[CurvePrimitive::X] = primitive.start.x;
    params[CurvePrimitive::Y] = primitive.start.y;
    params[CurvePrimitive::ANGLE] = primitive.startAngle;
    params[CurvePrimitive::LENGTH] = primitive.length;

    CurvePrimitivePtr out;

    switch(primitive.type)
    {
    case BasicPrimitive::LINE:
        out = new Line();
        break;
    case BasicPrimitive::ARC:
        out = new Arc();
        params[CurvePrimitive::CURVATURE] = primitive.startCurvature;
        break;
    case BasicPrimitive::CLOTHOID:
        out = new Clothoid();
        params[CurvePrimitive::CURVATURE] = primitive.startCurvature;
        params[CurvePrimitive::DCURVATURE] = primitive.curvatureDerivative;
        break;
    };

    out->setParams(params);

    return out;
}

void BasicPrimitive::eval(double s, Point *outPos, Point *outDer, Point *outDer2) const
{
    CurvePrimitivePtr curvePrimitive = _toCurvePrimitive(*this);

    Eigen::Vector2d pos, der, der2;
    curvePrimitive->eval(s, &pos, &der, &der2);

    if(outPos)
        *outPos = Point(pos[0], pos[1]);
    if(outDer)
        *outDer = Point(der[0], der[1]);
    if(outDer2)
        *outDer2 = Point(der2[0], der2[1]);
}

vector<BasicBezier> toBezierSpline(const vector<BasicPrimitive> &curve, double tolerance)
{
    VectorC<CurvePrimitiveConstPtr> curvePrimitives;
    for(int i = 0; i < (int)curve.size(); ++i)
        curvePrimitives.push_back(_toCurvePrimitive(curve[i]));

    PrimitiveSequence seq(curvePrimitives);

    BezierSplinePtr spline = seq.toBezierSpline(tolerance);

    vector<BasicBezier> out(spline->primitives().size());

    for(int i = 0; i < (int)spline->primitives().size(); ++i)
    {
        for(int j = 0; j < 4; ++j)
        {
            Vector2d pt = spline->primitives()[i].controlPoint(j);
            out[i].controlPoint[j].x = pt[0];
            out[i].controlPoint[j].y = pt[1];
        }
    }

    return out;
}

END_NAMESPACE_Cornu


