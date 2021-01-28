/*--
    GraphConstructor.cpp  

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

#include "GraphConstructor.h"
#include "PrimitiveFitter.h"
#include "Fitter.h"
#include "Resampler.h"
#include "Polyline.h"
#include "CurvePrimitive.h"
#include "Preprocessing.h"
#include "TwoCurveCombine.h"
#include "Oversketcher.h"

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

struct PrimitiveCacheData
{
    double x, y, angle, curvature, param;

    static PrimitiveCacheData make(CurvePrimitiveConstPtr curve, double param)
    {
        PrimitiveCacheData out;
        Vector2d pos = curve->pos(param);
        out.param = param;
        out.x = pos[0];
        out.y = pos[1];
        out.angle = curve->angle(param);
        out.curvature = curve->curvature(param);
        return out;
    }
};

class PrimitiveCache
{
public:
    PrimitiveCache(const PrimitiveCache &other)
    {
        _numVals = other._numVals;
        for(int i = 0; i < 3; ++i)
        {
            _startData[i] = other._startData[i];
            _endData[i] = other._endData[i];
        }
    }

    PrimitiveCache(const Fitter &fitter, const FitPrimitive &primitive)
    {
        const VectorC<Vector2d> &pts = fitter.output<RESAMPLING>()->output->pts();
        CurvePrimitiveConstPtr curve = primitive.curve;

        _numVals = min(3, primitive.numPts);

        for(int i = 0; i < 3; ++i)
        {
            if(i >= primitive.numPts)
                break;

            int realStartIdx = primitive.isFixed() ? 0 : primitive.startIdx;
            double startParam = (i == 0 ? 0. : curve->project(pts[realStartIdx + i]));
            _startData[i] = PrimitiveCacheData::make(curve, startParam);

            int realEndIdx = primitive.isFixed() ? (pts.size() - 1) : primitive.endIdx;
            double endParam = (i == 0 ? curve->length() : curve->project(pts[realEndIdx - i]));
            _endData[i] = PrimitiveCacheData::make(curve, endParam);
        }
    }

    const PrimitiveCacheData &start(int offs) const { assert(offs < _numVals); return _startData[offs]; }
    const PrimitiveCacheData &end(int offs) const { assert(offs < _numVals); return _endData[offs]; }

private:
    int _numVals;
    PrimitiveCacheData _startData[3], _endData[3];
};

class CostEvaluator : public smart_base
{
public:
    CostEvaluator(const Fitter &fitter) :
        _primitives(fitter.output<PRIMITIVE_FITTING>()->primitives),
        _corners(fitter.output<RESAMPLING>()->corners)
    {
        for(int i = 0; i < 3; ++i)
        {
            _curveCost[i] = fitter.params().get(Parameters::ParameterType(i + Parameters::LINE_COST));
            _continuityCost[i] = fitter.params().get(Parameters::ParameterType(i + Parameters::G0_COST));
        }
        _errorCostFactor = fitter.params().get(Parameters::ERROR_COST);
        _inflectionCost = fitter.params().get(Parameters::INFLECTION_COST);
        _lengthScale = 1. / (fitter.output<SCALE_DETECTION>()->scale * fitter.params().get(Parameters::PIXEL_SIZE));
        _shortnessCostFactor = fitter.params().get(Parameters::SHORTNESS_COST);
        _shortnessThreshold = fitter.scaledParameter(Parameters::SHORTNESS_THRESHOLD);

        for(int i = 0; i < (int)_primitives.size(); ++i)
            _primitiveCache.push_back(PrimitiveCache(fitter, _primitives[i]));
    }

    double vertexCost(int p) const
    {
        if(_primitives[p].isFixed())
            return 0.;

        //complexity
        double out = _curveCost[_primitives[p].curve->getType()];

        //error
        out += _errorCost(_primitives[p].error);

        //inflection
        if(_primitives[p].startCurvSign != _primitives[p].endCurvSign)
            out += _inflectionCost;

        //shortness
        double len = _primitives[p].curve->length();
        if(_continuityCost[0] == Parameters::infinity)
        {
            //figure out how much we expect the length to decrease when we join things up
            double lenDecrease = 0;
            if(!_corners[_primitives[p].startIdx])
                lenDecrease += _primitiveCache[p].start(1).param;
            if(!_corners[_primitives[p].endIdx])
                lenDecrease += len - _primitiveCache[p].end(1).param;
            if(_continuityCost[1] != Parameters::infinity)
                lenDecrease *= 0.5;
            len -= lenDecrease;
        }
        double st = _shortnessThreshold;
        out += _shortnessCostFactor * (sqrt(st * st + SQR(max(0., st - len))) - st);

        return out;
    }

    double edgeCost(int p1, int p2, int continuity, double error1 = -1., double error2 = -1.) const
    {
        double out = 0.;

        //continuity
        if(!_corners[_primitives[p1].endIdx]) //no primitive spans a corner, so the start index is also there
            out += _continuityCost[continuity];

        //inflection
        if(continuity > 1 && _primitives[p1].endCurvSign != _primitives[p2].startCurvSign)
            out += _inflectionCost;

        //error
        double err1 = _primitives[p1].error;
        double err2 = _primitives[p2].error;

        if(error1 < 0.) //we need to predict the error
        {
            double extra1, extra2;
            _getExtraError(p1, p2, continuity, extra1, extra2);
        
            out += _errorCost(extra1 + err1) - _errorCost(err1);
            out += _errorCost(extra2 + err2) - _errorCost(err2);
        }
        else //the error has been passed in
        {
            out += max(0., _errorCost(error1) - _errorCost(err1));
            out += max(0., _errorCost(error2) - _errorCost(err2));
        }

        return out;
    }

private:
    double _errorCost(double error) const
    {
        return _errorCostFactor * (error * SQR(_lengthScale)); //simple for now
    }

    void _getExtraError(int p1, int p2, int continuity, double &outExtra1, double &outExtra2) const
    {
        //outExtra1 = outExtra2 = 0.;  return;
        int offset = continuity;
        if(continuity > 0 && _primitives[p1].endIdx == _primitives[p2].startIdx)
            offset = 0; //one of the curve is a start or an end curve
        Vector3d diffs = _getDiffs(p1, p2, offset);
        for(int i = continuity + 1; i < 3; ++i)
            diffs[i] = 0.; //don't count more than necessary

        double len1 = _primitives[p1].curve->length();
        double len2 = _primitives[p2].curve->length();

        outExtra1 = diffs[0] * 0.5 + len1 * diffs[1] * 0.25 + SQR(len1) * diffs[2] * 0.125;
        outExtra2 = diffs[0] * 0.5 + len2 * diffs[1] * 0.25 + SQR(len2) * diffs[2] * 0.125;
    }

    Vector3d _getDiffs(int p1, int p2, int offset) const
    {
        Vector3d out;
        const PrimitiveCache &cache1 = _primitiveCache[p1];
        const PrimitiveCache &cache2 = _primitiveCache[p2];

        double minDistSq = Parameters::infinity;
        double minAngleDiff = TWOPI;
        double maxAngleDiff = -TWOPI;
        double minCurvatureDiff = Parameters::infinity;
        double maxCurvatureDiff = -minCurvatureDiff;
        for(int i = 0; i <= offset; ++i)
        {
            const PrimitiveCacheData &d1 = cache1.end(offset - i);
            const PrimitiveCacheData &d2 = cache2.start(i);

            //position
            minDistSq = min(minDistSq, Vector2d(d1.x - d2.x, d1.y - d2.y).squaredNorm());

            //angle
            double angleDiff = AngleUtils::toRange(d1.angle - d2.angle, -PI);
            minAngleDiff = min(angleDiff, minAngleDiff);
            maxAngleDiff = max(angleDiff, maxAngleDiff);

            //curvature
            double cDiff = d1.curvature - d2.curvature;
            minCurvatureDiff = min(cDiff, minCurvatureDiff);
            maxCurvatureDiff = max(cDiff, maxCurvatureDiff);
        }

        out[0] = sqrt(minDistSq);
        out[1] = (minAngleDiff * maxAngleDiff < 0.) ? 0. : min(fabs(minAngleDiff), fabs(maxAngleDiff));
        out[2] = (minCurvatureDiff * maxCurvatureDiff < 0.) ? 0. : min(fabs(minCurvatureDiff), fabs(maxCurvatureDiff));

        return out;
    }

    struct PointData
    {
        double x, y, angle, curvature;
    };

    const vector<FitPrimitive> &_primitives;
    vector<PrimitiveCache> _primitiveCache;
    const VectorC<bool> &_corners;

    double _curveCost[3];
    double _continuityCost[3];
    double _errorCostFactor;
    double _inflectionCost;
    double _lengthScale;
    double _shortnessCostFactor;
    double _shortnessThreshold;
};

class DefaultGraphConstructor : public Algorithm<GRAPH_CONSTRUCTION>
{
public:
    string name() const { return "Default"; }

protected:
    void _run(const Fitter &fitter, AlgorithmOutput<GRAPH_CONSTRUCTION> &out)
    {
        const vector<FitPrimitive> &primitives = fitter.output<PRIMITIVE_FITTING>()->primitives;
        PolylineConstPtr poly = fitter.output<RESAMPLING>()->output;
        VectorC<bool> corners = fitter.output<RESAMPLING>()->corners;
        smart_ptr<const AlgorithmOutput<OVERSKETCHING> > osOutput = fitter.output<OVERSKETCHING>();
        const VectorC<Vector2d> &pts = poly->pts();
        bool closed = poly->isClosed();

        out.costEvaluator = new CostEvaluator(fitter);

        //create vertices
        out.vertices.resize(primitives.size());
        for(int i = 0; i < (int)primitives.size(); ++i)
        {
            out.vertices[i].primitiveIdx = i;            
            out.vertices[i].source = out.vertices[i].target = false;
            if(!closed)
            {
                if(osOutput->startCurve)
                    out.vertices[i].source = primitives[i].isStartCurve();
                else
                    out.vertices[i].source = (primitives[i].startIdx == 0);

                if(osOutput->endCurve)
                    out.vertices[i].target = primitives[i].isEndCurve();
                else
                    out.vertices[i].target = (primitives[i].endIdx + 1 == pts.size());
            }

            out.vertices[i].cost = (float)out.costEvaluator->vertexCost(i);

            if(out.vertices[i].source && out.vertices[i].target) //one primitive over the entire curve--create dummy edge
            {
                Edge e;
                e.continuity = -1;
                e.startVtx = e.endVtx = i;
                e.cost = out.vertices[i].cost;
                out.vertices[i].edges.push_back((int)out.edges.size());
                out.edges.push_back(e);
            }
        }

        //create edges
        VectorC<vector<int> > curvesStartingAt(pts.size(), pts.circular());
        for(int i = 0; i < (int)primitives.size(); ++i)
            if(!primitives[i].isStartCurve())
                curvesStartingAt[primitives[i].startIdx].push_back(i);

        for(int i = 0; i < (int)primitives.size(); ++i)
        {
            if(primitives[i].isEndCurve()) //no edges from end curves
                continue;

            int endIdx = primitives[i].endIdx;
            int curve1len = primitives[i].numPts - 1;

            for(int continuity = 0; continuity <= 2; ++continuity)
            {
                int offset = continuity;
                int startIdx = endIdx - offset;
                if(!closed && startIdx < 0)
                    continue;
                if(curve1len <= offset * 2) //if the first curve is already too short
                    continue;

                bool firstCurveConstrained = (primitives[i].curve->getType() < continuity) || primitives[i].isFixed();
                int minType = firstCurveConstrained ? continuity : 0;

                for(int j = 0; j < (int)curvesStartingAt[startIdx].size(); ++j)
                {
                    int k = curvesStartingAt[startIdx][j]; //index of the second primitive
                    int curve2len = primitives[k].numPts - 1;
                    if(curve2len <= offset * 2)
                        continue;

                    bool secondCurveConstrained = (primitives[k].curve->getType() < continuity) || primitives[k].isFixed();
                    if(firstCurveConstrained && secondCurveConstrained)
                        continue;

                    //create edge
                    Edge e;
                    e.startVtx = i;
                    e.endVtx = k;
                    e.continuity = continuity;
                    e.cost = (float)out.costEvaluator->edgeCost(i, k, continuity);
                    e.cost += out.vertices[i].cost * (out.vertices[i].source ? 1.f : 0.5f);
                    e.cost += out.vertices[k].cost * (out.vertices[k].target ? 1.f : 0.5f);
                    if(e.cost >= Parameters::infinity)
                        continue;
                    out.vertices[i].edges.push_back((int)out.edges.size());
                    out.edges.push_back(e);

                    if(e.cost != e.cost)
                        Debugging::get()->printf("Error! Nan cost for edge");
                }
            }
        }

        Debugging::get()->printf("Graph vertices = %d edges = %d", out.vertices.size(), out.edges.size());
    }
};

float Edge::validatedCost(const Fitter &fitter) const
{
    if(continuity < 0) //dummy edge
        return cost;

    Combination comb;
    comb = twoCurveCombine(startVtx, endVtx, continuity, fitter);

#if 0
    if(fitter.output<GRAPH_CONSTRUCTION>()->vertices[startVtx].source)
    {
        Debugging::get()->drawCurve(comb.c1, Debugging::Color(1, 0, 0), "Combined");
        Debugging::get()->drawCurve(comb.c2, Debugging::Color(0.8, 0.5, 0), "Combined");
    }
#endif
    
    float newCost;
     newCost = (float)fitter.output<GRAPH_CONSTRUCTION>()->costEvaluator->edgeCost(startVtx, endVtx, continuity, comb.err1, comb.err2);

    //only increase cost
    return max(newCost, cost);
}

void Algorithm<GRAPH_CONSTRUCTION>::_initialize()
{
    new DefaultGraphConstructor();
}

END_NAMESPACE_Cornu


