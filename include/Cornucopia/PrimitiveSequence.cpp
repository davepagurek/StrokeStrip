/*--
    PrimitiveSequence.cpp  

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

#include "PrimitiveSequence.h"
#include "Bezier.h"

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

PrimitiveSequence::PrimitiveSequence(const VectorC<CurvePrimitiveConstPtr> &primitives) : _primitives(primitives)
{
    assert(_primitives.size() > 0);
    _lengths.resize(_primitives.size() + 1, 0);
    for(int i = 0; i < (int)_primitives.size(); ++i)
        _lengths[i + 1] = _lengths[i] + _primitives[i]->length();
}

int PrimitiveSequence::paramToIdx(double param, double *outParam) const
{
    int idx = (int)min(std::upper_bound(_lengths.begin(), _lengths.end(), param) - _lengths.begin(), (ptrdiff_t)_lengths.size() - 1) - 1;
    if(outParam)
        *outParam = param - _lengths[idx];
    return idx;
}

void PrimitiveSequence::eval(double s, Vec *pos, Vec *der, Vec *der2) const
{
    if(_primitives.circular())
    {
        s = fmod(s, _lengths.back());
        if(s < 0.)
            s += _lengths.back();
    }
    double cParam;
    int idx = paramToIdx(s, &cParam);
    _primitives[idx]->eval(cParam, pos, der, der2);
}

double PrimitiveSequence::project(const Vector2d &point) const
{
    double bestS = 0.;
    double minDistSq = 1e50;

    for(int i = 0; i < _primitives.size(); ++i)
    {
        double localS = _primitives[i]->project(point);
        Vector2d pt = _primitives[i]->pos(localS);
        double distSq = (pt - point).squaredNorm();
        if(distSq < minDistSq)
        {
            minDistSq = distSq;
            bestS = _lengths[i] + localS;
        }
    }

    return bestS;
}

PrimitiveSequencePtr PrimitiveSequence::trimmed(double from, double to) const
{
    double len = length();

    const double tol = 1e-10;

    if(isClosed())
    {
        if(to - from > len) //trim too long--trim around the middle
        {
            double mid = 0.5 * (to + from);
            from = mid - 0.5 * len + tol;
            to = mid + 0.5 * len - tol;
        }

        //get the arguments into range
        from = fmod(from, len);
        if(from < 0)
            from += len;
        to = fmod(to, len);
        if(to < 0)
            to += len;
    }
    else //if not closed
    {
        from = max(0., from);
        to = min(len, to);
        if(from > to)
            swap(from, to);
    }

    VectorC<CurvePrimitiveConstPtr> out(0, NOT_CIRCULAR);

    double startParamRemainder, endParamRemainder;
    int startIdx = paramToIdx(from, &startParamRemainder);

    if(_primitives[startIdx]->length() - startParamRemainder < tol)
    {
        startIdx = _primitives.toLinearIdx(startIdx + 1);
        startParamRemainder = 0;
    }

    int endIdx = paramToIdx(to, &endParamRemainder);
    if(endParamRemainder < tol)
    {
        endIdx = _primitives.toLinearIdx(endIdx - 1);
        endParamRemainder = _primitives[endIdx]->length();
    }

    if(startIdx == endIdx && startParamRemainder < endParamRemainder)
    {
        out.push_back(_primitives[startIdx]->trimmed(startParamRemainder, endParamRemainder));
    }
    else
    {
        out.push_back(_primitives[startIdx]->trimmed(startParamRemainder, _primitives[startIdx]->length()));

        for(VectorC<CurvePrimitiveConstPtr>::Circulator circ = _primitives.circulator(startIdx + 1); !circ.done() ; ++circ)
        {
            if(circ.index() == endIdx)
                break;

            out.push_back(*circ);
        }

        out.push_back(_primitives[endIdx]->trimmed(0, endParamRemainder));
    }

    return new PrimitiveSequence(out);
}

PrimitiveSequencePtr PrimitiveSequence::flipped() const
{
    VectorC<CurvePrimitiveConstPtr> out = _primitives;
    reverse(out.begin(), out.end());

    for(int i = 0; i < (int)out.size(); ++i)
        out.flatAt(i) = out.flatAt(i)->flipped();

    return new PrimitiveSequence(out);
}

BezierSplinePtr PrimitiveSequence::toBezierSpline(double tolerance) const
{
    BezierSpline::PrimitiveVector segments;

    const int tolerancePts = 3; //compute difference at this many points
    const int maxBeziers = 5; //every primitive will be represented by at most this many beziers

    for(int i = 0; i < (int)_primitives.size(); ++i) //loop over primitives
    {
        int j;
        for(j = 1; j < maxBeziers; ++j) //determine how many bezier segments we need
        {
            double step = _primitives[i]->length() / j;
            bool errorTooHigh = false;

            for(int k = 0; k < j; ++k) //for each bezier segment
            {
                double start = step * k;
                double end = start + step;
                CubicBezier cur = CubicBezier::hermite(_primitives[i]->pos(start), _primitives[i]->pos(end),
                                                       _primitives[i]->der(start) * step, _primitives[i]->der(end) * step);

                for(int m = 0; m < tolerancePts; ++m) //for each tolerance test point
                {
                    double t = double(m + 1) / double(tolerancePts + 1);
                    Vector2d splinePt;
                    cur.eval(t, &splinePt);
                    if((splinePt - _primitives[i]->pos(start + step * t)).squaredNorm() > tolerance * tolerance)
                    {
                        errorTooHigh = true;
                        break;
                    }
                }

                if(errorTooHigh)
                    break;
            }

            if(!errorTooHigh)
                break; //j now contains the right number of segments
        }

        double step = _primitives[i]->length() / j;
        for(int k = 0; k < j; ++k) //for each bezier segment
        {
            double start = step * k;
            double end = start + step;
            CubicBezier cur = CubicBezier::hermite(_primitives[i]->pos(start), _primitives[i]->pos(end),
                                                    _primitives[i]->der(start) * step, _primitives[i]->der(end) * step);
            segments.push_back(cur);
        }
    }

    return new BezierSpline(segments);
}

END_NAMESPACE_Cornu


