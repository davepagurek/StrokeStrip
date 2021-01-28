/*--
    Oversketcher.cpp  

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

#include "Oversketcher.h"
#include "Fitter.h"
#include "Preprocessing.h"
#include "PrimitiveSequence.h"
#include "PiecewiseLinearUtils.h"
#include "Polyline.h"

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

class DefaultOversketcher : public Algorithm<OVERSKETCHING>
{
public:
    string name() const { return "Default"; }

protected:
    void _run(const Fitter &fitter, AlgorithmOutput<OVERSKETCHING> &out)
    {
        PolylineConstPtr curve = fitter.output<CURVE_CLOSING>()->output;
        PrimitiveSequenceConstPtr base = fitter.oversketchBase();

        out.finallyClose = curve->isClosed();
        out.startContinuity = out.endContinuity = 2;
        out.parameters = fitter.output<CURVE_CLOSING>()->parameters;

        if(!base)
        {
            out.output = curve;
            return;
        }

        const double threshold = fitter.params().get(Parameters::OVERSKETCH_THRESHOLD);
        const double peelback = 2. * threshold;
        double startParam = base->project(curve->startPos());
        double endParam = base->project(curve->endPos());
        Vector2d startOnBase = base->pos(startParam);
        Vector2d endOnBase = base->pos(endParam);
        double startDist = (curve->startPos() - startOnBase).norm();
        double endDist = (curve->endPos() - endOnBase).norm();

        //determine whether to flip the base curve--that will determine which part of it remains
        //let's say don't ever flip for now

        //figure out which ends are "close" -- i.e., their projections on the base curve will be used
        bool startClose, endClose;

        //if the base is closed, both ends must be "close"
        if(base->isClosed())
        {
            if(startDist > threshold || endDist > threshold)
            {
                out.output = curve;
                return;
            }
            startClose = endClose = true;
        }
        else
        {
            startClose = startDist < threshold && startDist + peelback < startParam;
            endClose = endDist < threshold && endDist + peelback < base->length() - endParam;

            if(!(startClose || endClose)) //no point oversketching if neither end is close
            {
                out.output = curve;
                return;
            }
        }

        //compute transition points
        double startBaseTransition = 0, startSketchTransition = 0;
        double endBaseTransition = 0, endSketchTransition = 0;

        double remainingLength = startParam - endParam; //how much room we have to transition along the base curve
        if(!startClose || !endClose || (!base->isClosed() && remainingLength <= 0.))
            remainingLength = base->length() * 2.;
        if(base->isClosed() && remainingLength <= 0.)
            remainingLength += base->length();

        double underHalf = 0.49;

        if(startClose)
        {
            startSketchTransition = min(startDist, curve->length() * underHalf);
            startBaseTransition = startParam - min(startDist + peelback, remainingLength * underHalf);
        }
        else
            startSketchTransition = 0;

        if(endClose)
        {
            endSketchTransition = max(curve->length() - endDist, curve->length() * underHalf);
            endBaseTransition = endParam + min(endDist + peelback, remainingLength * underHalf);
        }
        else
            endSketchTransition = curve->length();

        PiecewiseLinearMonotone prevToCur(PiecewiseLinearMonotone::POSITIVE);

        //construct the transitions
        VectorC<Vector2d> cur = curve->trimmed(startSketchTransition, endSketchTransition)->pts();
        VectorC<Vector2d> pre, post;
        if(startClose)
        {
            VectorC<Vector2d> prepre = buildTransition(base, base,
                                                       startBaseTransition, (startBaseTransition + startParam) * 0.5,
                                                       startBaseTransition, (startBaseTransition + startParam) * 0.5, (int)threshold);
            pre = buildTransition(base, curve, (startBaseTransition + startParam) * 0.5, startParam, 0, startSketchTransition, (int)threshold);
            pre.insert(pre.begin(), prepre.begin(), prepre.end() - 1);
            pre.pop_back();
            prevToCur.add(0, length(prepre));
        }
        prevToCur.add(startSketchTransition, length(pre));
        prevToCur.add(endSketchTransition, length(pre) + length(cur));
        if(endClose)
        {
            cur.pop_back();
            post = buildTransition(curve, base, endSketchTransition, curve->length(), endParam, (endBaseTransition + endParam) * 0.5, (int)threshold);
            VectorC<Vector2d> postpost = buildTransition(base, base,
                                                         (endBaseTransition + endParam) * 0.5, endBaseTransition,
                                                         (endBaseTransition + endParam) * 0.5, endBaseTransition, (int)threshold);
            prevToCur.add(curve->length(), length(pre) + length(cur) + length(post));
            post.insert(post.end(), postpost.begin() + 1, postpost.end());
        }
        pre.insert(pre.end(), cur.begin(), cur.end());
        pre.insert(pre.end(), post.begin(), post.end());
        out.output = new Polyline(pre);
        prevToCur.batchEval(out.parameters);

        //construct the trimmed toAppend and toPrepend curves
        if(base->isClosed() || (startClose && endClose && startBaseTransition > endBaseTransition))
        {
            out.toAppend = base->trimmed(endBaseTransition, startBaseTransition);
            out.finallyClose = true;
            out.startCurve = out.toAppend->primitives().back();
            out.endCurve = out.toAppend->primitives()[0];
        }
        else
        {
            out.finallyClose = false;
            if(startClose)
            {
                out.toPrepend = base->trimmed(0, startBaseTransition);
                out.startCurve = out.toPrepend->primitives().back();
            }
            if(endClose)
            {
                out.toAppend = base->trimmed(endBaseTransition, base->length());
                out.endCurve = out.toAppend->primitives()[0];
            }
        }

        Debugging::get()->drawCurve(base, Debugging::Color(0, 0, 0), "Base Curve", 1., Debugging::DASHED);
        Debugging::get()->drawCurve(out.output, Debugging::Color(1, 0, 1), "Oversketch Modified");
        if(out.startCurve)
            Debugging::get()->drawPrimitive(out.startCurve, "Start Curve", 0, 2);
        if(out.endCurve)
            Debugging::get()->drawPrimitive(out.endCurve, "End Curve", 0, 2);
        if(out.toAppend)
            Debugging::get()->drawCurve(out.toAppend, Debugging::Color(1, 0, 0), "To Append");
        if(out.toPrepend)
            Debugging::get()->drawCurve(out.toPrepend, Debugging::Color(1, 0, 0), "To Prepend");
    }

    static VectorC<Vector2d> buildTransition(CurveConstPtr from, CurveConstPtr to, double fromStart, double fromEnd, double toStart, double toEnd, int numSteps)
    {
        VectorC<Vector2d> out(numSteps, NOT_CIRCULAR);
        for(int i = 0; i < numSteps; ++i)
        {
            double t = double(i) / double(numSteps - 1);
            double fromT = fromStart + t * (fromEnd - fromStart);
            double toT = toStart + t * (toEnd - toStart);
            double smoothT = t * t * (3. - 2. * t);
            out[i] = from->pos(fromT) * (1. - smoothT) + to->pos(toT) * smoothT;
        }

        return out;
    }

    static double length(const VectorC<Vector2d> &pts)
    {
        double result = 0;
        for(int i = 0; i < pts.endIdx(1); ++i)
            result += (pts[i] - pts[i + 1]).norm();
        return result;
    }
};

void Algorithm<OVERSKETCHING>::_initialize()
{
    new DefaultOversketcher();
}

END_NAMESPACE_Cornu


