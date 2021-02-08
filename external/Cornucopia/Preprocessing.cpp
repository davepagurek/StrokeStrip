/*--
    Preprocessing.cpp  

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

#include "Preprocessing.h"
#include "Fitter.h"
#include "Polyline.h"
#include "Debugging.h"
#include "Line.h"
#include "PiecewiseLinearUtils.h"
#include <iostream>

#include <Eigen/Geometry>

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

class NoScaleDetector : public Algorithm<SCALE_DETECTION>
{
public:
    string name() const { return "None"; }

protected:
    void _run(const Fitter &fitter, AlgorithmOutput<SCALE_DETECTION> &out)
    {
        out.scale = 1.;
    }
};

class AdaptiveScaleDetector : public Algorithm<SCALE_DETECTION>
{
public:
    string name() const { return "Adaptive"; }

protected:
    void _run(const Fitter &fitter, AlgorithmOutput<SCALE_DETECTION> &out)
    {
        double pixel = fitter.params().get(Parameters::PIXEL_SIZE);

        //compute the diagonal length of the bounding box of the curve
        const VectorC<Vector2d> &pts = fitter.originalSketch()->pts();
        AlignedBox<double, 2> boundingBox(pts[0]);
        for(int i = 1; i < (int)pts.size(); ++i)
            boundingBox.extend(pts[i]);
        double diag = boundingBox.diagonal().norm();

        const double smallCurvePixels = fitter.params().get(Parameters::SMALL_CURVE_PIXELS);
        const double largeCurvePixels = fitter.params().get(Parameters::LARGE_CURVE_PIXELS);
        const double maxRescale = fitter.params().get(Parameters::MAX_RESCALE);

        out.scale = 1.;

        if(diag < pixel * smallCurvePixels)
            out.scale = max(1. / maxRescale, diag / (pixel * smallCurvePixels));
        if(diag > pixel * largeCurvePixels)
            out.scale = min(maxRescale, diag / (pixel * largeCurvePixels));
    }
};

void Algorithm<SCALE_DETECTION>::_initialize()
{
    new AdaptiveScaleDetector();
    new NoScaleDetector();
};

//========================Preliminary Resampling=============================

class DefaultPrelimResampling : public Algorithm<PRELIM_RESAMPLING>
{
public:
    string name() const { return "Default"; }

protected:
    void _run(const Fitter &fitter, AlgorithmOutput<PRELIM_RESAMPLING> &out)
    {
        const VectorC<Vector2d> &pts = fitter.originalSketch()->pts();
        VectorC<Vector2d> outPts(0, NOT_CIRCULAR);
        PiecewiseLinearMonotone origToCur(PiecewiseLinearMonotone::POSITIVE);

        const double radius = fitter.scaledParameter(Parameters::MIN_PRELIM_LENGTH);

        //First mark the points on the original curve we definitely want to keep (e.g., corners) using Douglas-Peucker
        vector<bool> keep = _markDouglasPeucker(pts, fitter.scaledParameter(Parameters::DP_CUTOFF));

        //Go through the points and resample them at the rate of radius, discarding those that are too close
        //to the previous output point, and subdiving the line segments between points that are too far.
        //Note that we don't want to simply resample the curve by the arclength parameterization, because
        //noisy regions with too much arclength will get too many samples.
        outPts.push_back(pts[0]);
        origToCur.add(0, 0);
        double lengthSoFar = 0;
        for(int idx = 1; idx < pts.size(); ++idx)
        {
            Vector2d curResampled = outPts.back();
            Vector2d curPt = pts[idx];

            double distSqToCur = (curPt - curResampled).squaredNorm();

            //If it's the last point or was marked for keeping, just output it
            if(idx + 1 == pts.size() || keep[idx])
            {
                lengthSoFar += sqrt(distSqToCur);
                origToCur.add(fitter.originalSketch()->idxToParam(idx), lengthSoFar);
                if(distSqToCur > 1e-16)
                    outPts.push_back(curPt);
                continue;
            }

            //If this point is within radius of the last output point (i.e., too close), discard it.
            if(distSqToCur < SQR(radius))
                continue;

            //The previous point has to be within radius of the previous output point
            Vector2d prevPt = pts[idx - 1];
            
            //Find the intersection of the line segment between the previous and the current point with the circle
            //centered at the last output point whose radius is "radius".
            ParametrizedLine<double, 2> line = ParametrizedLine<double, 2>::Through(prevPt, curPt);

            double projectedLineParam = line.direction().dot(curResampled - line.origin());
            Vector2d closestOnLine = line.origin() + line.direction() * projectedLineParam;
            double distSqToLine = (curResampled - closestOnLine).squaredNorm();
            if(distSqToLine < 1e-16)
            {
                lengthSoFar += sqrt(distSqToCur);
                origToCur.add(fitter.originalSketch()->idxToParam(idx), lengthSoFar);
                outPts.push_back(curPt);
                continue;
            }
            double y = sqrt(max(0., radius * radius - distSqToLine));
            Vector2d newPt = closestOnLine + y * line.direction(); //y is positive, so this will find the correct point

            lengthSoFar += (curResampled - newPt).norm();
            origToCur.add(fitter.originalSketch()->idxToParam(idx - 1) + projectedLineParam + y, lengthSoFar);
            outPts.push_back(newPt);
            --idx;
        }

        out.output = new Polyline(outPts);
        out.parameters.resize(pts.size());
        for(int i = 0; i < pts.size(); ++i)
        {
            double paramOrig = fitter.originalSketch()->idxToParam(i);
            double paramNew;
            if(!origToCur.eval(paramOrig, paramNew))
                Debugging::get()->printf("Evaluation error!");
            out.parameters[i] = paramNew;
            //Debugging::get()->drawLine(pts[i], out.output->pos(paramNew), Vector3d(1, 0, 1), "Correspondence");
        }

        for(int i = 0; i < (int)outPts.size(); ++i)
            Debugging::get()->drawPoint(outPts[i], Vector3d(0, (i % 10 == 0) ? 0.6 : 0, 1), "Prelim resampled");
        Debugging::get()->drawCurve(out.output, Vector3d(0, 0, 1), "Prelim resampled curve");
    }

private:
    vector<bool> _markDouglasPeucker(const VectorC<Vector2d> &pts, double cutoff)
    {
        vector<bool> out(pts.size(), false);
        out.front() = out.back() = true;

        _dpHelper(pts, out, cutoff, 0, pts.size());
        return out;
    }

    void _dpHelper(const VectorC<Vector2d> &pts, vector<bool> &out, double cutoff, int start, int end)
    {
        Line cur(pts[start], pts[end - 1]);
        double maxDist = cutoff;
        int mid = -1;
        for(int i = start; i < end; ++i)
        {
            double dist = cur.distanceTo(pts[i]);
            if(dist > maxDist)
            {
                maxDist = dist;
                mid = i;
            }
        }
        if(mid < 0)
            return;
        out[mid] = true;
        _dpHelper(pts, out, cutoff, start, mid + 1);
        _dpHelper(pts, out, cutoff, mid, end);
    }
};

class NoPrelimResampling : public Algorithm<PRELIM_RESAMPLING>
{
public:
    string name() const { return "None"; }

protected:
    void _run(const Fitter &fitter, AlgorithmOutput<PRELIM_RESAMPLING> &out)
    {
        out.output = fitter.originalSketch();
        out.parameters.resize(fitter.originalSketch()->pts().size());
        for(int i = 0; i < (int)out.parameters.size(); ++i)
            out.parameters[i] = fitter.originalSketch()->idxToParam(i);
    }
};

void Algorithm<PRELIM_RESAMPLING>::_initialize()
{
    new DefaultPrelimResampling();
    new NoPrelimResampling();
}

//========================Curve Closer=============================

class OldCurveCloser : public Algorithm<CURVE_CLOSING>
{
public:
    string name() const { return "Old"; }

protected:
    void _run(const Fitter &fitter, AlgorithmOutput<CURVE_CLOSING> &out)
    {
        out.output = fitter.output<PRELIM_RESAMPLING>()->output;
        out.parameters = fitter.output<PRELIM_RESAMPLING>()->parameters;
        out.closed = false;
        
        if(fitter.oversketchBase())
            return; //if we're oversketching, the curve is not closed

        const double tol = fitter.scaledParameter(Parameters::CLOSEDNESS_THRESHOLD);
        const double tolSq = SQR(tol);
        const VectorC<Vector2d> &pts = fitter.output<PRELIM_RESAMPLING>()->output->pts();

        if(pts.size() < 3)
            return; //definitely not closed

        Vector2d start = pts[0], end = pts.back();

        //find point farthest away from start-end line segment
        Line startEnd(start, end);

        if(startEnd.length() > tol) //start and end too far
            return;

        double maxDistSq = 0;
        int farthest = 1;
        for(int i = 1; i + 1 < (int)pts.size(); ++i) {
            double t = startEnd.project(pts[i]);
            double distSq = (pts[i] - startEnd.pos(t)).squaredNorm();
            if(distSq > maxDistSq) {
                maxDistSq = distSq;
                farthest = i;
            }
        }

        if(maxDistSq < tolSq)
            return;

        //now find the points closest to each other that are within tol of each
        //other and of the segment connecting the start and end points
        int closest0 = 0, closest1 = (int)pts.size() - 1;
        double minDistSq = tolSq;
        for(int i = 0; i < farthest; ++i) {
            double t = startEnd.project(pts[i]);
            if((pts[i] - startEnd.pos(t)).squaredNorm() > tolSq)
                break;
            for(int j = (int)pts.size() - 1; j > farthest; --j) {
                double t2 = startEnd.project(pts[j]);
                if((pts[j] - startEnd.pos(t2)).squaredNorm() > tolSq)
                    break;
                double distSq = (pts[i] - pts[j]).squaredNorm();
                if(distSq > minDistSq)
                    continue;
                minDistSq = distSq;
                closest0 = i;
                closest1 = j;
            }
        }

        if(minDistSq == tolSq)
            return;
        //if we are here, curve is closed
        out.closed = true;

        //close and adjust curve
        VectorC<Vector2d>::Base newPts = pts;

        if(closest1 + 1 < (int)pts.size())
            newPts.erase(newPts.begin() + closest1 + 1, newPts.end());
        if(closest0 > 0)
            newPts.erase(newPts.begin(), newPts.begin() + closest0);

        PolylinePtr polyptr = new Polyline(VectorC<Vector2d>(newPts, NOT_CIRCULAR));
        Vector2d adjust = newPts[0] - newPts.back();
        const double adjustPower = 3.;
        for(int i = 0; i < (int)newPts.size(); ++i) {
            double param = polyptr->idxToParam(i) / polyptr->length();
            double adjustFactor = 0.5 * pow(fabs(param - 0.5) * 2., adjustPower) * (param > 0.5 ? 1. : -1.);
            newPts[i] += adjust * adjustFactor;
        }
        newPts.pop_back(); //the last point is duplicate

        out.output = new Polyline(VectorC<Vector2d>(newPts, CIRCULAR));

        //adjust parameters
        PiecewiseLinearMonotone prevToCur(PiecewiseLinearMonotone::POSITIVE);
        PolylineConstPtr prevOutput = fitter.output<PRELIM_RESAMPLING>()->output;

        for(int i = 0; i < (int)pts.size(); ++i)
        {
            double prevParam = prevOutput->idxToParam(i);
            if(i < closest0)
                prevToCur.add(prevParam, 0);
            else if(i > closest1)
                prevToCur.add(prevParam, out.output->length());
            else
                prevToCur.add(prevParam, out.output->idxToParam(i - closest0));
        }
        prevToCur.batchEval(out.parameters);

#if 0
        for(int i = 0; i < (int)out.parameters.size(); ++i)
            Debugging::get()->drawLine(fitter.originalSketch()->pts()[i], out.output->pos(out.parameters[i]), Vector3d(1, 0, 1), "Closed Correspondence");
#endif

        Debugging::get()->drawCurve(out.output, Debugging::Color(0., 0., 0.), "Closed", 2., Debugging::DOTTED);
    }
};

class UnscaledCurveCloser : public Algorithm<CURVE_CLOSING>
{
public:
    string name() const { return "Scaled"; }

protected:
    void _run(const Fitter &fitter, AlgorithmOutput<CURVE_CLOSING> &out)
    {
        out.output = fitter.output<PRELIM_RESAMPLING>()->output;
        out.parameters = fitter.output<PRELIM_RESAMPLING>()->parameters;
        out.closed = false;
        
        if(fitter.oversketchBase())
            return; //if we're oversketching, the curve is not closed

		const double tol = fitter.params().get(Parameters::CLOSEDNESS_THRESHOLD);
        const double tolSq = SQR(tol);
        const VectorC<Vector2d> &pts = fitter.output<PRELIM_RESAMPLING>()->output->pts();

        if(pts.size() < 3)
            return; //definitely not closed

        Vector2d start = pts[0], end = pts.back();

        //find point farthest away from start-end line segment
        Line startEnd(start, end);

        if(startEnd.length() > tol) //start and end too far
            return;

        double maxDistSq = 0;
        int farthest = 1;
        for(int i = 1; i + 1 < (int)pts.size(); ++i) {
            double t = startEnd.project(pts[i]);
            double distSq = (pts[i] - startEnd.pos(t)).squaredNorm();
            if(distSq > maxDistSq) {
                maxDistSq = distSq;
                farthest = i;
            }
        }

        if(maxDistSq < tolSq)
            return;

        //now find the points closest to each other that are within tol of each
        //other and of the segment connecting the start and end points
        int closest0 = 0, closest1 = (int)pts.size() - 1;
        double minDistSq = tolSq;
        for(int i = 0; i < farthest; ++i) {
            double t = startEnd.project(pts[i]);
            if((pts[i] - startEnd.pos(t)).squaredNorm() > tolSq)
                break;
            for(int j = (int)pts.size() - 1; j > farthest; --j) {
                double t2 = startEnd.project(pts[j]);
                if((pts[j] - startEnd.pos(t2)).squaredNorm() > tolSq)
                    break;
                double distSq = (pts[i] - pts[j]).squaredNorm();
                if(distSq > minDistSq)
                    continue;
                minDistSq = distSq;
                closest0 = i;
                closest1 = j;
            }
        }

        if(minDistSq == tolSq)
            return;
        //if we are here, curve is closed
        out.closed = true;

        //close and adjust curve
        VectorC<Vector2d>::Base newPts = pts;

        if(closest1 + 1 < (int)pts.size())
            newPts.erase(newPts.begin() + closest1 + 1, newPts.end());
        if(closest0 > 0)
            newPts.erase(newPts.begin(), newPts.begin() + closest0);

        PolylinePtr polyptr = new Polyline(VectorC<Vector2d>(newPts, NOT_CIRCULAR));
        Vector2d adjust = newPts[0] - newPts.back();
        const double adjustPower = 3.;
        for(int i = 0; i < (int)newPts.size(); ++i) {
            double param = polyptr->idxToParam(i) / polyptr->length();
            double adjustFactor = 0.5 * pow(fabs(param - 0.5) * 2., adjustPower) * (param > 0.5 ? 1. : -1.);
            newPts[i] += adjust * adjustFactor;
        }
        newPts.pop_back(); //the last point is duplicate

        out.output = new Polyline(VectorC<Vector2d>(newPts, CIRCULAR));

        //adjust parameters
        PiecewiseLinearMonotone prevToCur(PiecewiseLinearMonotone::POSITIVE);
        PolylineConstPtr prevOutput = fitter.output<PRELIM_RESAMPLING>()->output;

        for(int i = 0; i < (int)pts.size(); ++i)
        {
            double prevParam = prevOutput->idxToParam(i);
            if(i < closest0)
                prevToCur.add(prevParam, 0);
            else if(i > closest1)
                prevToCur.add(prevParam, out.output->length());
            else
                prevToCur.add(prevParam, out.output->idxToParam(i - closest0));
        }
        prevToCur.batchEval(out.parameters);

#if 0
        for(int i = 0; i < (int)out.parameters.size(); ++i)
            Debugging::get()->drawLine(fitter.originalSketch()->pts()[i], out.output->pos(out.parameters[i]), Vector3d(1, 0, 1), "Closed Correspondence");
#endif

        Debugging::get()->drawCurve(out.output, Debugging::Color(0., 0., 0.), "Closed", 2., Debugging::DOTTED);
    }
};


void Algorithm<CURVE_CLOSING>::_initialize()
{
    new OldCurveCloser();
	new UnscaledCurveCloser();
}

END_NAMESPACE_Cornu


