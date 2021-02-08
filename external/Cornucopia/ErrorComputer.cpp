/*--
    ErrorComputer.cpp  

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

#include "ErrorComputer.h"
#include "Resampler.h"
#include "Fitter.h"
#include "Polyline.h"
#include "CurvePrimitive.h"

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

class L2ErrorComputer : public ErrorComputer
{
public:
    L2ErrorComputer(const Fitter &fitter)
        : _pts(fitter.output<RESAMPLING>()->output->pts())
    {
        PolylineConstPtr poly = fitter.output<RESAMPLING>()->output;
        CircularType circular = _pts.circular();
        _weightsLeft = VectorC<double>(_pts.size(), circular);
        _weightsRight = VectorC<double>(_pts.size(), circular);
        _weightLeftRoots = VectorC<double>(_pts.size(), circular);
        _weightRightRoots = VectorC<double>(_pts.size(), circular);
        _weightRoots = VectorC<double>(_pts.size(), circular);

        for(int i = 0; i < _pts.size(); ++i)
        {
            VectorC<Vector2d>::Circulator prev = --_pts.circulator(i);
            VectorC<Vector2d>::Circulator next = ++_pts.circulator(i);
            _weightsLeft[i] = prev.done() ? 0 : poly->lengthFromTo(prev.index(), i);
            _weightsRight[i] = next.done() ? 0 : poly->lengthFromTo(i, next.index());
            _weightLeftRoots[i] = sqrt(_weightsLeft[i]);
            _weightRightRoots[i] = sqrt(_weightsRight[i]);
            _weightRoots[i] = sqrt(_weightsLeft[i] + _weightsRight[i]);
        }
    }

    double computeError(CurvePrimitiveConstPtr curve, int from, int to,
                        bool firstToEndpoint, bool lastToEndpoint, bool reversed) const
    {
        double error = 0;

        if(from < 0 || to >= (int)_pts.size())
            return 0.;

        bool first = true;
        for(VectorC<Vector2d>::Circulator circ = _pts.circulator(from); ; ++circ)
        {
            int idx = circ.index();
            bool last = (idx == to);

            bool toFirstEndpoint = first && firstToEndpoint;
            bool toLastEndpoint = last && lastToEndpoint;

            const Vector2d &pt = _pts.flatAt(idx);

            double s;
            if(toLastEndpoint)
                s = reversed ? 0 : curve->length();
            else if(toFirstEndpoint)
                s = reversed ? curve->length() : 0;
            else
                s = curve->project(pt);

            double distSq = (curve->pos(s) - pt).squaredNorm();
            double weight = 0;
            if(!toFirstEndpoint)
                weight += _weightsLeft.flatAt(idx);
            if(!toLastEndpoint)
                weight += _weightsRight.flatAt(idx);

            error += weight * distSq;

            first = false;
            if(last)
                break;
        }

        return error;
    }

    void computeErrorVector(CurvePrimitiveConstPtr curve, int from, int to, VectorXd &outError, MatrixXd *outErrorDer,
                            bool firstToEndpoint, bool lastToEndpoint, bool reversed) const
    {
        int numParams = (int)curve->params().size();
        int numOutputs = 2 * (_pts.numElems(from, to) + 1); //to is inclusive
        outError.resize(numOutputs); 
        if(outErrorDer)
            outErrorDer->resize(numOutputs, numParams);

        if(from < 0 || to >= (int)_pts.size())
        {
            outError.setZero();
            outErrorDer->setZero();
            return;
        }
         
        CurvePrimitive::ParamDer der, tanDer;
        bool first = true;
        int vecIdx = 0;
        for(VectorC<Vector2d>::Circulator circ = _pts.circulator(from); ; ++circ, vecIdx += 2)
        {
            int idx = circ.index();
            bool last = (idx == to);

            bool toFirstEndpoint = first && firstToEndpoint;
            bool toLastEndpoint = last && lastToEndpoint;

            const Vector2d &pt = _pts.flatAt(idx);
            double weightRoot = 0;
            if(toFirstEndpoint)
                weightRoot = _weightRightRoots.flatAt(idx);
            else if(toLastEndpoint)
                weightRoot = _weightLeftRoots.flatAt(idx);
            else
                weightRoot = _weightRoots.flatAt(idx);

            double s;
            if(toLastEndpoint)
                s = reversed ? 0 : curve->length();
            else if(toFirstEndpoint)
                s = reversed ? curve->length() : 0;
            else
                s = curve->project(pt);

            Vector2d err = curve->pos(s) - pt;
            outError.segment<2>(vecIdx) = err * weightRoot;

            if(outErrorDer)
            {
                curve->derivativeAt(s, der, tanDer);
                Vector2d tangent = curve->der(s);
                RowVectorXd ds = RowVectorXd::Zero(numParams); 

                const double tol = 1e-10;

                if(s + tol >= curve->length())
                    ds(CurvePrimitive::LENGTH) = 1.;
                else if(s > tol)
                {
                    double dfds = 1. + curve->der2(s).dot(err);
                    if(fabs(dfds) < tol)
                        dfds = (dfds < 0. ? -tol : tol);
                    ds = -(err.transpose() * tanDer + tangent.transpose() * der) / dfds;
                }

                outErrorDer->block(vecIdx, 0, 2, der.cols()) = (der + tangent * ds) * weightRoot;
            }

            first = false;
            if(last)
                break;
        }
    }

    double computeErrorForCost(CurvePrimitiveConstPtr curve, int from, int to,
                               bool firstToEndpoint, bool lastToEndpoint, bool reversed) const
    {
        return computeError(curve, from, to, firstToEndpoint, lastToEndpoint, reversed) / curve->length();
    }

protected:
    const VectorC<Vector2d> &_pts;
    VectorC<double> _weightsLeft, _weightsRight, _weightLeftRoots, _weightRightRoots, _weightRoots;
};

class LInfErrorComputer : public L2ErrorComputer
{
public:
    LInfErrorComputer(const Fitter &fitter)
        : L2ErrorComputer(fitter) {}

    double computeErrorForCost(CurvePrimitiveConstPtr curve, int from, int to,
                               bool firstToEndpoint, bool lastToEndpoint, bool reversed) const
    {
        double error = 0;

        if(from < 0 || to >= (int)_pts.size())
            return 0.;

        bool first = true;
        for(VectorC<Vector2d>::Circulator circ = _pts.circulator(from); ; ++circ)
        {
            int idx = circ.index();
            bool last = (idx == to);

            bool toFirstEndpoint = first && firstToEndpoint;
            bool toLastEndpoint = last && lastToEndpoint;

            const Vector2d &pt = _pts.flatAt(idx);

            double s;
            if(toLastEndpoint)
                s = reversed ? 0 : curve->length();
            else if(toFirstEndpoint)
                s = reversed ? curve->length() : 0;
            else
                s = curve->project(pt);

            double distSq = (curve->pos(s) - pt).squaredNorm();

            error = max(error, distSq);

            first = false;
            if(last)
                break;
        }

        return error;
    }    
};

class ErrorComputerCreator : public Algorithm<ERROR_COMPUTER>
{
public:
    ErrorComputerCreator(bool lInf = true)
        : _lInf(lInf) {}

    string name() const { return _lInf ? "L-Infinity" : "L2"; }

protected:
    void _run(const Fitter &fitter, AlgorithmOutput<ERROR_COMPUTER> &out)
    {
        out.errorComputer = _lInf ? new LInfErrorComputer(fitter) : new L2ErrorComputer(fitter);
    }
private:
    bool _lInf;
};

void Algorithm<ERROR_COMPUTER>::_initialize()
{
    new ErrorComputerCreator(true);
    new ErrorComputerCreator(false);
}


END_NAMESPACE_Cornu


