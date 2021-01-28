/*--
    PrimitiveFitter.cpp  

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

#include "PrimitiveFitter.h"
#include "Resampler.h"
#include "Fitter.h"
#include "Polyline.h"
#include "PrimitiveFitUtils.h"
#include "ErrorComputer.h"
#include "Solver.h"
#include "Oversketcher.h"

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

class OneCurveProblem : public LSProblem
{
public:
    OneCurveProblem(const FitPrimitive &primitive, ErrorComputerConstPtr errorComputer)
        : _primitive(primitive), _errorComputer(errorComputer)  {}

    //overrides
    double error(const VectorXd &x, LSEvalData *)
    {
        setParams(x);
        return _errorComputer->computeError(_primitive.curve, _primitive.startIdx, _primitive.endIdx);
    }

    LSEvalData *createEvalData()
    {
        return new LSDenseEvalData();
    }
    void eval(const VectorXd &x, LSEvalData *data)
    {
        LSDenseEvalData *curveData = static_cast<LSDenseEvalData *>(data);
        setParams(x);
        MatrixXd &errDer = curveData->errDerRef();
        _errorComputer->computeErrorVector(_primitive.curve, _primitive.startIdx, _primitive.endIdx,
                                           curveData->errVectorRef(), &errDer);

        _primitive.curve->toEndCurvatureDerivative(errDer);
    }

    VectorXd params() const
    {
        if(_primitive.curve->getType() != CurvePrimitive::CLOTHOID)
            return _primitive.curve->params();
        VectorXd out = _primitive.curve->params();
        out(CurvePrimitive::DCURVATURE) = out(CurvePrimitive::CURVATURE) + out(CurvePrimitive::LENGTH) * out(CurvePrimitive::DCURVATURE);
        return out;
    }

    void setParams(const VectorXd &x)
    {
        if(_primitive.curve->getType() != CurvePrimitive::CLOTHOID)
        {
            _primitive.curve->setParams(x);
            return;
        }

        VectorXd xm = x;
        xm(CurvePrimitive::DCURVATURE) = (xm(CurvePrimitive::DCURVATURE) - xm(CurvePrimitive::CURVATURE)) / xm(CurvePrimitive::LENGTH);
        _primitive.curve->setParams(xm);
    }

private:
    FitPrimitive _primitive;
    ErrorComputerConstPtr _errorComputer;
};

class DefaultPrimitiveFitter : public Algorithm<PRIMITIVE_FITTING>
{
public:
    DefaultPrimitiveFitter(bool adjust) : _adjust(adjust) {}

    string name() const { return _adjust ? "Adjust" : "Default"; }

private:
    bool _adjust;

protected:

    void _run(const Fitter &fitter, AlgorithmOutput<PRIMITIVE_FITTING> &out)
    {
        const VectorC<bool> &corners = fitter.output<RESAMPLING>()->corners;
        PolylineConstPtr poly = fitter.output<RESAMPLING>()->output;
        ErrorComputerConstPtr errorComputer = fitter.output<ERROR_COMPUTER>()->errorComputer;
        smart_ptr<const AlgorithmOutput<OVERSKETCHING> > osOutput = fitter.output<OVERSKETCHING>();

        const VectorC<Vector2d> &pts = poly->pts();

        const double errorThreshold = fitter.scaledParameter(Parameters::ERROR_THRESHOLD);
        std::string typeNames[3] = { "Lines", "Arcs", "Clothoids" };
        bool inflectionAccounting = fitter.params().get(Parameters::INFLECTION_COST) > 0.;

        if(osOutput->startCurve)
        {
            FitPrimitive fit;
            fit.curve = osOutput->startCurve->clone();
            fit.startIdx = -1;
            fit.endIdx = 0;
            fit.numPts = 1;
            fit.error = 0.;
            fit.fixed = true;
            fit.startCurvSign = (fit.curve->startCurvature() >= 0) ? 1 : -1;
            fit.endCurvSign = (fit.curve->endCurvature() >= 0) ? 1 : -1;
            out.primitives.push_back(fit);

            for(int i = 1; i < pts.size(); ++i)
            {
                fit.curve = fit.curve->clone();
                fit.endIdx = i;

                //extend the curve up to the next point
                fit.curve->trim(0, fit.curve->length() + 2 * (pts[i] - pts[i - 1]).norm());
                fit.curve->trim(0, fit.curve->project(pts[i]));

                fit.endCurvSign = (fit.curve->endCurvature() >= 0) ? 1 : -1;
                fit.error = errorComputer->computeErrorForCost(fit.curve, 0, fit.endIdx, false);

                fit.numPts++;

                if(fit.error > errorThreshold * errorThreshold)
                    break;

                //Debugging::get()->drawPrimitive(fit.curve, "Start Curves",  0);

                out.primitives.push_back(fit);

                if(corners[i])
                    break; //don't pull curve past corner
            }
        }

        if(osOutput->endCurve)
        {
            FitPrimitive fit;
            fit.curve = osOutput->endCurve->clone();
            fit.startIdx = (int)pts.size() - 1;
            fit.endIdx = fit.startIdx + 1;
            fit.numPts = 1;
            fit.error = 0.;
            fit.fixed = true;
            fit.startCurvSign = (fit.curve->startCurvature() >= 0) ? 1 : -1;
            fit.endCurvSign = (fit.curve->endCurvature() >= 0) ? 1 : -1;
            out.primitives.push_back(fit);

            for(int i = fit.startIdx - 1; i >= 0; --i)
            {
                fit.curve = fit.curve->clone();
                fit.startIdx = i;

                //extend the curve up to the previous point
                fit.curve->trim(-2 * (pts[i] - pts[i + 1]).norm(), fit.curve->length());
                fit.curve->trim(fit.curve->project(pts[i]), fit.curve->length());

                fit.startCurvSign = (fit.curve->startCurvature() >= 0) ? 1 : -1;
                fit.error = errorComputer->computeErrorForCost(fit.curve, fit.startIdx, (int)pts.size() - 1, true, false);

                fit.numPts++;

                if(fit.error > errorThreshold * errorThreshold)
                    break;

                //Debugging::get()->drawPrimitive(fit.curve, "End Curves",  0);

                out.primitives.push_back(fit);

                if(corners[i])
                    break; //don't pull curve past corner
            }
        }

        for(int i = 0; i < pts.size(); ++i) //iterate over start points
        {
            FitterBasePtr fitters[3];
            fitters[0] = new LineFitter();
            fitters[1] = new ArcFitter();
            fitters[2] = new ClothoidFitter();

            for(int type = 0; type <= 2; ++type) //iterate over lines, arcs, clothoids
            {
                int fitSoFar = 0;

                bool needType = fitter.params().get(Parameters::ParameterType(Parameters::LINE_COST + type)) < Parameters::infinity;

                for(VectorC<Vector2d>::Circulator circ = pts.circulator(i); !circ.done(); ++circ)
                {
                    ++fitSoFar;

                    if(!needType && (type == 2 || fitSoFar >= 3 + type)) //if we don't need primitives of this type
                        break;

                    fitters[type]->addPoint(*circ);
                    if(fitSoFar >= 2 + type) //at least two points per line, etc.
                    {
                        CurvePrimitivePtr curve = fitters[type]->getPrimitive();
                        Vector3d color(0, 0, 0);
                        color[type] = 1;

                        FitPrimitive fit;
                        fit.curve = curve;
                        fit.startIdx = i;
                        fit.endIdx = circ.index();
                        fit.numPts = fitSoFar;
                        fit.startCurvSign = (curve->startCurvature() >= 0) ? 1 : -1;
                        fit.endCurvSign = (curve->endCurvature() >= 0) ? 1 : -1;

                        if(_adjust)
                            adjustPrimitive(fit, fitter);

                        fit.error = errorComputer->computeErrorForCost(curve, i, fit.endIdx);

                        double length = poly->lengthFromTo(i, fit.endIdx);
                        if(fit.error > errorThreshold * errorThreshold)
                            break;

                        //Debugging::get()->drawCurve(curve, color, typeNames[type]);
                        out.primitives.push_back(fit);

                        if(type == 0 && inflectionAccounting) //line with "opposite" curvature
                        {
                            fit.startCurvSign = -fit.startCurvSign;
                            fit.endCurvSign = -fit.endCurvSign;
                            out.primitives.push_back(fit);
                        }

                        //if different start and end curvatures
                        if(fit.startCurvSign != fit.endCurvSign && inflectionAccounting)
                        {
                            double start = poly->idxToParam(i);
                            double end = poly->idxToParam(fit.endIdx);
                            CurvePrimitivePtr startNoCurv = static_pointer_cast<ClothoidFitter>(fitters[2])->getCurveWithZeroCurvature(0);
                            CurvePrimitivePtr endNoCurv = static_pointer_cast<ClothoidFitter>(fitters[2])->getCurveWithZeroCurvature(end - start);

                            fit.curve = startNoCurv;
                            fit.startCurvSign = fit.endCurvSign = (startNoCurv->endCurvature() > 0. ? 1 : -1);

                            if(_adjust)
                                adjustPrimitive(fit, fitter);

                            fit.error = errorComputer->computeErrorForCost(fit.curve, i, fit.endIdx);

                            if(fit.error < errorThreshold * errorThreshold)
                            {
                                out.primitives.push_back(fit);
                                //Debugging::get()->drawCurve(fit.curve, color, typeNames[type]);
                            }

                            fit.curve = endNoCurv;
                            fit.startCurvSign = fit.endCurvSign = (endNoCurv->startCurvature() > 0. ? 1 : -1);

                            if(_adjust)
                                adjustPrimitive(fit, fitter);

                            fit.error = errorComputer->computeErrorForCost(fit.curve, i, fit.endIdx);

                            if(fit.error < errorThreshold * errorThreshold)
                            {
                                out.primitives.push_back(fit);
                                //Debugging::get()->drawCurve(fit.curve, color, typeNames[type]);
                            }
                        }
                    }
                    if(fitSoFar > 1 && corners[circ.index()])
                        break;
                }
            }
        }
    }

    void adjustPrimitive(const FitPrimitive &primitive, const Fitter &fitter)
    {
        ErrorComputerConstPtr errorComputer = fitter.output<ERROR_COMPUTER>()->errorComputer;
        bool inflectionAccounting = fitter.params().get(Parameters::INFLECTION_COST) > 0.;

        vector<LSBoxConstraint> constraints;

        //minimum length constraint
        constraints.push_back(LSBoxConstraint(CurvePrimitive::LENGTH, primitive.curve->length() * 0.5, 1));

        //curvature sign constraints
        if(inflectionAccounting)
        {
            if(primitive.curve->getType() >= CurvePrimitive::ARC)
                constraints.push_back(LSBoxConstraint(CurvePrimitive::CURVATURE, 0., primitive.startCurvSign));
            if(primitive.curve->getType() == CurvePrimitive::CLOTHOID)
                constraints.push_back(LSBoxConstraint(CurvePrimitive::DCURVATURE, 0., primitive.endCurvSign));
        }

        //solve
        OneCurveProblem problem(primitive, errorComputer);
        LSSolver solver(&problem, constraints);
        solver.setDefaultDamping(fitter.params().get(Parameters::CURVE_ADJUST_DAMPING));
        solver.setMaxIter(1);
        problem.setParams(solver.solve(problem.params()));
    }
};

void Algorithm<PRIMITIVE_FITTING>::_initialize()
{
    new DefaultPrimitiveFitter(false);
    new DefaultPrimitiveFitter(true);
}

END_NAMESPACE_Cornu


