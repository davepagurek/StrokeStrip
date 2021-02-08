/*--
    TwoCurveCombine.cpp  

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

#include "TwoCurveCombine.h"
#include "Resampler.h"
#include "PrimitiveFitter.h"
#include "CurvePrimitive.h"
#include "Fitter.h"
#include "Polyline.h"
#include "ErrorComputer.h"
#include "AngleUtils.h"
#include "Solver.h"
#include <map>
#include <cstdio> //TMP
#include <iostream> //TMP

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

class CombinedCurve
{
public:
    CombinedCurve(FitPrimitive p[2], int continuity, const Fitter &fitter)
        : _continuity(continuity), _evalCount(0)
    {
        _errorComputer = fitter.output<ERROR_COMPUTER>()->errorComputer;
        int sampledPts = fitter.output<RESAMPLING>()->output->pts().size();
        double adjustmentPoint = fitter.params().get(Parameters::TWO_CURVE_CURVATURE_ADJUST);

        CurvePrimitive::ParamVec v[2];
        for(int i = 0; i < 2; ++i)
        {
            _c[i] = p[i].curve;
            v[i] = _c[i]->params();
            _type[i] = _c[i]->getType();
            _from[i] = p[i].startIdx;
            _to[i] = p[i].endIdx;
            _angleWeight[i] = _c[i]->length();
            _curvatureWeight[i] = _c[i]->length() * _c[i]->length();
            _origEndAngle[i] = _c[i]->endAngle();
            _origEndCurvature[i] = _c[i]->endCurvature();
        }

        //adjust the curves to enforce continuity
        Vector2d midPt;

        double w0 = _c[0]->length() / (_c[0]->length() + _c[1]->length());
        double w1 = _c[1]->length() / (_c[0]->length() + _c[1]->length());

        if(p[0].isFixed())
        {
            w0 = 1;
            w1 = 0;
            midPt = _c[0]->startPos();
        }
        else if(p[1].isFixed())
        {
            w0 = 0;
            w1 = 1;
            midPt = _c[1]->startPos();
        }
        else
        {
            midPt = 0.5 * (_c[0]->startPos() + _c[1]->startPos());
        }

        for(int i = 0; i < 2; ++i)
        {
            v[i][CurvePrimitive::X] = midPt[0];
            v[i][CurvePrimitive::Y] = midPt[1];
        }

        if(_continuity >= 1)
        {
            double angle1 = _c[1]->startAngle();
            double avgAngle = (w1 * angle1 + w0 * AngleUtils::toRange(_c[0]->startAngle() + PI, angle1 - PI));
            for(int i = 0; i < 2; ++i)
                v[i][CurvePrimitive::ANGLE] = avgAngle + (i ? 0. : PI);
        }

        if(_continuity == 2)
        {
            if(_type[0] == CurvePrimitive::LINE)
            {
                //keep the end curvature
                v[1][CurvePrimitive::DCURVATURE] += adjustmentPoint * v[1][CurvePrimitive::CURVATURE] / v[1][CurvePrimitive::LENGTH];
                v[1][CurvePrimitive::CURVATURE] = 0.;
            }
            else if(_type[1] == CurvePrimitive::LINE)
            {
                //keep the end curvature
                v[0][CurvePrimitive::DCURVATURE] += adjustmentPoint * v[0][CurvePrimitive::CURVATURE] / v[0][CurvePrimitive::LENGTH];
                v[0][CurvePrimitive::CURVATURE] = 0.;
            }
            else
            {
                double avgCurv = w1 * _c[1]->startCurvature() - w0 * _c[0]->startCurvature();
                //keep the end curvatures
                if(_type[0] == CurvePrimitive::CLOTHOID)
                    v[0][CurvePrimitive::DCURVATURE] += adjustmentPoint * (v[0][CurvePrimitive::CURVATURE] + avgCurv) / v[0][CurvePrimitive::LENGTH];
                if(_type[1] == CurvePrimitive::CLOTHOID)
                    v[1][CurvePrimitive::DCURVATURE] += adjustmentPoint * (v[1][CurvePrimitive::CURVATURE] - avgCurv) / v[1][CurvePrimitive::LENGTH];
                v[1][CurvePrimitive::CURVATURE] = avgCurv;
                v[0][CurvePrimitive::CURVATURE] = -avgCurv;
            }
        }

        _c[0]->setParams(v[0]);
        _c[1]->setParams(v[1]);

        //figure out the ranges over which error is computed --they should not overlap
        if(_continuity == 1)
            swap(_to[0], _from[1]); //they overlap by 1 originally
        if(_continuity == 2) //they overlap by 2
        {
            _from[1] = (_from[1] + 2) % sampledPts;
            _to[0] = (_to[0] + sampledPts - 2) % sampledPts;
        }

        initMapping();
    }

    void initMapping()
    {
        int paramIdx = 0;

        _mapping.reserve(8 + _type[0] + _type[1]);
        _mapping.push_back(MappingElement(paramIdx, 0, CurvePrimitive::X));
        _mapping.push_back(MappingElement(paramIdx++, 1, CurvePrimitive::X));
        _mapping.push_back(MappingElement(paramIdx, 0, CurvePrimitive::Y));
        _mapping.push_back(MappingElement(paramIdx++, 1, CurvePrimitive::Y));
        _mapping.push_back(MappingElement(paramIdx++, 0, CurvePrimitive::LENGTH));
        _mapping.push_back(MappingElement(paramIdx++, 1, CurvePrimitive::LENGTH));
        _mapping.push_back(MappingElement(paramIdx, 1, CurvePrimitive::ANGLE));
        if(_continuity == 0)
        {
            ++paramIdx;
            _mapping.push_back(MappingElement(paramIdx++, 0, CurvePrimitive::ANGLE));
        }
        else
        {
            _mapping.push_back(MappingElement(paramIdx++, 0, CurvePrimitive::ANGLE, MappingElement::PlusPi));
        }

        //now the nontrivial part -- curvature
        if(_continuity < 2) //if there's G0 or G1 continuity, we just have separate curvatures
        {
            if(_type[0] >= CurvePrimitive::ARC)
                _mapping.push_back(MappingElement(paramIdx++, 0, CurvePrimitive::CURVATURE));
            if(_type[1] >= CurvePrimitive::ARC)
                _mapping.push_back(MappingElement(paramIdx++, 1, CurvePrimitive::CURVATURE));
        }
        else //curvature continuous
        {
            if(_type[0] == CurvePrimitive::LINE)
                _mapping.push_back(MappingElement(paramIdx, 1, CurvePrimitive::CURVATURE, MappingElement::Zero));
            else if(_type[1] == CurvePrimitive::LINE)
                _mapping.push_back(MappingElement(paramIdx, 0, CurvePrimitive::CURVATURE, MappingElement::Zero));
            else //both are at least arcs
            {
                _mapping.push_back(MappingElement(paramIdx, 1, CurvePrimitive::CURVATURE));
                _mapping.push_back(MappingElement(paramIdx++, 0, CurvePrimitive::CURVATURE, MappingElement::Negative));
            }
        }

        if(_type[0] == CurvePrimitive::CLOTHOID)
            _mapping.push_back(MappingElement(paramIdx++, 0, CurvePrimitive::DCURVATURE, MappingElement::EndCurvature));
        if(_type[1] == CurvePrimitive::CLOTHOID)
            _mapping.push_back(MappingElement(paramIdx++, 1, CurvePrimitive::DCURVATURE, MappingElement::EndCurvature));
    }

    int numParams() const
    {
        return 6 + _type[0] + _type[1] - _continuity; //6 is the DOFs of two lines connected at an endpoint
    }

    void setParams(const VectorXd &params)
    {
        CurvePrimitive::ParamVec v[2];
        v[0].resize(_c[0]->numParams());
        v[1].resize(_c[1]->numParams());

        for(int i = 0; i < (int)_mapping.size(); ++i)
        {
            const MappingElement &elem = _mapping[i];
            double val;
            switch(elem.type)
            {
            case MappingElement::Normal:       val = params[elem.paramIdx];       break;
            case MappingElement::Negative:     val = -params[elem.paramIdx];      break;
            case MappingElement::PlusPi:       val = PI + params[elem.paramIdx];  break;
            case MappingElement::Zero:         val = 0;                           break;
            case MappingElement::EndCurvature:
                //by this point v[elem.curveIdx][CurvePrimitive::CURVATURE]) should already be set
                val = (params[elem.paramIdx] - v[elem.curveIdx][CurvePrimitive::CURVATURE]) / v[elem.curveIdx][CurvePrimitive::LENGTH];
                break;
            }
            v[elem.curveIdx][elem.parameter] = val;
        }

        _c[0]->setParams(v[0]);
        _c[1]->setParams(v[1]);
    }

    void getParams(VectorXd &params) const
    {
        CurvePrimitive::ParamVec v[2];
        v[0] = _c[0]->params();
        v[1] = _c[1]->params();
        params.resize(numParams());

        for(int i = 0; i < (int)_mapping.size(); ++i)
        {
            const MappingElement &elem = _mapping[i];
            if(elem.type == MappingElement::Normal)
                params[elem.paramIdx] = v[elem.curveIdx][elem.parameter];
            else if(elem.type == MappingElement::EndCurvature)
                params[elem.paramIdx] = _c[elem.curveIdx]->endCurvature();
        }
    }

    double computeErrorForCost(int idx) const
    {
        double endCost = SQR(_angleWeight[idx] * AngleUtils::toRange(_c[idx]->endAngle() - _origEndAngle[idx], -PI));
        endCost += SQR(_curvatureWeight[idx] * (_c[idx]->endCurvature() - _origEndCurvature[idx]));
        bool firstToEndpoint = (idx == 0 || _continuity == 0);
        bool lastToEndpoint = (idx == 1 || _continuity == 0);
        bool reverse = (idx == 0);
        return endCost + _errorComputer->computeErrorForCost(_c[idx], _from[idx], _to[idx], firstToEndpoint, lastToEndpoint, reverse);
    }

    double computeError() const
    {
        return _errorComputer->computeError(_c[0], _from[0], _to[0], true, _continuity == 0, true) +
               _errorComputer->computeError(_c[1], _from[1], _to[1], _continuity == 0, true);
    }

    void computeErrorVector(VectorXd &outError, MatrixXd &outErrorDer) const
    {
#if 0
        char name[100];
        sprintf(name, "Curves %d", const_cast<CombinedCurve *>(this)->_evalCount++);
        Debugging::get()->drawCurve(_c[0], Vector3d(1, 0, 0), name);
        Debugging::get()->drawCurve(_c[1], Vector3d(0, 0, 1), name);
#endif
        VectorXd err[2];
        MatrixXd errDer[2];

        //perhaps it should be whether that point is a corner, rather than continuity
        _errorComputer->computeErrorVector(_c[0], _from[0], _to[0], err[0], errDer + 0, true, _continuity == 0, true);
        _errorComputer->computeErrorVector(_c[1], _from[1], _to[1], err[1], errDer + 1, _continuity == 0, true);

#if 1
        CurvePrimitive::EndDer endDer;
        for(int c = 0; c < 2; ++c)
        {
            size_t offs = err[c].size();
            err[c].conservativeResize(offs + 2);
            errDer[c].conservativeResize(offs + 2, NoChange);

            err[c][offs] = _angleWeight[c] * AngleUtils::toRange(_c[c]->endAngle() - _origEndAngle[c], -PI);
            err[c][offs + 1] = _curvatureWeight[c] * (_c[c]->endCurvature() - _origEndCurvature[c]);
            
            _c[c]->derivativeAtEnd(2, endDer);
            endDer.row(2) *= _angleWeight[c];
            endDer.row(3) *= _curvatureWeight[c];
            errDer[c].block(offs, 0, 2, errDer[c].cols()) = endDer.block(2, 0, 2, errDer[c].cols());
        }
#endif

        _c[0]->toEndCurvatureDerivative(errDer[0]);
        _c[1]->toEndCurvatureDerivative(errDer[1]);

        outError.resize(err[0].size() + err[1].size());
        outError.segment(0, err[0].size()) = err[0];
        outError.segment(err[0].size(), err[1].size()) = err[1];

        outErrorDer = MatrixXd::Zero(outError.size(), numParams());

        for(int i = 0; i < (int)_mapping.size(); ++i)
        {
            const MappingElement &elem = _mapping[i];
            MatrixXd::Index offset = elem.curveIdx ? err[0].size() : 0;
            MatrixXd::Index colSize = err[elem.curveIdx].size();

            switch(elem.type)
            {
            case MappingElement::Normal:    
            case MappingElement::PlusPi:
            case MappingElement::EndCurvature:
                outErrorDer.col(elem.paramIdx).segment(offset, colSize) += errDer[elem.curveIdx].col(elem.parameter);
                break;
            case MappingElement::Negative:
                outErrorDer.col(elem.paramIdx).segment(offset, colSize) -= errDer[elem.curveIdx].col(elem.parameter);
                break;
            case MappingElement::Zero:
                break;
            }
        }
    }

    CurvePrimitivePtr getCurve(int idx) const
    {
        return _c[idx];
    }

    int getParamIndex(int curveIdx, CurvePrimitive::Param parameter) const
    {
        int mappingIdx = _getMappingIndex(curveIdx, parameter);
        if(mappingIdx >= 0)
            return _mapping[mappingIdx].paramIdx;
        return -1;
    }

    int getNormalParamIndex(int curveIdx, CurvePrimitive::Param parameter) const
    {
        int mappingIdx = _getMappingIndex(curveIdx, parameter);
        if(mappingIdx >= 0 && (_mapping[mappingIdx].type == MappingElement::Normal || _mapping[mappingIdx].type == MappingElement::EndCurvature))
            return _mapping[mappingIdx].paramIdx;
        return -1;
    }

private:
    int _getMappingIndex(int curveIdx, CurvePrimitive::Param parameter) const
    {
        for(int i = 0; i < (int)_mapping.size(); ++i)
            if(_mapping[i].parameter == parameter && _mapping[i].curveIdx == curveIdx)
                return i;
        return -1; //should not get here
    }

    struct MappingElement
    {
        enum MappingType
        {
            Normal, Negative, PlusPi, Zero, EndCurvature
        };

        MappingElement(int pIdx, int cIdx, CurvePrimitive::Param param, MappingType inType = Normal)
            : paramIdx(pIdx), curveIdx(cIdx), parameter(param), type(inType) {}

        int paramIdx;
        int curveIdx; //0 or 1
        CurvePrimitive::Param parameter;

        MappingType type;
    };

    vector<MappingElement> _mapping;
    CurvePrimitivePtr _c[2];
    int _from[2], _to[2];
    CurvePrimitive::PrimitiveType _type[2];
    ErrorComputerConstPtr _errorComputer;
    int _continuity;
    int _evalCount;
    double _angleWeight[2];
    double _curvatureWeight[2];
    double _origEndAngle[2];
    double _origEndCurvature[2];
};

class TwoCurveProblem : public LSProblem
{
public:
    TwoCurveProblem(CombinedCurve &curves) : _curves(curves) {}

    //overrides
    double error(const Eigen::VectorXd &x, LSEvalData *)
    {
        _curves.setParams(x);
        return _curves.computeError();
    }
    LSEvalData *createEvalData()
    {
        return new LSDenseEvalData();
    }
    void eval(const Eigen::VectorXd &x, LSEvalData *data)
    {
        LSDenseEvalData *twoCurveData = static_cast<LSDenseEvalData *>(data);
        _curves.setParams(x);
        _curves.computeErrorVector(twoCurveData->errVectorRef(), twoCurveData->errDerRef());
        //cout << "Err = " << twoCurveData->error() << endl;
    }

private:
    CombinedCurve &_curves;
};

Combination twoCurveCombine(int p1, int p2, int continuity, const Fitter &fitter)
{
    const vector<FitPrimitive> &primitives = fitter.output<PRIMITIVE_FITTING>()->primitives;
    const VectorC<Vector2d> &pts = fitter.output<RESAMPLING>()->output->pts();
    ErrorComputerConstPtr errorComputer = fitter.output<ERROR_COMPUTER>()->errorComputer;

    Combination out;

    out.c1 = primitives[p1].curve->clone();
    out.c2 = primitives[p2].curve->clone();

    //trim c1 and c2
    Vector2d trimPt = 0.5 * (out.c1->endPos() + out.c2->startPos());
    if(!primitives[p1].isFixed())
        out.c1->trim(0, out.c1->project(trimPt));
    if(!primitives[p2].isFixed())
        out.c2->trim(out.c2->project(trimPt), out.c2->length());

    out.c1->flip();

    FitPrimitive primArray[2] = { primitives[p1], primitives[p2] };
    primArray[0].curve = out.c1;
    primArray[1].curve = out.c2;

    CombinedCurve combined(primArray, continuity, fitter);

    VectorXd x;
    combined.getParams(x);
    VectorXd errVec;
    MatrixXd errVecDer;

    vector<LSBoxConstraint> constraints;

    //see if any curve is a start curve or end curve and therefore should not be allowed to move
    for(int curveIdx = 0; curveIdx < 2; ++curveIdx)
    {
        if(primArray[curveIdx].isFixed())
        {
            for(int i = 0; i < primArray[curveIdx].curve->numParams(); ++i)
            {
                int paramIdx = combined.getParamIndex(curveIdx, CurvePrimitive::Param(i));
                constraints.push_back(LSBoxConstraint(paramIdx, x[paramIdx], 0));
            }
        }
    }

#if 0
    static int cnt = 0;
    bool origDrawn = !constraints.empty() && !(cnt++);
    if(origDrawn)
    {
        Debugging::get()->drawCurve(primitives[p1].curve, Vector3d(1, 0, 0), "Curves Orig");
        Debugging::get()->drawCurve(primitives[p2].curve, Vector3d(0, 0, 1), "Curves Orig");
    }
#endif

    //minimum length
    if(constraints.empty()) //if neither curve is fixed
    {
        for(int curveIdx = 0; curveIdx < 2; ++curveIdx)
            constraints.push_back(LSBoxConstraint(combined.getParamIndex(curveIdx, CurvePrimitive::LENGTH),
                                                  combined.getCurve(curveIdx)->length() * 0.5, 1));

        bool inflectionAccounting = fitter.params().get(Parameters::INFLECTION_COST) > 0.;

        if(inflectionAccounting)
        {
            bool constantCurvature = primArray[0].startCurvSign == primArray[0].endCurvSign &&
                                     primArray[0].startCurvSign == primArray[1].startCurvSign &&
                                     primArray[0].startCurvSign == primArray[1].endCurvSign;

            for(int curveIdx = 0; curveIdx < 2; ++curveIdx)
            {
                int curStartSign = curveIdx == 0 ? -primArray[0].endCurvSign : primArray[1].startCurvSign;
                int curEndSign = curveIdx == 0 ? -primArray[0].startCurvSign : primArray[1].endCurvSign;

                if(primArray[curveIdx].curve->getType() == CurvePrimitive::CLOTHOID)
                {
                    //cout << "Idx = " << curveIdx << " end constraint = " << curEndSign << endl;
                    constraints.push_back(LSBoxConstraint(combined.getParamIndex(curveIdx, CurvePrimitive::DCURVATURE), 0., curEndSign));
                }
                int startCurv = combined.getNormalParamIndex(curveIdx, CurvePrimitive::CURVATURE);
                if(continuity == 2 && !constantCurvature)
                    startCurv = -1; //don't add a curvature constraint
                if(startCurv >= 0)
                {
                    //cout << "Idx = " << curveIdx << " start constraint = " << curStartSign  << endl;
                    constraints.push_back(LSBoxConstraint(startCurv, 0., curStartSign));
                }
            }
        }
    }

    TwoCurveProblem problem(combined);
    LSSolver solver(&problem, constraints);
    solver.setDefaultDamping(fitter.params().get(Parameters::CURVE_ADJUST_DAMPING));
    solver.setMaxIter(5);
    //solver.verifyDerivatives(x);
	bool is_valid;
    x = solver.solve(x, &is_valid);
	if (is_valid) {
		combined.setParams(x);
		out.err1 = combined.computeErrorForCost(0);
		out.err2 = combined.computeErrorForCost(1);
	}
	else {
		out.err1 = std::numeric_limits<double>::infinity();
		out.err2 = std::numeric_limits<double>::infinity();
		//cout << "Err = " << (out.err1 + out.err2) << endl;
	}

#if 0
    if(origDrawn)
    {
        Debugging::get()->drawCurve(combined.getCurve(0), Vector3d(1, 0, 0), "Curves Final");
        Debugging::get()->drawCurve(combined.getCurve(1), Vector3d(0, 0, 1), "Curves Final");
    }
#endif

    return out;
}

END_NAMESPACE_Cornu


