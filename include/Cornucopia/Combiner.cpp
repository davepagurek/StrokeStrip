/*--
    Combiner.cpp  

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

#include "Combiner.h"
#include "GraphConstructor.h"
#include "PrimitiveFitter.h"
#include "PrimitiveSequence.h"
#include "PathFinder.h"
#include "Preprocessing.h"
#include "Resampler.h"
#include "Polyline.h"
#include "Solver.h"
#include "ErrorComputer.h"
#include "Fitter.h"
#include "Oversketcher.h"
#include "PiecewiseLinearUtils.h"

#include <iterator>
#include <cstdio>

#include <Eigen/LU>
#include <Eigen/Cholesky>

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

#define SPARSE 1

class MulticurveDenseEvalData : public LSEvalData
{
public:
    //overrides
    double error() const { return _con.squaredNorm(); }

    void solveForDelta(double damping, Eigen::VectorXd &out, std::set<LSBoxConstraint> &constraints)
    {
        size_t vars = _errDer.cols();

        size_t size = vars + _con.size() + constraints.size();
        MatrixXd lhs = MatrixXd::Zero(size, size);
        VectorXd rhs = VectorXd::Zero(size);

        lhs.block(0, 0, vars, vars) = _errDer.transpose() * _errDer;
        rhs.segment(0, vars) = -_errDer.transpose() * _err;

        lhs.block(vars, 0, _conDer.rows(), _conDer.cols()) = _conDer;
        lhs.block(0, vars, _conDer.cols(), _conDer.rows()) = _conDer.transpose();
        rhs.segment(vars, _con.size()) = -_con;

        int cnt = 0;
        for(set<LSBoxConstraint>::const_iterator it = constraints.begin(); it != constraints.end(); ++it, ++cnt)
            lhs(vars + _conDer.rows() + cnt, it->index) = lhs(it->index, vars + _conDer.rows() + cnt) = 1.;

        //lhs += damping * MatrixXd::Identity(size, size);
        lhs.block(0, 0, vars, vars) += damping * MatrixXd::Identity(vars, vars);
        VectorXd result = lhs.lu().solve(rhs);

        out = result.segment(0, vars);

#if 0
        printf("Solve err = %lf\n", (lhs * result - rhs).norm());
        if(_conDer.size() > 0)
            printf("Con Solve err = %lf\n", (_conDer * out + _con).norm());
#endif

        //check which constraints we don't need
        cnt = 0;
        for(set<LSBoxConstraint>::iterator it = constraints.begin(); it != constraints.end(); ++cnt)
        {
            set<LSBoxConstraint>::iterator next = it;
            ++next;
            if(result(vars + _conDer.rows() + cnt) * it->sign > 0)
            {
                //printf("Unsetting constraint on variable at index %d\n", it->index);
                constraints.erase(it);
            }
            it = next;
        }
    }

    //for derivative verification, combine error and constraints
    VectorXd errVec() const
    {
        VectorXd out(_err.size() + _con.size());
        out << _err, _con;
        return out;
    }

    MatrixXd errVecDer() const
    {
        MatrixXd out(_err.size() + _con.size(), _errDer.cols());
        out << _errDer, _conDer;
        return out;
    }

    VectorXd &errVectorRef() { return _err; }
    MatrixXd &errDerRef() { return _errDer; }

    VectorXd &conVectorRef() { return _con; }
    MatrixXd &conDerRef() { return _conDer; }

private:
    Eigen::VectorXd _err;
    Eigen::MatrixXd _errDer;
    Eigen::VectorXd _con;
    Eigen::MatrixXd _conDer;
};

class MulticurveSparseEvalData : public LSEvalData
{
public:
    typedef Matrix<double, Dynamic, Dynamic, 0, 6, 6> BlockType;
    typedef vector<BlockType, aligned_allocator<BlockType> > BlockVectorType;
    typedef LLT<BlockType> BlockCholType;
    typedef vector<BlockCholType, aligned_allocator<BlockCholType> > BlockCholVectorType;

    //overrides
    double error() const { return _con.squaredNorm(); }

    void solveForDelta(double damping, Eigen::VectorXd &out, std::set<LSBoxConstraint> &constraints)
    {
        _computeIndices();

        size_t vars = _blockIndices.back();
        size_t cons = _con.size() + constraints.size();

        size_t size = vars + cons;
        VectorXd rhs = VectorXd::Zero(size);
        rhs.segment(0, vars) = _err;
        rhs.segment(vars, _con.size()) = -_con;

        BlockCholVectorType cholBlocks(_errDerBlocks.size());
        for(size_t i = 0; i < cholBlocks.size(); ++i)
            cholBlocks[i] = BlockCholType(_errDerBlocks[i] + damping * MatrixXd::Identity(_blockSizes[i], _blockSizes[i]));

        MatrixXd C = MatrixXd::Zero(cons, vars);
        C.block(0, 0, _con.size(), vars) = _conDer;

        int cnt = 0;
        for(set<LSBoxConstraint>::const_iterator it = constraints.begin(); it != constraints.end(); ++it, ++cnt)
            C(_conDer.rows() + cnt, it->index) = 1.;

        MatrixXd K(C.cols(), C.rows());
        _solveCholL(K, cholBlocks, C.transpose());
        LLT<MatrixXd> CholK(K.transpose() * K);

        VectorXd mid(size), result(size);

        //solve for mid
        _solveCholL(mid.segment(0, vars), cholBlocks, rhs.segment(0, vars));
        mid.segment(vars, cons) = CholK.matrixL().solve(rhs.segment(vars, cons) - K.adjoint() * mid.segment(0, vars));

        //solve for result
        result.segment(vars, cons) = -CholK.matrixU().solve(mid.segment(vars, cons));
        _solveCholU(result.segment(0, vars), cholBlocks, mid.segment(0, vars) - K * result.segment(vars, cons));

        out = result.segment(0, vars);
#if 0
        if(_conDer.size() > 0)
            printf("Con Solve err = %lf\n", (_conDer * out + _con).norm());
#endif

        //check which constraints we don't need
        cnt = 0;
        for(set<LSBoxConstraint>::iterator it = constraints.begin(); it != constraints.end(); ++cnt)
        {
            set<LSBoxConstraint>::iterator next = it;
            ++next;
            if(result(vars + _conDer.rows() + cnt) * it->sign > 0)
            {
                //printf("Unsetting constraint on variable at index %d\n", it->index);
                constraints.erase(it);
            }
            it = next;
        }
    }

    VectorXd &errVectorRef() { return _err; }
    BlockVectorType &errDerBlocksRef() { return _errDerBlocks; }

    VectorXd &conVectorRef() { return _con; }
    MatrixXd &conDerRef() { return _conDer; }

private:
    void _computeIndices()
    {
        _blockIndices.resize(_errDerBlocks.size() + 1);
        _blockSizes.resize(_errDerBlocks.size());
        _blockIndices[0] = 0;

        for(int i = 0; i < (int)_errDerBlocks.size(); ++i)
        {
            _blockSizes[i] = _errDerBlocks[i].rows();
            _blockIndices[i + 1] = _blockIndices[i] + _blockSizes[i];
        }
    }

    MatrixXd _blocksToMatrix(const BlockVectorType &blocks)
    {
        size_t sz = _blockIndices.back();
        MatrixXd out = MatrixXd::Zero(sz, sz);

        for(int i = 0; i < (int)blocks.size(); ++i)
            out.block(_blockIndices[i], _blockIndices[i], _blockSizes[i], _blockSizes[i]) = blocks[i];

        return out;
    }

    template<class M1, class M2> void _solveCholL(const M1 &out, const BlockCholVectorType &chols, const M2 &rhs)
    {
        for(int i = 0; i < (int)chols.size(); ++i)
            const_cast<M1 &>(out).block(_blockIndices[i], 0, _blockSizes[i], out.cols()) =
                chols[i].matrixL().solve(rhs.block(_blockIndices[i], 0, _blockSizes[i], rhs.cols()));
    }

    template<class M1, class M2> void _solveCholU(const M1 &out, const BlockCholVectorType &chols, const M2 &rhs)
    {
        for(int i = 0; i < (int)chols.size(); ++i)
            const_cast<M1 &>(out).block(_blockIndices[i], 0, _blockSizes[i], out.cols()) =
                chols[i].matrixU().solve(rhs.block(_blockIndices[i], 0, _blockSizes[i], rhs.cols()));
    }

    vector<size_t> _blockIndices, _blockSizes;
    BlockVectorType _errDerBlocks;

    Eigen::VectorXd _err;
    Eigen::VectorXd _con;
    Eigen::MatrixXd _conDer;
};

class MulticurveProblem : public LSProblem
{
public:
    MulticurveProblem(const Fitter &fitter)
        : _primitives(fitter.output<PRIMITIVE_FITTING>()->primitives), _iter(0)
    {
        smart_ptr<const AlgorithmOutput<GRAPH_CONSTRUCTION> > graph = fitter.output<GRAPH_CONSTRUCTION>();
        const vector<int> &path = fitter.output<PATH_FINDING>()->path;
        _errorComputer = fitter.output<ERROR_COMPUTER>()->errorComputer;
        _closed = fitter.output<CURVE_CLOSING>()->closed;
        _inflectionAccounting = fitter.params().get(Parameters::INFLECTION_COST) > 0.;

        for(int i = 0; i < (int)path.size(); ++i)
        {
            _primIdcs.push_back(graph->edges[path[i]].startVtx);
            _continuities.push_back(graph->edges[path[i]].continuity);
        }
        if(!_closed)
            _primIdcs.push_back(graph->edges[path.back()].endVtx);
        
        _curves = VectorC<CurvePrimitivePtr>((int)_primIdcs.size(), _closed ? CIRCULAR : NOT_CIRCULAR);
        _curveRanges = VectorC<pair<int, int> >((int)_primIdcs.size(), _curves.circular());
        for(int i = 0; i < (int)_primIdcs.size(); ++i)
        {
            _curveRanges[i] = make_pair(_primitives[_primIdcs[i]].startIdx, _primitives[_primIdcs[i]].endIdx);
            _curves[i] = _primitives[_primIdcs[i]].curve->clone();
        }

        //trim curves and ranges
        const int sampledPts = fitter.output<RESAMPLING>()->output->pts().size();
        for(int i = 0; i < (int)_continuities.size(); ++i)
        {
            if(_continuities[i] == 0)
                continue;

            //trim
            Vector2d trimPt = 0.5 * (_curves[i]->endPos() + _curves[i + 1]->startPos());
            if(!_primitives[_primIdcs[i]].isFixed())
                _curves[i]->trim(0, _curves[i]->project(trimPt));
            if(_curves.circular() || !_primitives[_primIdcs[i + 1]].isFixed())
                _curves[i + 1]->trim(_curves[i + 1]->project(trimPt), _curves[i + 1]->length());

            //the ranges over which error is computed should not overlap
            if(_continuities[i] == 1)
                swap(_curveRanges[i].second, _curveRanges[i + 1].first); //they overlap by 1 originally
            if(_continuities[i] == 2) //they overlap by 2
            {
                _curveRanges[i + 1].first = (_curveRanges[i + 1].first + 2) % sampledPts;
                _curveRanges[i].second = (_curveRanges[i].second + sampledPts - 2) % sampledPts;
            }
        }
    }

    vector<LSBoxConstraint> getConstraints() const
    {
        vector<LSBoxConstraint> out;
        VectorXd curParams = params();

        int curVar = 0;
        for(int i = 0; i < _curves.size(); ++i)
        {
            int primIdx = _primIdcs[i];
            if(_primitives[primIdx].isFixed())
            {
                int n = _curves[i]->numParams();
                for(int j = 0; j < n; ++j, ++curVar)
                    out.push_back(LSBoxConstraint(curVar, curParams[curVar], 0));

                continue;
            }

            //length must be at least half initial
            out.push_back(LSBoxConstraint(curVar + CurvePrimitive::LENGTH, _curves[i]->length() * 0.5, 1));

            //inflections
            //Debugging::get()->printf("idx = %d Type = %d, startSign = %d, endSign = %d", i, _curves[i]->getType(), _primitives[primIdx].startCurvSign, _primitives[primIdx].endCurvSign);
            if(_inflectionAccounting && _curves[i]->getType() >= CurvePrimitive::ARC)
            {
                int pi = _curves.toLinearIdx(i - 1);
                int ni = _curves.toLinearIdx(i + 1);
                bool c2toPrev = (pi >= 0 && _continuities[pi] == 2);
                bool c2toNext = (ni < _curves.size() && _continuities[i] == 2);
                bool line2Prev = c2toPrev && _curves[pi]->getType() == CurvePrimitive::LINE;
                bool line2Next = c2toNext && _curves[ni]->getType() == CurvePrimitive::LINE;

                int startSign = _primitives[primIdx].startCurvSign;
                int endSign = _primitives[primIdx].endCurvSign;

                if(!c2toPrev || (endSign == _primitives[_primIdcs[pi]].startCurvSign && startSign == endSign && !line2Prev))
                    out.push_back(LSBoxConstraint(curVar + CurvePrimitive::CURVATURE, 0., startSign));
                
                if(c2toPrev && c2toNext && startSign != endSign) //constrain the "easier one" of the endpoints
                {
                    double startViolation = -_curves[i]->startCurvature() * startSign;
                    double endViolation = -_curves[i]->endCurvature() * endSign;

                    //if we're C2 with a line, don't constrain that curvature
                    if(line2Prev)
                        startViolation = Parameters::infinity;
                    if(line2Next)
                        endViolation = Parameters::infinity;

                    if(startViolation < endViolation)
                        out.push_back(LSBoxConstraint(curVar + CurvePrimitive::CURVATURE, 0., startSign));
                    else
                        out.push_back(LSBoxConstraint(curVar + CurvePrimitive::DCURVATURE, 0., endSign));
                }

                if(!c2toNext && _curves[i]->getType() == CurvePrimitive::CLOTHOID)
                    out.push_back(LSBoxConstraint(curVar + CurvePrimitive::DCURVATURE, 0., endSign));
            }

            curVar += _curves[i]->numParams();
        }

        return out;
    }

#if SPARSE
    typedef MulticurveSparseEvalData EvalDataType;
#else
    typedef MulticurveDenseEvalData EvalDataType;
#endif

    int _iter;
    LSEvalData *createEvalData() { return new EvalDataType(); }
    void eval(const Eigen::VectorXd &x, LSEvalData *data)
    {
        setParams(x);
        EvalDataType *evalData = static_cast<EvalDataType *>(data);
        _evalError(evalData);
        _evalConstraints(evalData);
        //printf("Err: obj = %lf con = %lf\n", evalData->errVectorRef().norm(), evalData->conVectorRef().norm()); 
    }

    void setParams(const Eigen::VectorXd &x)
    {
        int curIdx = 0;
        for(int i = 0; i < (int)_curves.size(); ++i)
        {
            if(_curves[i]->getType() != CurvePrimitive::CLOTHOID)
            {
                _curves[i]->setParams(x.segment(curIdx, _curves[i]->numParams()));
            }
            else
            {
                VectorXd xm = x.segment(curIdx, _curves[i]->numParams());
                xm(CurvePrimitive::DCURVATURE) = (xm(CurvePrimitive::DCURVATURE) - xm(CurvePrimitive::CURVATURE)) / xm(CurvePrimitive::LENGTH);
                _curves[i]->setParams(xm);
            }
            curIdx += _curves[i]->numParams();
        }

        if(false)
        {
            char name[100];
            sprintf(name, "Out%d", _iter);
            for(int i = 0; i < _curves.size(); ++i)
                Debugging::get()->drawPrimitive(_curves[i], name, i, 2.);
        }

        ++_iter;
    }

    VectorXd params() const
    {
        int totParams = 0;
        for(int i = 0; i < (int)_curves.size(); ++i)
            totParams += _curves[i]->numParams();

        VectorXd out(totParams);

        int curIdx = 0;
        for(int i = 0; i < (int)_curves.size(); ++i)
        {
            if(_curves[i]->getType() != CurvePrimitive::CLOTHOID)
            {
                out.segment(curIdx, _curves[i]->numParams()) = _curves[i]->params();
            }
            else
            {
                VectorXd xm = _curves[i]->params(); 
                xm(CurvePrimitive::DCURVATURE) = xm(CurvePrimitive::CURVATURE) + xm(CurvePrimitive::DCURVATURE) * xm(CurvePrimitive::LENGTH);
                out.segment(curIdx, _curves[i]->numParams()) = xm;
            }
            curIdx += _curves[i]->numParams();
        }

        return out;
    }

    VectorC<CurvePrimitiveConstPtr> curves() const
    {
        VectorC<CurvePrimitiveConstPtr> out;
        copy(_curves.begin(), _curves.end(), back_inserter(out));
        return out;
    }

    double objective() const
    {
        double out = 0;

        for(int i = 0; i < (int)_curves.size(); ++i)
        {
            int csz = (int)_continuities.size();
            bool firstCorner = (!_closed && i == 0) || (_continuities[(i + csz - 1) % csz] == 0);
            bool lastCorner = (!_closed && i + 1 == (int)_curves.size()) || (_continuities[i] == 0);

            out += _errorComputer->computeError(_curves[i], _curveRanges[i].first, _curveRanges[i].second, firstCorner, lastCorner);
        }

        return out;
    }

private:
    void _evalError(EvalDataType *evalData)
    {
        vector<VectorXd> errVecs(_curves.size());
        vector<MatrixXd> errVecDers(_curves.size());

        size_t numErr = 0, numVar = 0;
        for(int i = 0; i < (int)_curves.size(); ++i)
        {
            int csz = (int)_continuities.size();
            bool firstCorner = (!_closed && i == 0) || (_continuities[(i + csz - 1) % csz] == 0);
            bool lastCorner = (!_closed && i + 1 == (int)_curves.size()) || (_continuities[i] == 0);

            _errorComputer->computeErrorVector(_curves[i], _curveRanges[i].first, _curveRanges[i].second,
                errVecs[i], &(errVecDers[i]), firstCorner, lastCorner);
            
            _curves[i]->toEndCurvatureDerivative(errVecDers[i]);

            numErr += errVecs[i].size();
            numVar += _curves[i]->numParams();
        }

        VectorXd &outErr = evalData->errVectorRef();

#if SPARSE
        EvalDataType::BlockVectorType &outErrDerBlocks = evalData->errDerBlocksRef();
        outErrDerBlocks.resize(_curves.size());
        outErr = VectorXd::Zero(numVar);
#else
        outErr = VectorXd::Zero(numErr);
        MatrixXd &outErrDer = evalData->errDerRef();
        outErrDer = MatrixXd::Zero(numErr, numVar);
#endif

        size_t curErr = 0, curVar = 0;
        for(int i = 0; i < (int)_curves.size(); ++i)
        {
            size_t nErr = errVecs[i].size();
            size_t nVar = errVecDers[i].cols();
#if SPARSE
            outErrDerBlocks[i] = errVecDers[i].transpose() * errVecDers[i];
            outErr.segment(curVar, nVar) = -errVecDers[i].transpose() * errVecs[i];
#else
            outErrDer.block(curErr, curVar, nErr, nVar) = errVecDers[i];
            outErr.segment(curErr, nErr) = errVecs[i];
#endif
            curErr += nErr;
            curVar += nVar;
        }
    }

    void _evalConstraints(EvalDataType *evalData)
    {
        VectorXd &outCon = evalData->conVectorRef();
        MatrixXd &outConDer = evalData->conDerRef();

        vector<VectorXd> conVecs(_continuities.size());
        vector<MatrixXd> conVecDers(_continuities.size());

        CurvePrimitive::EndDer endDer;

        size_t numCon = 0, numVar = 0;
        for(int i = 0; i < (int)_continuities.size(); ++i)
        {
            conVecs[i].resize(2 + _continuities[i]);
            conVecs[i].head<2>() = _curves[i]->endPos() - _curves[i + 1]->startPos();
            if(_continuities[i] >= 1)
                conVecs[i][2] = AngleUtils::toRange(_curves[i]->endAngle() - _curves[i + 1]->startAngle(), -PI);
            if(_continuities[i] == 2)
                conVecs[i][3] = _curves[i]->endCurvature() - _curves[i + 1]->startCurvature();

            _curves[i]->derivativeAtEnd(_continuities[i], endDer);
            conVecDers[i] = endDer;
            _curves[i]->toEndCurvatureDerivative(conVecDers[i]);

            numCon += conVecs[i].size();
            numVar += _curves[i]->numParams();
        }

        if(!_closed)
            numVar += _curves.back()->numParams();

        outCon = VectorXd::Zero(numCon);
        outConDer = MatrixXd::Zero(numCon, numVar);

        size_t curCon = 0, curVar = 0;
        for(int i = 0; i < (int)_continuities.size(); ++i)
        {
            size_t nCon = conVecs[i].size();
            size_t nVar = conVecDers[i].cols();
            outCon.segment(curCon, nCon) = conVecs[i];
            outConDer.block(curCon, curVar, nCon, nVar) = conVecDers[i];

            //now the derivatives for the second curve
            outConDer(curCon + 0, (curVar + nVar + CurvePrimitive::X) % numVar) = -1.;
            outConDer(curCon + 1, (curVar + nVar + CurvePrimitive::Y) % numVar) = -1.;
            if(nCon > 2)
                outConDer(curCon + 2, (curVar + nVar + CurvePrimitive::ANGLE) % numVar) = -1.;
            if(nCon > 3 && _curves[i + 1]->getType() != CurvePrimitive::LINE)
                outConDer(curCon + 3, (curVar + nVar + CurvePrimitive::CURVATURE) % numVar) = -1.;

            curCon += nCon;
            curVar += nVar;
        }
    }

    VectorC<CurvePrimitivePtr> _curves;
    const vector<FitPrimitive> &_primitives;
    vector<int> _primIdcs;
    vector<int> _continuities; //continuity[i] is between curves i and i + 1
    VectorC<pair<int, int> > _curveRanges;
    bool _closed;
    ErrorComputerConstPtr _errorComputer;
    bool _inflectionAccounting;
};

class DefaultCombiner : public Algorithm<COMBINING>
{
public:
    string name() const { return "Default"; }

protected:
    void _run(const Fitter &fitter, AlgorithmOutput<COMBINING> &out)
    {
        smart_ptr<const AlgorithmOutput<GRAPH_CONSTRUCTION> > graph = fitter.output<GRAPH_CONSTRUCTION>();
        const vector<FitPrimitive> &primitives = fitter.output<PRIMITIVE_FITTING>()->primitives;
        const vector<int> &path = fitter.output<PATH_FINDING>()->path;
        bool closed = fitter.output<CURVE_CLOSING>()->closed;

        if(path.empty())
            return; //no path

        VectorC<CurvePrimitiveConstPtr> outV;

        //if a single primitive
        if(graph->edges[path[0]].continuity == -1)
        {
            outV = VectorC<CurvePrimitiveConstPtr>(1, NOT_CIRCULAR);
            outV[0] = primitives[graph->edges[path[0]].startVtx].curve;
        }
        else //solve the nonlinear problem
        {
            MulticurveProblem problem(fitter);
            vector<LSBoxConstraint> constraints = problem.getConstraints();
            LSSolver solver(&problem, constraints);
            solver.setDefaultDamping(fitter.params().get(Parameters::COMBINE_DAMPING));
            solver.setMaxIter(50);
            solver.setIncreaseDampingAfter(5);
            solver.setDampingIncreaseFactor(1.5);

            VectorXd result = solver.solve(problem.params());
            problem.setParams(result);
            Debugging::get()->printf("Final objective = %lf", sqrt(problem.objective()));

            outV = problem.curves();
        }

        //==== track what happens to parameters ====
        out.parameters = fitter.output<RESAMPLING>()->parameters;
        PolylineConstPtr resampledCurve = fitter.output<RESAMPLING>()->output;
        const VectorC<Vector2d> &resampled = resampledCurve->pts();

        vector<int> finalPrimitives; //gather the indices of the graph vertices corresponding to the primitives
        for(int i = 0; i < (int)path.size(); ++i)
            finalPrimitives.push_back(graph->edges[path[i]].startVtx);
        if(outV.size() > (int)finalPrimitives.size())
            finalPrimitives.push_back(graph->edges[path.back()].endVtx);

        assert(outV.size() == finalPrimitives.size());

        vector<double> idxToParam(resampled.size()); //idx is the index into the resampled array
        vector<double> idxToDistSq(resampled.size(), 1e10);

        double lenSoFar = 0;
        for(int i = 0; i < outV.size(); ++i) //for each primitive see what projects to it
        {
            const FitPrimitive &primitive = primitives[graph->vertices[finalPrimitives[i]].primitiveIdx];
            for(int j = primitive.startIdx; ; ++j) //project each associated resampled point onto this primitive
            {
                if(j == (int)resampled.size()) //be careful with starts and ends of oversketched primitives
                {
                    if(closed) j = 0;
                    else
                        break;
                }

                if(j < 0)
                    continue;

                double proj = outV[i]->project(resampled[j]);
                double distSq = (resampled[j] - outV[i]->pos(proj)).squaredNorm();
                if(distSq < idxToDistSq[j])
                {
                    idxToDistSq[j] = distSq;
                    idxToParam[j] = proj + lenSoFar;
                }

                if(j == primitive.endIdx)
                    break;
            }
            lenSoFar += outV[i]->length();
        }
        double outputLength = lenSoFar;

        if(!closed)
        {
            idxToParam[0] = 0;
            idxToParam.back() = outputLength;
        }

        int minParamSample = min_element(idxToParam.begin(), idxToParam.end()) - idxToParam.begin();

        PiecewiseLinearMonotone prevToFinal(PiecewiseLinearMonotone::POSITIVE);
        //populate prevToFinal -- don't forget duplicating first point if closed
        prevToFinal.add(0, idxToParam[minParamSample]);
        lenSoFar = 0;
        for(VectorC<Vector2d>::Circulator ci = resampled.circulator(minParamSample + 1); !ci.done(); ++ci)
        {            
            lenSoFar += (*ci - *(ci - 1)).norm();
            double finalParam = 0;
            if(ci.index() == minParamSample)
                finalParam += outputLength;
            else
                idxToParam[ci.index()] = max(idxToParam[ci.index()], idxToParam[(ci - 1).index()]); //ensure monotonicity
            finalParam += idxToParam[ci.index()];
            prevToFinal.add(lenSoFar, finalParam);
        }
        //adjust parameters into range
        for(int i = 0; i < (int)out.parameters.size(); ++i)
        {
            out.parameters[i] -= resampledCurve->idxToParam(minParamSample);
            if(out.parameters[i] < 0)
                out.parameters[i] += resampledCurve->length();
        }
        prevToFinal.batchEval(out.parameters);

        //==== combine with what needs to be done w.r.t. oversketching ====
        smart_ptr<const AlgorithmOutput<OVERSKETCHING> > osOutput = fitter.output<OVERSKETCHING>();
        VectorC<CurvePrimitiveConstPtr> outFinal(0, osOutput->finallyClose ? CIRCULAR : NOT_CIRCULAR);

        if(osOutput->toPrepend)
        {
            outFinal.insert(outFinal.end(), osOutput->toPrepend->primitives().begin(), osOutput->toPrepend->primitives().end() - 1);
            for(int i = 0; i < (int)out.parameters.size(); ++i)
                out.parameters[i] += (osOutput->toPrepend->length() - osOutput->toPrepend->primitives().back()->length());
        }
        outFinal.insert(outFinal.end(), outV.begin(), outV.end());

        if(osOutput->toAppend)
        {
            if(!osOutput->finallyClose)
            {
                outFinal.insert(outFinal.end(), osOutput->toAppend->primitives().begin() + 1, osOutput->toAppend->primitives().end());
            }
            else
            {
                if(osOutput->toAppend->primitives().size() >= 2)
                {
                    outFinal.insert(outFinal.end(), osOutput->toAppend->primitives().begin() + 1, osOutput->toAppend->primitives().end() - 1);
                }
                else //start and end curve is the same curve -- its original length is the one toAppend curve
                {
                    //get rid of last curve and possibly extend the first one
                    double lastCurveLen = outFinal.back()->length();
                    double firstCurveLen = outFinal[0]->length();
                    double origLen = osOutput->toAppend->primitives()[0]->length();
                    outFinal.pop_back();
                    //now extend the first curve to the combined length
                    outFinal[0] = outFinal[0]->trimmed(origLen - lastCurveLen, firstCurveLen);
                    for(int i = 0; i < (int)out.parameters.size(); ++i)
                        out.parameters[i] -= (origLen - lastCurveLen);
                }
            }
        }

        out.output = new PrimitiveSequence(outFinal);

#if 1
        for(int i = 0; i < (int)out.parameters.size(); ++i)
        {
            double paramOrig = fitter.originalSketch()->idxToParam(i);
            Debugging::get()->drawLine(fitter.originalSketch()->pts()[i], out.output->pos(out.parameters[i]), Vector3d(1, 0, 1), "Correspondence");
        }
#endif

    }
};

void Algorithm<COMBINING>::_initialize()
{
    new DefaultCombiner();
}

END_NAMESPACE_Cornu


