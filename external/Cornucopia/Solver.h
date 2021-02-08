/*--
    Solver.h  

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

#ifndef CORNUCOPIA_SOLVER_H_INCLUDED
#define CORNUCOPIA_SOLVER_H_INCLUDED

#include "defs.h"
#include <vector>
#include <set>
#include <Eigen/Core>

NAMESPACE_Cornu

struct LSBoxConstraint
{
    LSBoxConstraint(int inIndex, double inValue, int inSign) : index(inIndex), value(inValue), sign(inSign) {}
    bool operator<(const LSBoxConstraint &other) const { return index < other.index; }

    int index; //of variable
    double value;
    int sign; //variable greater than value if sign is positive, less than value if negative, equal if zero
};

class LSEvalData
{
public:
    virtual ~LSEvalData() {}

    virtual double error() const = 0;
    virtual void solveForDelta(double damping, Eigen::VectorXd &out, std::set<LSBoxConstraint> &constraints) = 0;

    //debugging functions for derivative check
    virtual Eigen::VectorXd errVec() const { return Eigen::VectorXd(); }
    virtual Eigen::MatrixXd errVecDer() const { return Eigen::MatrixXd(); }
};

class LSProblem
{
public:
    virtual ~LSProblem() {}

    virtual double error(const Eigen::VectorXd &x, LSEvalData *data) { eval(x, data); return data->error(); }
    virtual LSEvalData *createEvalData() = 0;
    virtual void eval(const Eigen::VectorXd &x, LSEvalData *data) = 0;
};

class LSSolver
{
public:
    LSSolver(LSProblem *problem, const std::vector<LSBoxConstraint> &constraints);

	Eigen::VectorXd solve(const Eigen::VectorXd &guess, bool *is_valid = nullptr);
    void setDefaultDamping(double damping) { _damping = damping; }
    void setMaxIter(int maxIter) { _maxIter = maxIter; }
    void setIncreaseDampingAfter(int iter) { _increaseDampingAfter = iter; }
    void setDampingIncreaseFactor(double factor) { _dampingIncreaseFactor = factor; }

    bool verifyDerivatives(const Eigen::VectorXd &pt, double eps = 1e-6) const;

private:
    int _project(const Eigen::VectorXd &from, Eigen::VectorXd &x, const std::set<LSBoxConstraint> &activeSet); //returns the index of the constraint
    std::set<LSBoxConstraint> _clamp(Eigen::VectorXd &x);

    LSProblem *_problem;
    std::vector<LSBoxConstraint> _constraints;
    double _damping;
    int _maxIter;
    int _increaseDampingAfter;
    double _dampingIncreaseFactor;
};

class LSDenseEvalData : public LSEvalData
{
public:
    //overrides
    double error() const { return _err.squaredNorm(); }
    void solveForDelta(double damping, Eigen::VectorXd &out, std::set<LSBoxConstraint> &constraints);
    Eigen::VectorXd errVec() const { return _err; }
    Eigen::MatrixXd errVecDer() const { return _errDer; }

    Eigen::VectorXd &errVectorRef() { return _err; }
    Eigen::MatrixXd &errDerRef() { return _errDer; }
private:
    Eigen::VectorXd _err;
    Eigen::MatrixXd _errDer;
};


END_NAMESPACE_Cornu

#endif //CORNUCOPIA_SOLVER_H_INCLUDED
