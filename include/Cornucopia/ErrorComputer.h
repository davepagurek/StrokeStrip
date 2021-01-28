/*--
    ErrorComputer.h  

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

#ifndef CORNUCOPIA_ERRORCOMPUTER_H_INCLUDED
#define CORNUCOPIA_ERRORCOMPUTER_H_INCLUDED

#include "defs.h"
#include "Algorithm.h"

NAMESPACE_Cornu

CORNU_SMART_FORW_DECL(CurvePrimitive);

class ErrorComputer : public smart_base
{
public:
    virtual ~ErrorComputer() {}
    //Computes the summed squared distance from the samples between from and to (incl.) to the given curve, weighted by distance between samples.
    //If firstToEndpoint is true, computes the distance of the first sample to the start of the curve (end if reversed==true) instead of to the
    //projection.  Analogously for lastToEndpoint.
    virtual double computeError(CurvePrimitiveConstPtr curve, int from, int to,
                                bool firstToEndpoint = true, bool lastToEndpoint = true, bool reversed = false) const = 0;
    //Computes the individual error terms in a format suitable for optimization.
    virtual void computeErrorVector(CurvePrimitiveConstPtr curve, int from, int to,
                                    Eigen::VectorXd &outError, Eigen::MatrixXd *outErrorDer = NULL,
                                    bool firstToEndpoint = true, bool lastToEndpoint = true, bool reversed = false) const = 0;
    //Computes the error to be used in the graph weight--by default, the squared maximum distance to the curve
    virtual double computeErrorForCost(CurvePrimitiveConstPtr curve, int from, int to,
                                       bool firstToEndpoint = true, bool lastToEndpoint = true, bool reversed = false) const = 0;
};

CORNU_SMART_TYPEDEFS(ErrorComputer);

template<>
struct AlgorithmOutput<ERROR_COMPUTER> : public AlgorithmOutputBase
{
    ErrorComputerConstPtr errorComputer;
};

template<>
class Algorithm<ERROR_COMPUTER> : public AlgorithmBaseTemplate<ERROR_COMPUTER>
{
public:
    //override
    std::string stageName() const { return "Error Computer"; }

private:
    friend class AlgorithmBase;
    static void _initialize();
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_ERRORCOMPUTER_H_INCLUDED
