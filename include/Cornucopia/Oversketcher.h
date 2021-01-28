/*--
    Oversketcher.h  

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

#ifndef CORNUCOPIA_OVERSKETCHER_H_INCLUDED
#define CORNUCOPIA_OVERSKETCHER_H_INCLUDED

#include "defs.h"
#include "Algorithm.h"

NAMESPACE_Cornu

CORNU_SMART_FORW_DECL(Polyline);
CORNU_SMART_FORW_DECL(PrimitiveSequence);

template<>
struct AlgorithmOutput<OVERSKETCHING> : public AlgorithmOutputBase
{
    PolylineConstPtr output;
    std::vector<double> parameters; //parameters[i] is the parameter in output of the original point with index i
    CurvePrimitiveConstPtr startCurve;
    CurvePrimitiveConstPtr endCurve;
    int startContinuity;
    int endContinuity;
    PrimitiveSequenceConstPtr toAppend;
    PrimitiveSequenceConstPtr toPrepend;
    bool finallyClose;
};

template<>
class Algorithm<OVERSKETCHING> : public AlgorithmBaseTemplate<OVERSKETCHING>
{
public:
    //override
    std::string stageName() const { return "Oversketcher"; }

private:
    friend class AlgorithmBase;
    static void _initialize();
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_OVERSKETCHER_H_INCLUDED
