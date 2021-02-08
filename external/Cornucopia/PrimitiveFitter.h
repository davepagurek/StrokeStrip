/*--
    PrimitiveFitter.h  

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

#ifndef CORNUCOPIA_PRIMITIVEFITTER_H_INCLUDED
#define CORNUCOPIA_PRIMITIVEFITTER_H_INCLUDED

#include "defs.h"
#include "Algorithm.h"

NAMESPACE_Cornu

CORNU_SMART_FORW_DECL(CurvePrimitive);

struct FitPrimitive
{
    FitPrimitive() : fixed(false) {}

    CurvePrimitivePtr curve;
    int startIdx;
    int endIdx;
    int numPts;
    double error;
    int startCurvSign;
    int endCurvSign;
    bool fixed;

    bool isStartCurve() const { return startIdx == -1; }
    bool isEndCurve() const { return fixed && startIdx != -1; }
    bool isFixed() const { return fixed; }
};

template<>
struct AlgorithmOutput<PRIMITIVE_FITTING> : public AlgorithmOutputBase
{
    std::vector<FitPrimitive> primitives;
};

template<>
class Algorithm<PRIMITIVE_FITTING> : public AlgorithmBaseTemplate<PRIMITIVE_FITTING>
{
public:
    //override
    std::string stageName() const { return "Primitive Fitter"; }

private:
    friend class AlgorithmBase;
    static void _initialize();
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_PRIMITIVEFITTER_H_INCLUDED
