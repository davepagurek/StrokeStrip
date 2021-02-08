/*--
    Algorithm.h  

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

#ifndef CORNUCOPIA_ALGORITHM_H_INCLUDED
#define CORNUCOPIA_ALGORITHM_H_INCLUDED

#include "defs.h"
#include "smart_ptr.h"
#include <vector>

NAMESPACE_Cornu

//This header file defines the stages and base classes for algorithms in the fitting pipeline.
//Fitting proceeds in stages, each one having access to the outputs of all the previous stages
//through the Fitter class.

enum AlgorithmStage
{
    SCALE_DETECTION,
    PRELIM_RESAMPLING,
    CURVE_CLOSING,
    OVERSKETCHING,
    CORNER_DETECTION,
    RESAMPLING,
    ERROR_COMPUTER,
    PRIMITIVE_FITTING,
    GRAPH_CONSTRUCTION,
    PATH_FINDING,
    COMBINING,
    NUM_ALGORITHM_STAGES //must be last
};

struct AlgorithmOutputBase : public smart_base
{
};

CORNU_SMART_TYPEDEFS(AlgorithmOutputBase);

//The output of an algorithm stage.  It is specialized for every stage.
template<int AlgStage>
struct AlgorithmOutput : public AlgorithmOutputBase
{
};

class Fitter;

template<int AlgStage>
class Algorithm
{
};

class AlgorithmBase
{
public:
    virtual std::string name() const { return "Default"; }
    virtual std::string stageName() const = 0;
    virtual AlgorithmOutputBasePtr run(const Fitter &) = 0;

    static int numAlgorithmsForStage(AlgorithmStage stage) { return (int)_getAlgorithms()[stage].size(); }
    static AlgorithmBase *get(AlgorithmStage stage, int algorithm) { return _getAlgorithms()[stage][algorithm]; }

protected:
    static const std::vector<std::vector<AlgorithmBase *> > &_getAlgorithms();
    static void _addAlgorithm(int stage, AlgorithmBase *algorithm);

private:
    static bool _initializationFinished;
    static void _initialize();
    static std::vector<std::vector<AlgorithmBase *> > _algorithms;
};

template<int AlgStage>
class AlgorithmBaseTemplate : public AlgorithmBase
{
public:
    //override
    AlgorithmOutputBasePtr run(const Fitter &fitter)
    {
        smart_ptr<AlgorithmOutput<AlgStage> > out = new AlgorithmOutput<AlgStage>();
        _run(fitter, *out);
        return out;
    }

    static std::vector<std::string> names()
    {
        std::vector<std::string> out;
        for(int i = 0; i < (int)_getAlgorithms().size(); ++i)
            out.push_back(_getAlgorithms()[AlgStage][i]->name());

        return out;
    }

protected:
    AlgorithmBaseTemplate()
    {
        _addAlgorithm(AlgStage, this);
    }

    virtual void _run(const Fitter &fitter, AlgorithmOutput<AlgStage> &out) = 0;
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_ALGORITHM_H_INCLUDED
