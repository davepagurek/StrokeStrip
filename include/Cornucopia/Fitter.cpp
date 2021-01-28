/*--
    Fitter.cpp  

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

#include "Fitter.h"
#include "Preprocessing.h"
#include "Polyline.h"
#include "Resampler.h"
#include "Combiner.h"
#include "PrimitiveSequence.h"

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

void Fitter::run()
{
    Debugging::get()->clear();
    Debugging::get()->printf("============= Starting =============");
    Debugging::get()->drawCurve(_originalSketch, Vector3d(0, 0, 0), "Original Sketch", 2., Debugging::DOTTED);
    Debugging::get()->startTiming("Total");

    for(int i = 0; i < NUM_ALGORITHM_STAGES; ++i)
    {
        if(!(_outputs[i]))
        {
            std::string stageName = AlgorithmBase::get((AlgorithmStage)i, 0)->stageName();
            Debugging::get()->startTiming(stageName);
            _runStage((AlgorithmStage)i);
            if(Debugging::get()->getTimeElapsed(stageName) > 0.001) //only print significant times
                Debugging::get()->elapsedTime(stageName);
        }
    }
    Debugging::get()->elapsedTime("Total");

    if(Debugging::get()->isDebuggingOn() && finalOutput())
    {
        // Now output some the final curve and a normal field for debugging
        PrimitiveSequenceConstPtr out = finalOutput();
        for(int i = 0; i < out->primitives().size(); ++i)
        {
            Debugging::get()->drawPrimitive(out->primitives()[i], "Final Result Color", i, 3.);
            Debugging::get()->drawCurve(out->primitives()[i], Vector3d(0, 0, 0), "Final Result");
            Debugging::get()->drawCurvatureField(out->primitives()[i], Vector3d(1, 0, 0), "Normal Field");
        }
    }
}

void Fitter::_runStage(AlgorithmStage stage)
{
    _outputs[stage] = AlgorithmBase::get(stage, _params.getAlgorithm(stage))->run(*this);
}

void Fitter::_clearBefore(AlgorithmStage stage)
{
    for(int i = stage; i < NUM_ALGORITHM_STAGES; ++i)
        _outputs[i] = AlgorithmOutputBasePtr();
}

double Fitter::scale() const
{
    return output<SCALE_DETECTION>()->scale * _params.get(Parameters::PIXEL_SIZE);
}

double Fitter::scaledParameter(Parameters::ParameterType param) const
{
    return _params.get(param) * scale();
}

PrimitiveSequenceConstPtr Fitter::finalOutput() const
{
    return output<COMBINING>()->output;
}

const vector<double> &Fitter::originalSketchToFinalParameters() const
{
    return output<COMBINING>()->parameters;
}


END_NAMESPACE_Cornu


