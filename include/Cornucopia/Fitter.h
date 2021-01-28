/*--
    Fitter.h  

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

#ifndef CORNUCOPIA_FITTER_H_INCLUDED
#define CORNUCOPIA_FITTER_H_INCLUDED

#include "defs.h"
#include "Parameters.h"
#include "Algorithm.h"

NAMESPACE_Cornu

CORNU_SMART_FORW_DECL(Polyline);
CORNU_SMART_FORW_DECL(PrimitiveSequence);

class Fitter
{
public:
    Fitter() : _outputs(NUM_ALGORITHM_STAGES) {}

    const Parameters &params() const { return _params; }
    void setParams(const Parameters &params) { _params = params; _clearBefore(SCALE_DETECTION); }

    PolylineConstPtr originalSketch() const { return _originalSketch; }
    void setOriginalSketch(PolylineConstPtr originalSketch) { _originalSketch = originalSketch; _clearBefore(SCALE_DETECTION); }

    PrimitiveSequenceConstPtr oversketchBase() const { return _oversketchBase; }
    void setOversketchBase(PrimitiveSequenceConstPtr oversketchBase) { _oversketchBase = oversketchBase; _clearBefore(SCALE_DETECTION); }

    template<int AlgStage>
    smart_ptr<const AlgorithmOutput<AlgStage> > output() const
    {
        return static_pointer_cast<const AlgorithmOutput<AlgStage> >(_outputs[AlgStage]);
    }

    void run();

    PrimitiveSequenceConstPtr finalOutput() const; //returns null if fitting failed for some reason
    const std::vector<double> &originalSketchToFinalParameters() const; //returns a vector that for each original sketch point has the final parameter value

    double scale() const;  //returns the scale (pixel size * detected scale)
    double scaledParameter(Parameters::ParameterType param) const;

private:
    void _runStage(AlgorithmStage stage);
    void _clearBefore(AlgorithmStage stage);

    PrimitiveSequenceConstPtr _oversketchBase;
    PolylineConstPtr _originalSketch;
    Parameters _params;

    std::vector<AlgorithmOutputBasePtr> _outputs;
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_FITTER_H_INCLUDED
