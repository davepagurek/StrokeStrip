/*--
    Parameters.cpp  

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

#include "Parameters.h"
#include "Algorithm.h"

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

Parameters::Parameters(const string &name)
: _name(name), _algorithms(NUM_ALGORITHM_STAGES, 0)
{
    _initializeParameters();
    _values.resize(_parameters.size());
    for(int i = 0; i < (int)_parameters.size(); ++i)
        _values[i] = _parameters[i].defaultVal;
}

void Parameters::_initializeParameters()
{
    static bool initialized = false;
    if(initialized)
        return;
    initialized = true;

    _parameters.push_back(Parameter(LINE_COST, "Line cost", 0., 20., 7.5));
    _parameters.push_back(Parameter(ARC_COST, "Arc cost", 0., 30., 9.));
    _parameters.push_back(Parameter(CLOTHOID_COST, "Clothoid cost", 0., 50., 15.));
    _parameters.push_back(Parameter(G0_COST, "G0 cost", 0., 50., infinity));
    _parameters.push_back(Parameter(G1_COST, "G1 cost", 0., 50., infinity));
    _parameters.push_back(Parameter(G2_COST, "G2 cost", 0., 50., 0.));
    _parameters.push_back(Parameter(ERROR_COST, "Error cost", 0., 10., 1.));
    _parameters.push_back(Parameter(SHORTNESS_COST, "Shortness cost", 0., 10., 2.));
    _parameters.push_back(Parameter(INFLECTION_COST, "Inflection cost", 0., 100., 20.));

    _parameters.push_back(Parameter(INTERNAL_PARAMETERS_MARKER, "NOT A PARAMETER", 0., 0., 0.));
    _parameters.push_back(Parameter(PIXEL_SIZE, "Pixel size", 1.));
    _parameters.push_back(Parameter(SMALL_CURVE_PIXELS, "Small curve pixels", 200.));
    _parameters.push_back(Parameter(LARGE_CURVE_PIXELS, "Large curve pixels", 500.));
    _parameters.push_back(Parameter(MAX_RESCALE, "Max rescale", 2.));
    _parameters.push_back(Parameter(MIN_PRELIM_LENGTH, "Min prelim length", 2.));
    _parameters.push_back(Parameter(DP_CUTOFF, "Douglas-Peucker cutoff", 3.));
    _parameters.push_back(Parameter(CLOSEDNESS_THRESHOLD, "Closedness threshold", 20.));
    _parameters.push_back(Parameter(MINIMUM_CORNER_SPACING, "Min corner spacing", 5.));
    _parameters.push_back(Parameter(CORNER_NEIGHBORHOOD, "Corner neighborhood", 15.));
    _parameters.push_back(Parameter(DENSE_SAMPLING_STEP, "Dense sampling step", 1.));
    _parameters.push_back(Parameter(CORNER_SCALES, "Num corner scales (int)", 5.));
    _parameters.push_back(Parameter(CORNER_THRESHOLD, "Corner angle threshold", PI * 0.25));
    _parameters.push_back(Parameter(MAX_SAMPLING_INTERVAL, "Maximum sampling interval", 50.));
    _parameters.push_back(Parameter(CURVATURE_ESTIMATE_REGION, "Curvature estimate region", 20.));
    _parameters.push_back(Parameter(POINTS_PER_CIRCLE, "Points per circle", 20.));
    _parameters.push_back(Parameter(MAX_SAMPLE_RATE_SLOPE, "Max sample rate slope", 0.4));
    _parameters.push_back(Parameter(ERROR_THRESHOLD, "Error Threshold", 4.));
    _parameters.push_back(Parameter(SHORTNESS_THRESHOLD, "Shortness Threshold", 50.));
    _parameters.push_back(Parameter(TWO_CURVE_CURVATURE_ADJUST, "Two-Curve Adjustment Point", 2.));
    _parameters.push_back(Parameter(CURVE_ADJUST_DAMPING, "Curve Adjust Damping", 1.));
    _parameters.push_back(Parameter(REDUCE_GRAPH_EVERY, "Reduce Graph Every", 10.));
    _parameters.push_back(Parameter(COMBINE_DAMPING, "Combine Damping", 2.));
    _parameters.push_back(Parameter(OVERSKETCH_THRESHOLD, "Oversketch Threshold", 15.));
}

void Parameters::_initializePresets()
{
    static bool initialized = false;
    if(initialized)
        return;
    initialized = true;

    _initializeParameters();

    _presets.resize(NUM_PRESETS);

    //default
    Parameters defaultParams("Default (G2)");
    _presets[DEFAULT] = defaultParams;

    //Loose
    Parameters loose = defaultParams;
    loose._name = "Loose (G2)";
    loose.set(ERROR_COST, 0.2);
    loose.set(ERROR_THRESHOLD, 8.);
    _presets[LOOSE] = loose;

    //Accurate
    Parameters accurate = defaultParams;
    accurate._name = "Accurate (G2)";
    accurate.set(ERROR_COST, 5.);
    accurate.set(SHORTNESS_COST, 1.);
    accurate.set(CURVATURE_ESTIMATE_REGION, 10.);
    accurate.set(POINTS_PER_CIRCLE, 30.);
    accurate.set(ERROR_THRESHOLD, 2.);
    _presets[ACCURATE] = accurate;

    //polyline
    Parameters polyline = defaultParams;
    polyline._name = "Polyline";
    polyline.set(LINE_COST, 5.);
    polyline.set(ARC_COST, infinity);
    polyline.set(CLOTHOID_COST, infinity);
    polyline.set(G0_COST, 0.);
    polyline.set(G1_COST, infinity);
    polyline.set(G2_COST, infinity);
    polyline.set(SHORTNESS_COST, 0);
    _presets[POLYLINE] = polyline;

    //lines and arcs
    Parameters linesarcs = defaultParams;
    linesarcs._name = "Lines And Arcs";
    linesarcs.set(LINE_COST, 8.);
    linesarcs.set(ARC_COST, 12.);
    linesarcs.set(CLOTHOID_COST, infinity);
    linesarcs.set(G1_COST, 0.);
    linesarcs.set(G2_COST, infinity);
    linesarcs.set(SHORTNESS_COST, 0);
    _presets[LINES_AND_ARCS] = linesarcs;

    //Clothoid only
    Parameters co = defaultParams;
    co._name = "Clothoid Only";
    co.set(LINE_COST, infinity);
    co.set(ARC_COST, infinity);
    _presets[CLOTHOID_ONLY] = co;
}

const double Parameters::infinity = 1e30;
vector<Parameters::Parameter> Parameters::_parameters;
vector<Parameters> Parameters::_presets;

END_NAMESPACE_Cornu


