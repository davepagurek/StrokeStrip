/*--
    Parameters.h  

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

#ifndef CORNUCOPIA_PARAMETERS_H_INCLUDED
#define CORNUCOPIA_PARAMETERS_H_INCLUDED

#include <vector>
#include <string>

namespace Cornu
{

/*
 * The Parameters object controls the execution of the fitting.  It specifies what the costs are for the
 * shortest path graph, and it specifies many internal parameters.  It also specifies the choice of
 * different algorithms for fitting stages.  Using the Parameters object, you can make the fitter
 * use only arcs, increase the tolerance, etc.  The default constructor provides reasonable parameters
 * for the G2 clothoid fit.  Several presets (defined in void Parameters::_initializePresets) allow
 * other useful behaviors.
 */
class Parameters
{
public:
    enum ParameterType
    {
        //The costs are used for computing edge and vertex weights in the graph.
        //Each comment is the effect increasing the parameter has on the results
        LINE_COST, //Uses fewer lines
        ARC_COST, //Uses fewer arcs
        CLOTHOID_COST, //Uses fewer clothoids
        G0_COST, //Uses fewer G0 transitions
        G1_COST, //Uses fewer G1 transitions
        G2_COST, //Uses fewer G2 transitions
        ERROR_COST, //Favors better fit, but uses more primitives
        SHORTNESS_COST, //Avoids shorter segments
        INFLECTION_COST, //Avoids inflections
        //This is not a real parameter--below are parameters the user should not control
        INTERNAL_PARAMETERS_MARKER, 
        PIXEL_SIZE, //Increasing the pixel size is almost equivalent to scaling down the curve by the same factor.  Many of the other parameters are in these "pixels"
        SMALL_CURVE_PIXELS, //Curves below this size get scaled up
        LARGE_CURVE_PIXELS, //Curves above this size get scaled down
        MAX_RESCALE, //Maximum scale applied to the curve based on size -- set this to 1 to prevent rescaling
        MIN_PRELIM_LENGTH, //The rate at which the curve gets resampled the first time around.  Increasing this decreases accuracy but makes the algorithm less sensitive to noise.
        DP_CUTOFF, //The stopping criterion for the Douglas-Peucker part of preliminary resampling.  Increasing this decreases accuracy but makes the algorithm less sensitive to noise.
        CLOSEDNESS_THRESHOLD, //Roughly how far the endpoints of a curve need to be for the curve to be considered closed.
        MINIMUM_CORNER_SPACING, //Two corners will not be identified at smaller than this distance
        CORNER_NEIGHBORHOOD, //How much arclength on each side of a point on the curve is used to determine whether it is a corner
        DENSE_SAMPLING_STEP,
        CORNER_SCALES,
        CORNER_THRESHOLD,
        MAX_SAMPLING_INTERVAL, //When resampling, maximum distance between adjacent samples, even on a straight line
        CURVATURE_ESTIMATE_REGION, //How much arclength is used to estimate the local curvature at a point.  Decreasing this tends to make the estimated curvature larger.
        POINTS_PER_CIRCLE, //If the curvature were constant, this would be the sampling rate.  Increasing this increases result quality at the cost of performance
        MAX_SAMPLE_RATE_SLOPE, //How much the curvature-dependent sampling rate can change per unit arclength
        ERROR_THRESHOLD, //Primitves whose error is greater than are discarded.  Decreasing this increases performance, but may hurt quality and eventually result in a fit failure
        SHORTNESS_THRESHOLD, //Primitives below this length are considered "short" for the purposes of the shortness cost
        TWO_CURVE_CURVATURE_ADJUST, //When combining two curves and matching their curvature, how much to compensate with the curvature at the opposite endpoints
        CURVE_ADJUST_DAMPING, //How much regularization is added to the solver for edge validation--increasing this makes the solver more stable, but converge slower
        REDUCE_GRAPH_EVERY, //How many invalid paths are found before the A* heuristic is recomputed.  Setting this too high or too low hurts performance.
        COMBINE_DAMPING, //How much regularization is added to the solver for the final combine--increasing this makes the solver more stable, but converge slower
        OVERSKETCH_THRESHOLD //How far the endpoints need to be from the base curve for them to be considered on the curve
    };

    enum Preset
    {
        DEFAULT,
        LOOSE,
        ACCURATE,
        POLYLINE,
        LINES_AND_ARCS,
        CLOTHOID_ONLY,
        NUM_PRESETS //must be last
    };

    Parameters(const std::string &name = std::string());
    Parameters(const Parameters &parameters) : _name(parameters._name), _values(parameters._values), _algorithms(parameters._algorithms) {}
    Parameters(Preset preset) { *this = presets()[preset]; }

    void set(ParameterType param, double val) { _values[param] = val; }
    double get(ParameterType param) const { return _values[param]; }
    void setAlgorithm(int stage, int algorithm) { _algorithms[stage] = algorithm; }
    int getAlgorithm(int stage) const { return _algorithms[stage]; }

    const std::string &name() const { return _name; } 
    void setName(const std::string &name) { _name = name; }

    bool operator==(const Parameters &other) const { return _values == other._values && _algorithms == other._algorithms; }
    bool operator!=(const Parameters &other) const { return !((*this) == other); }

private:
    std::string _name;
    std::vector<double> _values;
    std::vector<int> _algorithms;

    //===== static description data =====
public:
    struct Parameter
    {
        Parameter(ParameterType inType, const std::string &inName, double inMin, double inMax, double inDefault, bool inInfinityAllowed = true)
            : type(inType), typeName(inName), min(inMin), max(inMax), defaultVal(inDefault), infinityAllowed(inInfinityAllowed) {}
        Parameter(ParameterType inType, const std::string &inName, double inValue) //for internal parameters
            : type(inType), typeName(inName), min(inValue), max(inValue), defaultVal(inValue), infinityAllowed(false) {}

        ParameterType type;
        std::string typeName;
        double min, max, defaultVal;
        bool infinityAllowed;
    };

    static const double infinity;

    static const std::vector<Parameter> &parameters() { _initializeParameters(); return _parameters; }
    static const std::vector<Parameters> &presets() { _initializePresets(); return _presets; }

private:
    static void _initializeParameters(); 
    static void _initializePresets(); 
    static std::vector<Parameter> _parameters;
    static std::vector<Parameters> _presets;
};

} //end of namespace Cornu

#endif //CORNUCOPIA_PARAMETERS_H_INCLUDED
