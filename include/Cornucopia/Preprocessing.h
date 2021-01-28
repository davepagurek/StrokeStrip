/*--
    Preprocessing.h  

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

#ifndef CORNUCOPIA_PREPROCESSING_H_INCLUDED
#define CORNUCOPIA_PREPROCESSING_H_INCLUDED

#include "defs.h"
#include "Algorithm.h"

NAMESPACE_Cornu

CORNU_SMART_FORW_DECL(Polyline)

template<>
struct AlgorithmOutput<SCALE_DETECTION> : public AlgorithmOutputBase
{
    AlgorithmOutput() : scale(1.) {}

    double scale;
};

template<>
struct AlgorithmOutput<PRELIM_RESAMPLING> : public AlgorithmOutputBase
{
    PolylineConstPtr output;
    std::vector<double> parameters; //parameters[i] is the parameter in output of the original point with index i
};

template<>
struct AlgorithmOutput<CURVE_CLOSING> : public AlgorithmOutputBase
{
    PolylineConstPtr output;
    bool closed;
    std::vector<double> parameters; //parameters[i] is the parameter in output of the original point with index i
};

template<>
class Algorithm<SCALE_DETECTION> : public AlgorithmBaseTemplate<SCALE_DETECTION>
{
public:
    //override
    std::string stageName() const { return "Scale Detector"; }

private:
    friend class AlgorithmBase;
    static void _initialize();
};

template<>
class Algorithm<PRELIM_RESAMPLING> : public AlgorithmBaseTemplate<PRELIM_RESAMPLING>
{
public:
    //override
    std::string stageName() const { return "Preliminary Resampling"; }

private:
    friend class AlgorithmBase;
    static void _initialize();
};

template<>
class Algorithm<CURVE_CLOSING> : public AlgorithmBaseTemplate<CURVE_CLOSING>
{
public:
    //override
    std::string stageName() const { return "Closedness Detector"; }

private:
    friend class AlgorithmBase;
    static void _initialize();
};


END_NAMESPACE_Cornu

#endif //CORNUCOPIA_PREPROCESSING_H_INCLUDED
