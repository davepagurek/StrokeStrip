/*--
    PiecewiseLinearUtils.h  

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

#ifndef CORNUCOPIA_PIECEWISELINEARUTILS_H_INCLUDED
#define CORNUCOPIA_PIECEWISELINEARUTILS_H_INCLUDED

#include "defs.h"
#include <vector>
#include <set>

NAMESPACE_Cornu

//Approximates a monotone function piecewise-linearly.  Supports adding points and solving for x given y.
//Used for solving for the adjustment to the sampling rate that leads to an integer number of samples and
//other resampling stuff.
class PiecewiseLinearMonotone
{
public:
    enum Sign { POSITIVE, NEGATIVE };

    PiecewiseLinearMonotone(Sign sign) : _sign(sign ? -1. : 1.) {}

    void add(double x, double y);
    bool eval(double x, double &outY) const; //returns true if evaluation is successful
    bool invert(double y, double &outX) const; //returns true if inversion is successful

    double minX() const;
    double maxX() const;

    bool batchEval(std::vector<double> &inXoutY) const; //returns false if any evaluation fails
private:
    struct PLPoint //allows comparison by x or by y -- order in the map should be the same
    {
        PLPoint(double inX, double inY) : x(inX), y(inY), compareByY(false) {}
        PLPoint(double inY) : x(-1), y(inY), compareByY(true) {}

        bool operator<(const PLPoint &other) const { if(compareByY || other.compareByY) return y < other.y; else return x < other.x; }

        double x, y;
        bool compareByY;
    };

    double _sign;
    std::set<PLPoint> _points;
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_PIECEWISELINEARUTILS_H_INCLUDED
