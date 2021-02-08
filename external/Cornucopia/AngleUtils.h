/*--
    AngleUtils.h  

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

#ifndef CORNUCOPIA_ANGLEUTILS_H_INCLUDED
#define CORNUCOPIA_ANGLEUTILS_H_INCLUDED

#include "defs.h"

#include <Eigen/Core>

NAMESPACE_Cornu

/*
    This class has a bunch of static methods for working with angles
*/
class AngleUtils
{
public:
    //signed angle from (1, 0) to v
    static double angle(const Eigen::Vector2d &v) { return atan2(v[1], v[0]); }
    //signed angle from v1 to v2
    static double angle(const Eigen::Vector2d &v1, const Eigen::Vector2d &v2) { return atan2(v1[0] * v2[1] - v1[1] * v2[0], v1.dot(v2)); }

    //These functions bring an angle into the range [0,2Pi] or [rangeStart, rangeStart+2Pi]
    //They assume we're not too far on the negative side of the range
    static double toRange(double angle) { return fmod(angle + 8 * PI, TWOPI); }
    static double toRange(double angle, double rangeStart) { return fmod(angle + 16 * PI - rangeStart, TWOPI) + rangeStart; }
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_ANGLEUTILS_H_INCLUDED
