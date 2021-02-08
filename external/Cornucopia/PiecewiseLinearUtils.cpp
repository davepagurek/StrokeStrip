/*--
    PiecewiseLinearUtils.cpp  

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

#include "PiecewiseLinearUtils.h"
#include <limits>

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

void PiecewiseLinearMonotone::add(double x, double y)
{
    const double tol = 1e-8;

    y *= _sign;

    set<PLPoint>::const_iterator it = _points.insert(PLPoint(x, y)).first;
    set<PLPoint>::const_iterator it2 = it;

    //self test
    if(it != _points.begin() && y + tol < (--it)->y)
        Debugging::get()->printf("ERROR: Not monotone w.r.t. prev!");
    if(++it2 != _points.end() && y - tol > it2->y)
        Debugging::get()->printf("ERROR: Not monotone w.r.t. next!");        
}

bool PiecewiseLinearMonotone::eval(double x, double &outY) const
{
    const double tol = 1e-8;

    if(_points.empty())
        return false;

    //check if we're near the ends of the range
    if(fabs(x - _points.begin()->x) < tol)
    {
        outY = _sign * _points.begin()->y;
        return true;
    }
    set<PLPoint>::const_iterator it = _points.end();
    --it; //last element
    if(fabs(x - it->x) < tol)
    {
        outY = _sign * it->y;
        return true;
    }

    //now look for the point for real
    it = _points.lower_bound(PLPoint(x, 0.));

    //check if we're out of range
    if(it == _points.begin())
    {
        outY = _sign * it->y;
        return false;
    }
    if(it == _points.end())
    {
        outY = _sign * (--it)->y;
        return false;
    }

    set<PLPoint>::const_iterator prev = it;
    --prev;

    if(it->x - prev->x < tol) //if the slope is infinite or somehow negative
    {
        outY = _sign * prev->y;
        return true;
    }

    outY = prev->y + (it->y - prev->y) * (x - prev->x) / (it->x - prev->x);
    outY *= _sign;
    return true;
}

bool PiecewiseLinearMonotone::invert(double y, double &outX) const
{
    y *= _sign;

    const double tol = 1e-8;

    if(_points.empty())
        return false;

    //check if we're near the ends of the range
    if(fabs(y - _points.begin()->y) < tol)
    {
        outX = _points.begin()->x;
        return true;
    }
    set<PLPoint>::const_iterator it = _points.end();
    --it; //last element
    if(fabs(y - it->y) < tol)
    {
        outX = it->x;
        return true;
    }

    //now look for the point for real
    it = _points.lower_bound(PLPoint(y));

    //check if we're out of range
    if(it == _points.begin())
    {
        outX = it->x;
        return false;
    }
    if(it == _points.end())
    {
        outX = (--it)->x;
        return false;
    }

    set<PLPoint>::const_iterator prev = it;
    --prev;

    if(it->y - prev->y < tol) //if the slope is small or somehow negative
    {
        outX = prev->x;
        return true;
    }

    outX = prev->x + (it->x - prev->x) * (y - prev->y) / (it->y - prev->y);
    return true;
}

double PiecewiseLinearMonotone::minX() const
{
    if(_points.empty())
        return numeric_limits<double>::max();
    return _points.begin()->x;
}

double PiecewiseLinearMonotone::maxX() const
{
    if(_points.empty())
        return -numeric_limits<double>::max();
    return (--_points.end())->x;
}

bool PiecewiseLinearMonotone::batchEval(vector<double> &inXoutY) const
{
    const double tol = 1e-8;

    if(_points.empty())
        return false;

    bool allGood = true;
    for(int i = 0; i < (int)inXoutY.size(); ++i)
    {
        if(!eval(inXoutY[i], inXoutY[i]))
        {
            Debugging::get()->printf("PiecewiseLinearMonotone evaluation error!");
            allGood = false;
        }
    }
    return allGood;
}

END_NAMESPACE_Cornu


