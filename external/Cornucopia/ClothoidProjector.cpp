/*--
    ClothoidProjector.cpp  

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

#include "Clothoid.h"
#include "Arc.h"
#include "Fresnel.h"

#include <deque>

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

//an arc that approximates a part of a clothoid
class _ApproxArc
{
public:
    _ApproxArc(double start, double length)
        : _start(start), _length(length)
    {
        //build an arc through three points on the clothoid to approximate it
        Vector2d p[3];
        for(int i = 0; i < 3; ++i)
            fresnelApprox(start + 0.5 * length * i, &(p[i][1]), &(p[i][0]));
        _arc = new Arc(p[0], p[1], p[2]);
    }

    bool test(const Vector2d &pt, double &minDistSq, double &minT, double from, double to) const
    {
        //Do a quick-reject test:
        //Let r be the arc radius and let d be the distance from the point to the center.
        //We want to reject if (r - d)^2 > minDistSq.
        //To avoid division and square roots (we don't have d, only d^2), we write this as:
        //((r^2 - d^2)/(r + d))^2 > minDistSq, or (r^2 - d^2)^2 > minDistSq (r + d)^2
        //Since (r + d)^2 < 2 * (r^2 + d^2), we can conservatively check (r^2 - d^2)^2 > 2 * minDistSq * (r^2 + d^2) 
        double dSq = (pt - _arc->center()).squaredNorm();
        double rSq = SQR(_arc->radius());
        if(SQR(dSq - rSq) > 2. * minDistSq * (dSq + rSq))
            return false;

        //Now really check
        double t = _arc->project(pt);
        t = min(max(t + _start, from), to) - _start;

        double distSq = (pt - _arc->pos(t)).squaredNorm();
        if(distSq >= minDistSq)
            return false;

        minT = _start + t;
        minDistSq = distSq;
        return true;
    }

private:
    ArcPtr _arc;
    double _start;
    double _length;
};

class Clothoid::_ClothoidProjectorImpl : public Clothoid::_ClothoidProjector
{
public:
    typedef Vector2d Vec;

    _ClothoidProjectorImpl() //initialize
        : _arcSpacing(0.1)
    {
        double t;
        for(t = 0; t * _arcSpacing < 0.9; t += _arcSpacing)
        {
            _arcs.push_front(_ApproxArc(-t - _arcSpacing, _arcSpacing));
            _arcs.push_back(_ApproxArc(t, _arcSpacing));
        }
        _maxArcParam = t;
    }

    double project(const Vec &pt, double from, double to) const
    {
        double minT, minDistSq;

        Vector2d startPt, endPt;
        fresnelApprox(from, &(startPt[1]), &(startPt[0]));
        fresnelApprox(to, &(endPt[1]), &(endPt[0]));

        //test start and end points
        minT = from;
        minDistSq = (pt - startPt).squaredNorm();
        double distSq = (pt - endPt).squaredNorm();
        if(distSq < minDistSq)
        {
            minDistSq = distSq;
            minT = to;
        }

        //A possible optimization for clothoids outside the precomputed range is to batch
        //the fresnel integral evaluations.  Not sure how much it'll help.
        int minArcIdx = (int)floor((_maxArcParam + from) / _arcSpacing);
        if(minArcIdx < 0)
        {
            minArcIdx = 0;
            //test against arcs before the precomputed ones
            double start = from;
            double stop = min(to, -_maxArcParam);
            int cnt = 0; //iteration count to prevent looping over a really spirally clothoid
            while(start + 1e-8 < stop && ++cnt < 100)
            {
                double len = min(stop - start, -1. / start);
                _ApproxArc(start, len).test(pt, minDistSq, minT, from, to);

                start += len;
            }
        }
        int maxArcIdx = (int)ceil((_maxArcParam + to) / _arcSpacing);
        if(maxArcIdx > (int)_arcs.size())
        {
            maxArcIdx = (int)_arcs.size();
            //test against arcs past the end of the precomputed ones
            double start = to;
            double stop = max(_maxArcParam, from);
            int cnt = 0; //iteration count to prevent looping over a really spirally clothoid
            while(start - 1e-8 > stop && ++cnt < 100)
            {
                double len = min(start - stop, 1. / start);
                _ApproxArc(start - len, len).test(pt, minDistSq, minT, from, to);

                start -= len;
            }
        }

        for(int i = minArcIdx; i < maxArcIdx; ++i)
            _arcs[i].test(pt, minDistSq, minT, from, to);

        minT = projectNewton(minT, pt, from, to);
        minT = projectNewton(minT, pt, from, to);
        return minT;
    }

private:
    double projectNewton(double guess, const Vec &pt, double from, double to) const
    {
        Vec p, der, der2;
        eval(guess, &p, &der, &der2);

        double dot = der.dot(pt - p);
        double dotDer = der2.dot(pt - p) - der.squaredNorm();

        if(dotDer >= -1e-30) //if distance will increase
            return guess;

        return max(from, min(to, guess - dot / dotDer));
    }

    void eval(double t, Vec *pos, Vec *der = NULL, Vec *der2 = NULL) const
    {
        if(pos)
            fresnelApprox(t, &((*pos)[1]), &((*pos)[0]));
        if(der || der2)
        {
            double s = sin(HALFPI * t * t), c = cos(HALFPI * t * t);
            if(der)
                *der = Vec(c, s);
            if(der2)
                *der2 = (PI * t) * Vec(-s, c);
        }
    }

    const double _arcSpacing;
    deque<_ApproxArc> _arcs;
    double _maxArcParam;
};

Clothoid::_ClothoidProjector *Clothoid::_clothoidProjector()
{
    static Clothoid::_ClothoidProjector *projector = new _ClothoidProjectorImpl();
    return projector;
}

END_NAMESPACE_Cornu


