/*--
    Debugging.h  

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

#ifndef CORNUCOPIA_DEBUGGING_H_INCLUDED
#define CORNUCOPIA_DEBUGGING_H_INCLUDED

#include "defs.h"
#include "smart_ptr.h"

#include <cstdarg>
#include <string>
#include <Eigen/Core>

NAMESPACE_Cornu

CORNU_SMART_FORW_DECL(Curve);
CORNU_SMART_FORW_DECL(CurvePrimitive);

//This class does nothing by default
class Debugging
{
public:
    typedef Eigen::Vector3d Color;
    typedef Eigen::Vector2d Vector2d;

    enum LineStyle
    {
        SOLID,
        DASHED,
        DOTTED
    };

    static Debugging *get() { return _currentDebugging; }

    virtual ~Debugging() {}

    virtual bool isDebuggingOn() const { return false; }

    virtual void printf(const char *fmt, ...)
    {
		return;

        const int sz = 200; //we don't need terribly long debug strings
        char buffer[sz];

        va_list ap;
        va_start(ap, fmt);

        vsnprintf(buffer, sz, fmt, ap);

        va_end(ap);

        std::printf("%s\n", buffer);
    }

    virtual void startTiming(const std::string &/*description*/) {}
    virtual void elapsedTime(const std::string &/*description*/) {} //prints the elapsed time
    virtual double getTimeElapsed(const std::string &/*description*/) { return 0.; } //in seconds

    virtual void clear(const std::string &/*groups*/ = "") {}

    virtual void drawPoint(const Vector2d &/*pos*/, const Color &/*color*/, const std::string &/*group*/ = "") {}
    virtual void drawLine(const Vector2d &/*p1*/, const Vector2d &/*p2*/, const Color &/*color*/, const std::string &/*group*/ = "", double /*thickness*/ = 1, LineStyle /*style*/ = SOLID) {}
    virtual void drawCurve(CurveConstPtr /*curve*/, const Color &/*color*/, const std::string &/*group*/ = "", double /*thickness*/ = 1, LineStyle /*style*/ = SOLID) {}
    virtual void drawCurvatureField(CurveConstPtr /*curve*/, const Color &/*color*/, const std::string &/*group*/ = "", double /*thickness*/ = 1, LineStyle /*style*/ = SOLID) {}

    //utility function that picks the right color
    void drawPrimitive(CurvePrimitiveConstPtr curve, const std::string &group, int idx = 0, double thickness = 1.);

protected:
    Debugging() {}
    static void set(Debugging *debugging);

private:
    static Debugging *_currentDebugging;
};

END_NAMESPACE_Cornu

#endif //CORNUCOPIA_DEBUGGING_H_INCLUDED
