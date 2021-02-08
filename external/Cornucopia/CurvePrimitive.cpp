/*--
    CurvePrimitive.cpp  

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

#include "CurvePrimitive.h"

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

bool CurvePrimitive::isValid() const
{
    if((int)_params.size() != numParams())
        return false;

    for(int i = 0; i < (int)_params.size(); ++i)
        if(!(_params[i] == _params[i])) //check for NaN
            return false;

    return isValidImpl();
}

END_NAMESPACE_Cornu


