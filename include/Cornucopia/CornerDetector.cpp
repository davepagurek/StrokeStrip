/*--
    CornerDetector.cpp  

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

#include "CornerDetector.h"
#include "Fitter.h"
#include "Oversketcher.h"
#include "Polyline.h"
#include "PrimitiveFitter.h"

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

class BaseCornerDetector : public Algorithm<CORNER_DETECTION>
{
public:
    string name() const { return "Default"; }

protected:
    void _run(const Fitter &fitter, AlgorithmOutput<CORNER_DETECTION> &out)
    {
        const VectorC<Vector2d> &pts = fitter.output<OVERSKETCHING>()->output->pts();

        out.corners.resize(pts.size(), false);
        out.corners.setCircular(pts.circular());

        VectorC<double> scores = cornerScores(fitter);

        const double spacing = fitter.scaledParameter(Parameters::MINIMUM_CORNER_SPACING);
        const double threshold = fitter.params().get(Parameters::CORNER_THRESHOLD);

        for(int i = 0; i < pts.size(); ++i)
        {
            bool valid = true;
            if(pts.circular() == NOT_CIRCULAR && (i == 0 || i + 1 == pts.size()))
                valid = false;

            for(int sign = -1; sign < 2; sign += 2) //walk in both directions
            {
                VectorC<Vector2d>::Circulator cur = pts.circulator(i);
                double dist = 0;
                for(cur += sign; !cur.done(); cur += sign)
                {
                    dist += ((*cur) - (*(cur - sign))).norm();
                    if(dist > spacing) //if we walked too far along the line, break
                        break;

                    //check if we have a point near i with a higher corner probability
                    if(cur.index() != i && scores[cur.index()] >= scores[i])
                    {
                        if(cur.index() < i || scores[cur.index()] > scores[i]) //if scores are equal, break ties by index
                        {
                            valid = false; //if so, i is not a corner
                            break;
                        }
                    }
                }
            }

            out.corners[i] = valid && scores[i] > threshold;

            if(out.corners[i])
                Debugging::get()->drawPoint(pts[i], Vector3d(1, 0, 0), "Corners");
        }
    }

    virtual VectorC<double> cornerScores(const Fitter &fitter) = 0;
};

class DefaultCornerDetector : public BaseCornerDetector
{
protected:
    virtual VectorC<double> cornerScores(const Fitter &fitter)
    {
        PolylineConstPtr input = fitter.output<OVERSKETCHING>()->output;
        const VectorC<Vector2d> &pts = input->pts();
        VectorC<double> out(pts.size(), pts.circular());

        double offset = fitter.scaledParameter(Parameters::CORNER_NEIGHBORHOOD);

        for(int i = 0; i < pts.size(); ++i)
        {
            double param = input->idxToParam(i);
            double cornerParam = offset; //corner parameter on the trimmed curve
            if(!input->isClosed())
                cornerParam = min(param, offset);
            out[i] = cornerScore(fitter, input->trimmed(param - offset, param + offset), cornerParam);
        }

        return out;
    }

    double cornerScore(const Fitter &fitter, PolylineConstPtr cornerPoly, double cornerParam)
    {
        double step = fitter.scaledParameter(Parameters::DENSE_SAMPLING_STEP);

        VectorC<Vector2d> resampled(0, NOT_CIRCULAR);

        for(double param = cornerParam; param >= 0.; param -= step)
            resampled.insert(resampled.begin(), cornerPoly->pos(param));
        int cornerIdx = resampled.size() - 1;
        for(double param = cornerParam + step; param <= cornerPoly->length(); param += step)
            resampled.push_back(cornerPoly->pos(param));

        int smoothingSteps = (int)fitter.params().get(Parameters::CORNER_SCALES);
        vector<VectorC<Vector2d> > smoothed(smoothingSteps, resampled);

        //laplacian smooth
        for(int i = 1; i < (int)smoothed.size(); ++i)
            for(int j = 1; j + 1 < (int)resampled.size(); ++j)
                smoothed[i][j] = 0.25 * (smoothed[i - 1][j - 1] + smoothed[i - 1][j + 1] + 2. * smoothed[i - 1][j]);

        double maxAngle = -1e10, minAngle = 1e10;

        for(int i = 1; i < (int)smoothed.size(); ++i)
        {
            int offs = i + 1;
            if(cornerIdx - offs < 0)
                return 0.;
            if(cornerIdx + offs >= resampled.size())
                return 0.;

            Vector2d v1 = smoothed[i][cornerIdx] - smoothed[i][cornerIdx - offs];
            Vector2d v2 = smoothed[i][cornerIdx + offs] - smoothed[i][cornerIdx];
            double angle = AngleUtils::angle(v1, v2);
            maxAngle = max(maxAngle, angle);
            minAngle = min(minAngle, angle);
        }

        double worstAngle = min(fabs(minAngle), fabs(maxAngle));

        return worstAngle;
    }
};

void Algorithm<CORNER_DETECTION>::_initialize()
{
    new DefaultCornerDetector();
}

END_NAMESPACE_Cornu


