#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <ostream>
#include <string>

namespace SVG {
	inline void begin(std::ostream& os, double x, double y, double w, double h) {
		os << "<svg viewBox = '" << x << " " << y << " " << w << " " << h << "' xmlns = 'http://www.w3.org/2000/svg'>" << std::endl;
	}

	inline void line(std::ostream& os, double x1, double y1, double x2, double y2, double w = 1., const std::string& color = "#000") {
		os << "<line stroke-width='" << w << "' x1='" << x1 << "' y1='" << y1 << "' x2='" << x2 << "' y2='" << y2 << "' fill='none' stroke='" << color << "' />" << std::endl;
	}

	inline void polyline(std::ostream& os, const std::vector<glm::dvec2>& pts, double w = 1., const std::string& color = "#000") {
		os << "<path stroke-width='" << w << "' fill='none' stroke='" << color << "' d='";
		for (size_t i = 0; i < pts.size(); ++i) {
			if (i == 0) {
				os << "M ";
			}
			else {
				os << " L ";
			}
			os << pts[i].x << " " << pts[i].y;
		}
		os << "' />";
	}

	inline void variable_polyline(std::ostream& os, std::vector<glm::dvec2>& pts, std::vector<double> w, const std::string& color = "#000") {
		std::vector<glm::dvec2> normals;
		for (size_t i = 0; i < pts.size()-1; ++i) {
			normals.push_back(normal(glm::normalize(pts[i+1] - pts[i])));
		}
		normals.push_back(normals.back());

		std::vector<glm::dvec2> bottom, top;
		for (size_t i = 0; i < pts.size()-1; ++i) {
			top.push_back(pts[i] - normals[i] * w[i] / 2.0);
			bottom.push_back(pts[i] + normals[i] * w[i] / 2.0);
		}
		std::reverse(bottom.begin(), bottom.end());

		std::vector<glm::dvec2> cap_start, cap_end;
		if (pts.front() != pts.back()) {
			auto rotate = [](glm::dvec2 vec, double theta) -> glm::dvec2 {
				glm::mat2 rotation = {
					{std::cos(theta), std::sin(theta)},
					{-std::sin(theta), std::cos(theta)}
				};
				return rotation*vec;
			};

			for (double t = 0.1; t <= 0.9; t += 0.1) {
				cap_start.push_back(pts.front() + rotate(normals.front(), t * M_PI + M_PI) * w.front()/2.0);

				cap_end.push_back(pts.back() + rotate(normals.back(), t * M_PI) * w.back()/2.0);
			}
		}

		std::vector<glm::dvec2> all_pts;
		for (auto pt : cap_start) all_pts.push_back(pt);
		for (auto pt : top) all_pts.push_back(pt);
		for (auto pt : cap_end) all_pts.push_back(pt);
		for (auto pt : bottom) all_pts.push_back(pt);

		os << "<path stroke='none' fill='" << color << "' d='";
		for (size_t i = 0; i < all_pts.size(); ++i) {
			if (i == 0) {
				os << "M ";
			}
			else {
				os << " L ";
			}
			os << all_pts[i].x << " " << all_pts[i].y;
		}
		os << " Z' />";
	}

	inline void end(std::ostream& os) {
		os << "</svg>" << std::endl;
	}
}
