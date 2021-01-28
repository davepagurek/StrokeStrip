#pragma once
#include <ostream>
#include <string>

namespace SVG {
	void begin(std::ostream& os, double x, double y, double w, double h) {
		os << "<svg viewBox = '" << x << " " << y << " " << w << " " << h << "' xmlns = 'http://www.w3.org/2000/svg'>" << std::endl;
	}

	void line(std::ostream& os, double x1, double y1, double x2, double y2, double w = 1., const std::string& color = "#000") {
		os << "<line stroke-width='" << w << "' x1='" << x1 << "' y1='" << y1 << "' x2='" << x2 << "' y2='" << y2 << "' fill='none' stroke='" << color << "' />" << std::endl;
	}

	void polyline(std::ostream& os, const std::vector<glm::dvec2>& pts, double w = 1., const std::string& color = "#000") {
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

	void end(std::ostream& os) {
		os << "</svg>" << std::endl;
	}
}