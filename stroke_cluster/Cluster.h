#pragma once

#include <vector>
#include <map>
#include <ostream>
#include <glm/glm.hpp>

struct Cluster {
	struct Stroke {
		std::vector<glm::dvec2> points;
		std::vector<double> u;
	};

	/*struct XSecPoint {

	};

	struct XSec {

	};*/

	std::vector<Stroke> strokes;
	//std::vector<XSec> xsecs;
};

glm::dvec2 tangent(const std::vector<glm::dvec2>& points, size_t i);
glm::dvec2 normal(const glm::dvec2& v);

struct Input {
	double thickness;
	double width;
	double height;
	std::map<int, Cluster> clusters;

	void param_svg(std::ostream& os) const;
	void orientation_svg(std::ostream& os) const;
};