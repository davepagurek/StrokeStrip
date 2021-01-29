#pragma once

#include <vector>
#include <glm/glm.hpp>

struct Intersection {
	double i;
	glm::dvec2 pt;
};
std::vector<Intersection> intersections(const std::vector<glm::dvec2>& polyline, glm::dvec2 from, glm::dvec2 to);