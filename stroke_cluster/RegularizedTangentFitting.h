#pragma once
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "MLSFitting.h"

double compute_edge_alignment(PointSample source, PointSample target);

std::vector<glm::dvec2> ls_solve(std::vector<glm::dvec2> const &samples,
				   std::vector<glm::dvec2> const &tangents, std::vector<double> const &weights, 
				   std::vector<double> const &pos_weights, 
				   double gamma, bool is_closed = false);