#pragma once

#include <vector>
#include <glm/glm.hpp>
#include <gurobi_c++.h>

struct Intersection {
	double i;
	glm::dvec2 pt;
};
std::vector<Intersection> intersections(const std::vector<glm::dvec2>& polyline, glm::dvec2 from, glm::dvec2 to);

template<typename T> T map(T in, T in_from, T in_to, T out_from, T out_to) {
	T mixed = (in - in_from) / (in_to - in_from);
	T res = mixed * (out_to - out_from) + out_from;
	if (res < std::min(out_from, out_to)) res = std::min(out_from, out_to);
	if (res > std::max(out_from, out_to)) res = std::max(out_from, out_to);
	return res;
}

double gaussian(double x, double sigma, double mu);

double poly_in_out(double t, double k);

GRBLinExpr l1_norm(GRBModel* model, const std::vector<GRBLinExpr>& x);
GRBQuadExpr l2_norm_sq(GRBModel* model, const std::vector<GRBLinExpr>& x);