#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <future>

#include "Utils.h"


std::vector<Intersection> intersections(const std::vector<glm::dvec2>& polyline, glm::dvec2 from, glm::dvec2 to) {
	bool vertical = std::abs(from.x - to.x) < 1e-4;
	std::vector<double> idx_b;
	std::vector<double> x;
	std::vector<double> y;
	const int numSegments = polyline.size() - 1;

	double x0a = from.x;
	double x1a = to.x;
	double y0a = from.y;
	double y1a = to.y;

	std::vector<std::future<Intersection>> futures;
	futures.reserve(numSegments);
	for (int s = 0; s < numSegments; s++) futures.push_back(std::async(std::launch::async, [=]() -> Intersection {
		double x0b = polyline[s].x;
		double x1b = polyline[s + 1].x;
		double y0b = polyline[s].y;
		double y1b = polyline[s + 1].y;

		bool other_vertical = std::abs(x0b - x0a) < 1e-4;

		// Ignore if outside bbox
		if (
			std::min(x0b, x1b) > std::max(x0a, x1a) ||
			std::max(x0b, x1b) < std::min(x0a, x1a) ||
			std::min(y0b, y1b) > std::max(y0a, y1a) ||
			std::max(y0b, y1b) < std::min(y0a, y1a)
			) return { -1, glm::dvec2() };

		double tb = 0;
		if (vertical) {
			if (other_vertical) {
				bool intersection = std::max(std::min(y0a, y1a), std::min(y0b, y1b)) < std::min(std::max(y0a, y1a), std::max(y0b, y1b));
				if (intersection) {
					tb = (std::max(std::min(y0a, y1a), std::min(y0b, y1b)) - y0b) / (y1b - y0b);
				}
				else {
					return { -1, glm::dvec2() };
				}
			}
			else {
				tb = (x0a - x0b) / (x1b - x0b);
				if (tb < 0 || tb > 1) return { -1, glm::dvec2() };
			}

		}
		else {
			double tb_num = y0a + (y1a - y0a)*(x0b - x0a) / (x1a - x0a) - y0b;
			double tb_denom = y1b - (y1a - y0a)*(x1b - x0b) / (x1a - x0a) - y0b;
			tb = tb_num / tb_denom;
			if (tb < 0 || tb > 1) return { -1, glm::dvec2() };

			double ta = (x0b + (x1b - x0b)*tb - x0a) / (x1a - x0a);
			if (ta < 0 || ta > 1) return { -1, glm::dvec2() };
		}

		return { static_cast<double>(s) + tb, glm::dvec2(x0b + tb * (x1b - x0b), y0b + tb * (y1b - y0b)) };
	}));

	std::vector<Intersection> intersections;
	for (auto& f : futures) {
		auto intersection = f.get();
		if (intersection.i != -1) {
			intersections.push_back(intersection);
		}
	}
	return intersections;
}

double gaussian(double x, double sigma, double mu) {
	return std::exp(-0.5 * std::pow((x - mu) / sigma, 2));
}

double poly_in_out(double t, double k) {
	double t2 = 2 * t;
	if (t2 <= 1) {
		return std::pow(t2, k) / 2.;
	}
	else {
		return (2. - std::pow(2. - t2, k)) / 2.;
	}
}

GRBLinExpr l1_norm(GRBModel* model, const std::vector<GRBLinExpr>& x) {
	// min ||x||_1
	//
	// ...is equivalent to:
	//
	// min t
	// s.t.
	// x_i <= y_i,
	// -x_i <= y_i,
	// \sum_i y_i = t

	GRBLinExpr sum_y = 0.0;
	auto t = model->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
	std::vector<GRBVar> y;
	y.reserve(x.size());
	for (auto& term : x) {
		y.push_back(model->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS));
		sum_y += y.back();
		model->addConstr(term <= y.back());
		model->addConstr(-term <= y.back());
	}
	model->addConstr(sum_y == t);

	return t;
}

GRBQuadExpr l2_norm_sq(GRBModel* model, const std::vector<GRBLinExpr>& x) {
	GRBQuadExpr result = 0.0;
	for (auto& term : x) {
		result += term * term;
	}

	return result;
}
