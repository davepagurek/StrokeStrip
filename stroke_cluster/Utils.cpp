#include "Utils.h"

#include <cmath>
#include <vector>
#include <future>


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

		// Ignore if outside bbox
		if (
			std::min(x0b, x1b) > std::max(x0a, x1a) ||
			std::max(x0b, x1b) < std::min(x0a, x1a) ||
			std::min(y0b, y1b) > std::max(y0a, y1a) ||
			std::max(y0b, y1b) < std::min(y0a, y1a)
			) return { -1, glm::dvec2() };

		double tb = 0;
		if (vertical) {
			tb = (x0a - x0b) / (x1b - x0b);
			if (tb < 0 || tb > 1) return { -1, glm::dvec2() };

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