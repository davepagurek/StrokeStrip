#pragma once

#include <vector>
#include <map>
#include <ostream>
#include <functional>
#include <glm/glm.hpp>

struct Cluster {
	struct Stroke {
		std::vector<glm::dvec2> points;
		std::vector<double> u;
	};

	struct XSecPoint {
		size_t stroke_idx;
		double i;
		glm::dvec2 point;
		glm::dvec2 to_next;
		glm::dvec2 tangent;
	};

	struct XSecConnection {
		size_t a_idx;
		size_t b_idx;
		double weight;
	};

	struct XSec {
		std::vector<XSecPoint> points;
		std::vector<XSecConnection> connections;
		size_t center_idx;
		double u;
		bool connector;

		double distance_weight(size_t i) const;
		glm::dvec2 avg_tangent() const;
		glm::dvec2 avg_point() const;
	};

	std::vector<Stroke> strokes;
	std::vector<XSec> xsecs;
	bool periodic = false;

	double max_u() const;
};

glm::dvec2 point(const std::vector<glm::dvec2>& points, double i);
glm::dvec2 tangent(const std::vector<glm::dvec2>& points, double i);
glm::dvec2 discrete_tangent(const std::vector<glm::dvec2>& points, size_t i);
glm::dvec2 normal(const glm::dvec2& v);
glm::dvec2 midpoint(const std::vector<glm::dvec2>& points);
double total_length(const std::vector<glm::dvec2>& points);

struct Input {
	double thickness;
	double width;
	double height;
	std::map<int, Cluster> clusters;

	void param_svg(std::ostream& os, bool rainbow = false) const;
	void orientation_svg(std::ostream& os, std::function<void(std::ostream&)> cb = [](std::ostream&) { }) const;
	void cluster_svg(std::ostream& os, std::function<void(std::ostream&)> cb = [](std::ostream&) {}) const;
};