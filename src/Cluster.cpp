#define _USE_MATH_DEFINES
#include <cmath>

#include "Cluster.h"
#include "SvgUtils.h"
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

std::vector<std::string> colors = {
	"#9E0142",
	"#D53E4F",
	"#F46D43",
	"#FDAE61",
	"#FEE08B",
	"#ABDDA4",
	"#66C2A5",
	"#294dcc",
	"#3288BD",
	"#5E4FA2",
	"#8DD3C7",
	"#BEBADA",
	"#FB8072",
	"#80B1D3",
	"#FDB462",
	"#B3DE69",
	"#053061",
	"#BC80BD",
	"#FFED6F",
	"#8E0152",
	"#DE77AE",
	"#5AAE61",
	"#8AEF91",
	"#00441B",
	"#FB9A99",
	"#E31A1C",
	"#FDBF6F",
	"#FF7F00",
	"#8a20cf",
	"#C51B7D",
	"#8C510A"
};

glm::dvec2 simple_discrete_tangent(const std::vector<glm::dvec2>& points, size_t i) {
	glm::dvec2 result(0.0, 0.0);
	if (i > 0) {
		result += glm::normalize(points[i] - points[i - 1]);
	} else {
		result += glm::normalize(points[i + 1] - points[i]);
	}
	return glm::normalize(result);
}

glm::dvec2 discrete_tangent(const std::vector<glm::dvec2>& points, size_t i) {
	glm::dvec2 result(0.0, 0.0);
	if (i > 0) {
		result += glm::normalize(points[i] - points[i - 1]);
	}
	if (i < points.size() - 1) {
		result += glm::normalize(points[i + 1] - points[i]);
	}
	return glm::normalize(result);
}

glm::dvec2 tangent(const std::vector<glm::dvec2>& points, double i) {
	glm::dvec2 tangent(0.0, 0.0);
	size_t a = std::floor(i);
	if (double(a) == i) {
		tangent = discrete_tangent(points, (size_t)i);
	}
	else {
		size_t b = std::ceil(i);
		double mix = i - double(a);
		tangent = glm::normalize((1 - mix) * discrete_tangent(points, a) + mix * discrete_tangent(points, b));
	}

	if (!glm::any(glm::isnan(tangent))) {
		return tangent;
	}
	return simple_discrete_tangent(points, (size_t)std::floor(i));
}

glm::dvec2 point(const std::vector<glm::dvec2>& points, double i) {
	size_t a = std::floor(i);
	size_t b = std::ceil(i);
	double mix = i - double(a);
	return (1 - mix) * points[a] + mix * points[b];
}

glm::dvec2 normal(const glm::dvec2& v) {
	return glm::dvec2(-v.y, v.x);
}

double total_length(const std::vector<glm::dvec2>& points) {
	double len = 0;
	for (size_t i = 1; i < points.size(); ++i) {
		len += glm::distance(points[i - 1], points[i]);
	}
	return len;
}

glm::dvec2 midpoint(const std::vector<glm::dvec2>& points) {
	return points[points.size() / 2];
}

void Input::param_svg(std::ostream& os, bool rainbow) const {
	double padding = thickness;
	double w = width * thickness + 2 * padding;
	double h = height * thickness + 2 * padding;
	double x = -w / 2;
	double y = -h / 2;

	SVG::begin(os, x, y, w, h);
	for (auto& kv : clusters) {
		double max_u = kv.second.max_u();

		for (auto& stroke : kv.second.strokes) {
			std::vector<glm::dvec2> scaled;
			scaled.reserve(stroke.points.size());
			for (auto& pt : stroke.points) {
				scaled.push_back(pt * thickness);
			}
			for (size_t i = 0; i < stroke.points.size() - 1; ++i) {
				std::stringstream ss;
				if (rainbow) {
					ss << "hsl(";
					ss << int(stroke.u[i] / max_u * 360) << ", ";
					ss << "90%, 50%)";
				}
				else {
					ss << "#";
					ss << std::setfill('0') << std::setw(2);
					ss << std::hex << int(stroke.u[i] / max_u * 255); // red
					ss << "00"; // green
					ss << std::setfill('0') << std::setw(2);
					ss << std::hex << int(255 - stroke.u[i] / max_u * 255); // blue
				}
				SVG::line(os, scaled[i].x, scaled[i].y, scaled[i + 1].x, scaled[i + 1].y, thickness, ss.str());
			}
		}
	}
	SVG::end(os);
}

void Input::orientation_svg(std::ostream& os, std::function<void(std::ostream&)> cb) const {
	double padding = thickness;
	double w = width * thickness + 2 * padding;
	double h = height * thickness + 2 * padding;
	double x = -w / 2;
	double y = -h / 2;

	SVG::begin(os, x, y, w, h);
	size_t stroke_num = 0;
	for (auto& kv : clusters) {
		for (auto& stroke : kv.second.strokes) {
			std::vector<glm::dvec2> scaled;
			scaled.reserve(stroke.points.size());
			for (auto& pt : stroke.points) {
				scaled.push_back(pt * thickness);
			}
			std::string color(colors[stroke_num % colors.size()]);
			SVG::polyline(os, scaled, thickness, color);

			for (size_t i = 2; i < stroke.points.size(); i += 5) {
				glm::dvec2 origin = scaled[i];
				glm::dvec2 next = origin + normal(tangent(scaled, i)) * 4. * thickness;
				SVG::line(os, origin.x, origin.y, next.x, next.y, thickness, color);
			}
			++stroke_num;
		}
	}
	cb(os);
	SVG::end(os);
}

void Input::cluster_svg(std::ostream& os, std::function<void(std::ostream&)> cb) const {
	double padding = thickness;
	double w = width * thickness + 2 * padding;
	double h = height * thickness + 2 * padding;
	double x = -w / 2;
	double y = -h / 2;

	SVG::begin(os, x, y, w, h);
	size_t cluster_num = 0;
	for (auto& kv : clusters) {
		std::string color(colors[cluster_num % colors.size()]);
		for (auto& stroke : kv.second.strokes) {
			std::vector<glm::dvec2> scaled;
			scaled.reserve(stroke.points.size());
			for (auto& pt : stroke.points) {
				scaled.push_back(pt * thickness);
			}
			SVG::polyline(os, scaled, thickness, color);
		}
		++cluster_num;
	}
	cb(os);
	SVG::end(os);
}

double Cluster::XSec::distance_weight(size_t i) const {
	if (points.size() == 1) return 1.0;

	double weight = 0.0;
	if (i == 0) {
		weight += glm::distance(points[0].point, points[1].point);
	}
	else {
		weight += glm::distance(points[i - 1].point, points[i].point);
	}

	if (i == points.size() - 1) {
		weight += glm::distance(points[points.size() - 2].point, points[points.size() - 1].point);
	}
	else {
		weight += glm::distance(points[i].point, points[i + 1].point);
	}

	return weight;
}

glm::dvec2 Cluster::XSec::avg_point() const {
	glm::dvec2 sum(0.0, 0.0);
	for (auto& p : points) {
		sum += p.point;
	}
	if (points.size() > 0) sum /= double(points.size());
	return sum;
}

glm::dvec2 Cluster::XSec::avg_tangent() const {
	if (points.size() == 1) return points.front().tangent;

	glm::dvec2 sum(0.0, 0.0);

	std::vector<double> dists(points.size(), 0.0);
	for (size_t i = 0; i < points.size() - 1; ++i) {
		dists[i] = glm::distance(points[i].point, points[i + 1].point);
	}
	for (size_t i = 0; i < points.size(); ++i) {
		double weight = 0.0;
		if (i == 0) {
			weight += dists[0];
		}
		else {
			weight += dists[i - 1];
		}
		if (i == points.size() - 1) {
			weight += dists.back();
		}
		else {
			weight += dists[i + 1];
		}
		weight += 0.1; // Regularize to avoid zeros

		sum += weight * points[i].tangent;
	}

	sum = glm::normalize(sum);
	if (glm::any(glm::isnan(sum))) {
		sum = glm::dvec2();
	}
	return sum;
}

double Cluster::max_u() const {
	double result = -std::numeric_limits<double>::infinity();
	for (auto& stroke : strokes) {
		for (double u : stroke.u) {
			result = std::max(result, u);
		}
	}
	return result;
}