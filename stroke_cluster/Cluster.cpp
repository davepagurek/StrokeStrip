#include "Cluster.h"
#include "SvgUtils.h"
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>

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

glm::dvec2 discrete_tangent(const std::vector<glm::dvec2>& points, size_t i) {
	glm::dvec2 result;
	if (i > 0) {
		result += glm::normalize(points[i] - points[i - 1]);
	}
	if (i < points.size() - 1) {
		result += glm::normalize(points[i + 1] - points[i]);
	}
	return glm::normalize(result);
}

glm::dvec2 tangent(const std::vector<glm::dvec2>& points, double i) {
	size_t a = std::floor(i);
	if (double(a) == i) return discrete_tangent(points, (size_t)i);

	size_t b = std::ceil(i);
	double mix = i - double(a);
	return glm::normalize((1 - mix) * discrete_tangent(points, a) + mix * discrete_tangent(points, b));
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

void Input::param_svg(std::ostream& os) const {
	double padding = thickness;
	double w = width * thickness + 2 * padding;
	double h = height * thickness + 2 * padding;
	double x = -w / 2;
	double y = -h / 2;

	SVG::begin(os, x, y, w, h);
	for (auto& kv : clusters) {
		for (auto& stroke : kv.second.strokes) {
			std::vector<glm::dvec2> scaled;
			scaled.reserve(stroke.points.size());
			for (auto& pt : stroke.points) {
				scaled.push_back(pt * thickness);
			}
			for (size_t i = 0; i < stroke.points.size() - 1; ++i) {
				std::stringstream ss;
				ss << "#";
				ss << std::setfill('0') << std::setw(2);
				ss << std::hex << int(stroke.u[i] * 255); // red
				ss << "00"; // green
				ss << std::hex << int(255 - stroke.u[i] * 255); // blue
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

			for (size_t i = 5; i < stroke.points.size(); i += 5) {
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

