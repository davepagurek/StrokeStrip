#include "DrawPrimitive.h"
#include "Renderer.h"

extern Renderer renderer;

void DrawDot::draw() const {
	renderer.draw_dot(dot, color);
}

std::string DrawDot::to_string() const {
	std::string str = "";

	str += "{";
	str += "\t#-2\t-2\n";

	str += "\t" + std::to_string(dot.x) + "\t" + std::to_string(dot.y) + "\t0\n";

	str += "}\n";

	return str;
}

void DrawCircle::draw() const {
	renderer.draw_circle(glm::dvec3(center.x, center.y, radius), color);
}

std::string DrawCircle::to_string() const {
	std::string str = "";

	str += "{";
	str += "\t#-3\t-3\n";

	str += "\t" + std::to_string(center.x) + "\t" + std::to_string(center.y) + "\t0\n";
	str += "\t" + std::to_string(radius) + "\t" + std::to_string(radius) + "\t0\n";

	str += "}\n";

	return str;
}

void DrawSegment::draw() const {
	renderer.draw_segment(std::make_pair(source, destination), color);
}

std::string DrawSegment::to_string() const {
	std::string str = "";

	str += "{";
	str += "\t#-4\t-4\n";

	str += "\t" + std::to_string(source.x) + "\t" + std::to_string(source.y) + "\t0\n";
	str += "\t" + std::to_string(destination.x) + "\t" + std::to_string(destination.y) + "\t0\n";

	str += "}\n";

	return str;
}

void DrawPath::draw() const {
	Sketch stroke;
	for (auto const &p : points) {
		stroke.points.emplace_back(SketchPoint(p.x, p.y), 0);
	}

	renderer.draw_stroke(stroke, color);
}

std::string DrawPath::to_string() const {
	std::string str = "";

	Sketch stroke;
	for (auto const &p : points) {
		stroke.points.emplace_back(SketchPoint(p.x, p.y), 0);
	}
	stroke.stroke_ind = -1;
	stroke.group_ind = -1;

	str += stroke.to_string();

	return str;
}