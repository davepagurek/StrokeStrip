#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>
#include <limits>

#include <algorithm>
#include <iostream>
#include <deque>

#include "StrokeCutting.h"
#include "StrokeOrientation.h"
#include "Cluster.h"
#include "Parameterization.h"
#include "Fitting.h"

void to_scap(std::string const &filename, Capture const &capture, int width, int height) {
	std::ofstream scap_ofs(filename);
	std::string out_buffer;
	out_buffer += capture.to_string();

	scap_ofs << "#" << width << "\t" << height << std::endl;
	scap_ofs.write(out_buffer.c_str(), out_buffer.size());
	scap_ofs.close();
}

Input from_capture(Capture capture) {
	Input input;
	auto& clusters = input.clusters;


	double min_x = std::numeric_limits<double>::infinity();
	double max_x = -std::numeric_limits<double>::infinity();
	double min_y = std::numeric_limits<double>::infinity();
	double max_y = -std::numeric_limits<double>::infinity();

	input.thickness = capture.thickness;

	double rate = 2.75;

	for (auto& polyline : capture.sketchedPolylines) {
		polyline.reparameterize(std::min(
			rate * capture.thickness,
			polyline.totalLen() / 3));

		for (auto& point : polyline.points) {
			min_x = std::min(min_x, point.first.x);
			max_x = std::max(max_x, point.first.x);
			min_y = std::min(min_y, point.first.y);
			max_y = std::max(max_y, point.first.y);
		}
	}
	glm::dvec2 center((max_x + min_x) / 2, (max_y + min_y) / 2);
	input.width = (max_x - min_x) / capture.thickness;
	input.height = (max_y - min_y) / capture.thickness;

	for (auto& polyline : capture.sketchedPolylines) {
		if (polyline.points.empty()) continue;
		clusters[polyline.group_ind].strokes.emplace_back();
		auto& stroke = clusters[polyline.group_ind].strokes.back();
		stroke.points.reserve(polyline.points.size());
		stroke.u.reserve(polyline.points.size());
		for (size_t i = 0; i < polyline.points.size(); ++i) {
			auto& point = polyline.points[i];

			// Recenter and normalize to stroke width
			stroke.points.push_back((glm::dvec2(point.first.x, point.first.y) - center) / capture.thickness);
			stroke.u.push_back(0);
		}
	}

	return input;
}

int main(int argc, char** argv) {
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " input.scap [args]" << std::endl;
		std::cout << "\t-d, --debug:   Write debug visuals" << std::endl;
		std::cout << "\t-c, --cut:     Cut sharp turns before processing" << std::endl;
		std::cout << "\t-r, --rainbow: Visualize parameters as a rainbow instead of red-blue" << std::endl;
		std::cout << "\t-w, --widths: Fit widths to clusters" << std::endl;
		std::cout << "\t-t, --taper: Force fitted widths to taper at ends" << std::endl;
		return -1;
	}
	std::string scap_filename(argv[1]);

	Context context;

	std::deque<std::string> args;
	for (size_t i = 2; i < argc; ++i) args.push_back(argv[i]);

	while (!args.empty()) {
		std::string arg = args.front();
		args.pop_front();

		if (arg == "-d" || arg == "--debug") {
			context.debug_viz = true;
		}
		else if (arg == "-c" || arg == "--cut") {
			context.cut = true;
		} else if (arg == "-r" || arg == "--rainbow") {
			context.rainbow = true;
		} else if (arg == "-w" || arg == "--widths") {
			context.widths = true;
		} else if (arg == "-t" || arg == "--taper") {
			context.taper_widths = true;
		}
	}

	/*std::string scap_filename = "D:\\strokestrip\\tests\\siggraph.scap";
	Context context;
	context.debug_viz = true;*/

	Input input;

	// 1. Preprocess
	{
		std::ifstream scap_ifs(scap_filename);
		std::stringstream buffer;
		buffer << scap_ifs.rdbuf();
		scap_ifs.close();

		Capture capture;
		capture.from_string(buffer.str());

		std::string canvas_size_line = buffer.str().substr(0, buffer.str().find_first_of("\n"));
		int width;
		int height;
		sscanf(canvas_size_line.c_str(), "#%d\t%d", &width, &height);

		// Assign missing ids
		int max_ind = -1;
		for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
			max_ind = std::max(max_ind, capture.getSketchedPolyline(i).stroke_ind);
		}

		if (context.cut) {
			preprocess_cluster(1, width, height, &capture);
		}
		input = from_capture(capture);
	}

	// 2. Orientation
	{
		StrokeOrientation orientation(context);
		orientation.orient_strokes(&input);
		if (context.debug_viz) {
			std::string final_output_name = scap_filename;
			final_output_name.erase(final_output_name.length() - 5, 5); // remove .scap
			final_output_name += "_orientation_debug.svg";
			std::ofstream orientation_svg(final_output_name);
			orientation.orientation_debug(orientation_svg, input);
		}
		orientation.flip_strokes(&input);
	}

	// 3. Parameterization
	{
		Parameterization param(context);
		param.parameterize(&input);
		if (context.debug_viz) {
			std::string final_output_name = scap_filename;
			final_output_name.erase(final_output_name.length() - 5, 5); // remove .scap
			final_output_name += "_param_debug.svg";
			std::ofstream debug_svg(final_output_name);
			param.debug_svg(debug_svg, input);
		}
		{
			std::string final_output_name = scap_filename;
			final_output_name.erase(final_output_name.length() - 5, 5); // remove .scap
			final_output_name += "_isolines.svg";
			std::ofstream isolines_svg(final_output_name);
			param.isolines_svg(isolines_svg, input);
		}
	}

	// 4. Fitting
	{
		Fitting fitting(context);
		auto fits = fitting.fit(&input);
		{
			std::string final_output_name = scap_filename;
			final_output_name.erase(final_output_name.length() - 5, 5); // remove .scap
			final_output_name += "_fit.svg";
			std::ofstream fit_svg(final_output_name);
			fitting.fit_svg(fit_svg, input, fits);
		}
	}

	{
		std::string final_output_name = scap_filename;
		final_output_name.erase(final_output_name.length() - 5, 5); // remove .scap
		final_output_name += "_params.svg";
		std::ofstream param_svg(final_output_name);
		input.param_svg(param_svg, context.rainbow);
	}

	{
		std::string final_output_name = scap_filename;
		final_output_name.erase(final_output_name.length() - 5, 5); // remove .scap
		final_output_name += "_orientation.svg";
		std::ofstream orientation_svg(final_output_name);
		input.orientation_svg(orientation_svg);
	}

	/*{
		Capture result = to_capture(input);
		to_scap(final_output_name, result, input.width * input.thickness, input.height * input.thickness);
	}*/

	return 0;
}
