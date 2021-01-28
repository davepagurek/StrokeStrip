#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>
#include <map>

#include <Windows.h>
#include <algorithm>

#include "StrokeCutting.h"

void to_scap(std::string const &filename, Capture const &capture, int width, int height) {
	std::ofstream scap_ofs(filename);
	std::string out_buffer;
	out_buffer += capture.to_string();

	scap_ofs << "#" << width << "\t" << height << std::endl;
	scap_ofs.write(out_buffer.c_str(), out_buffer.size());
	scap_ofs.close();
}

int main(int argc, char** argv) {
	std::string scap_filename(argv[argc - 1]);
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

	preprocess_cluster(1, width, height, &capture);
	std::string final_output_name = scap_filename;
	final_output_name.erase(final_output_name.length() - 5, 5); // remove .scap
	final_output_name += "_out.scap";
	to_scap(final_output_name, capture, width, height);

	return 0;
}