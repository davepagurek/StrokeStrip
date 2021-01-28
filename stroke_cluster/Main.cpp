#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>
#include <map>

#include <Windows.h>
#include <algorithm>

#include "StrokeCutting.h"

// === Globals ===

double epsilon_small;
double input_median_length;
double min_reparam_step_length = 5.;
size_t min_num_sample_per_stroke = 30;
size_t max_num_sample_per_stroke = 100;

/*double min_reparam_step_length = 5.;

size_t min_num_sample_per_stroke = 30;
size_t max_num_sample_per_stroke = 100;

std::string final_output_name = "";
int width = 0, height = 0;
std::string output_extension = ".scap";
bool to_step_file = false;
bool is_conservative = false;

int prox_type;

extern double post_parallel_threshold;
extern double far_stroke_distance_threshold;



double input_median_length;

unsigned long start_time;
unsigned long no_pre_start_time;

std::map<size_t, size_t> reindex_to_index;
std::map<size_t, size_t> index_to_reindex;
std::map<size_t, size_t> cutindex_to_index;
std::map<size_t, std::vector<size_t>> c_record;

using namespace std;*/

void preprocess_cluster(int cut_opt, int width, int height, Capture* inout_capture) {
	Capture& capture = *inout_capture;

	StrokeCut cut;
	std::map<size_t, size_t> cutindex_to_index;
	std::map<size_t, std::vector<size_t>> c_record;
	std::map<size_t, size_t> reindex_to_index;
	std::map<size_t, size_t> index_to_reindex;

	// cut
	//if (!precomputed) {
		capture = cut.cut_hooks_simple(capture);
		capture = cut.cut_RDP_sharp_turns(capture);

		cut.prepare_cornucopia(capture);
		capture = cut.cut_cornucopia(capture, cutindex_to_index, c_record, cut_opt != 0);
		capture = cut.cut_hooks_simple(capture);

		std::map<size_t, std::vector<size_t>> c_record2;
		std::map<size_t, size_t> cutindex_to_index2;
		//capture = cut.cut_spirals(capture, cutindex_to_index2, c_record2);

		for (auto const ind : cutindex_to_index2) {
			if (cutindex_to_index.count(ind.first) == 0)
				cutindex_to_index[ind.first] = ind.second;
		}

		// c_record2 => c_record
		{
			for (auto &record : c_record) {
				std::vector<size_t> new_cuts;
				for (auto &cut : record.second) {
					if (c_record2.count(cut) > 0) {
						for (auto &c2 : c_record2[cut]) {
							new_cuts.emplace_back(c2);
						}
					}
					else
						new_cuts.emplace_back(cut);
				}
				record.second = new_cuts;
			}
		}

		// reparameterize
		for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
			//capture.getSketchedPolyline(i).reparameterize(min_reparam_step_length);
			double step_length = capture.getSketchedPolyline(i).totalLen() / min_num_sample_per_stroke;
			step_length = std::min(step_length, min_reparam_step_length);
			capture.getSketchedPolyline(i).reparameterize(step_length);
		}
	//}

	// filter size
	// Debug
	//if (!precomputed) {
		Capture filtered_capture;
		//double long_len = 15;
		double long_len = std::min(5., 0.015 * std::sqrt(width * width + height * height));
		for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
			if (capture.getSketchedPolyline(i).totalLen() >= long_len) {
				filtered_capture.sketchedPolylines.emplace_back(capture.getSketchedPolyline(i));
			}
		}

		capture = filtered_capture;
	//}

	// compute median length
	std::vector<double> lengths;
	lengths.reserve(capture.sketchedPolylines.size());

	for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
		lengths.emplace_back(capture.getSketchedPolyline(i).totalLen());
	}

	std::nth_element(lengths.begin(), lengths.begin() + lengths.size() / 2, lengths.end());
	input_median_length = lengths[lengths.size() / 2];

	// Reorder stroke indices
	{
		size_t re_ind = 0;
		for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
			reindex_to_index.emplace(re_ind, capture.getSketchedPolyline(i).stroke_ind);
			index_to_reindex.emplace(capture.getSketchedPolyline(i).stroke_ind, re_ind);
			cut.update_cut_record(-1 * (int)capture.getSketchedPolyline(i).stroke_ind, re_ind);
			capture.getSketchedPolyline(i).stroke_ind = re_ind++;
		}
	}
}

/*void to_svg(std::string const &filename, Capture const &capture) {
	std::ofstream scap_ofs(filename);
	std::string svg_string = capture.to_svg_string(width, height);
	scap_ofs.write(svg_string.c_str(), svg_string.size());
	scap_ofs.close();
}*/

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

	// Assign thickness
	epsilon_small = capture.thickness;

	//std::cout << "epsilon_small: " << epsilon_small << std::endl;

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