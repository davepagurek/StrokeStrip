#pragma once

#include <vector>
#include <unordered_set>
#include <map>
#include <glm/glm.hpp>

#include "SketchInfo.h"

void preprocess_cluster(int cut_opt, int width, int height, Capture* inout_capture);

void RDP(std::vector<glm::dvec2> const &points, std::vector<size_t> &indices, double epsilon);
void RDP(Sketch &s, double epsilon);

void smooth_stroke(Sketch &s);

class StrokeCut {
public:
	StrokeCut() { }

	Capture cut_sharp_turns(Capture const &capture) const;
	Capture cut_spirals(Capture const &capture, std::map<size_t, size_t> &cutindex_to_index, std::map<size_t, std::vector<size_t>> &c_record) const;
	std::vector<Sketch> cut_evenly(std::vector<Sketch> const &strokes, double cut_distance) const;

	Capture cut_RDP(Capture const &capture) const;

	Capture cut_RDP_sharp_turns(Capture const &capture) const;

	Capture cut_hooks_simple(Capture const &capture) const;

	void prepare_cornucopia(Capture &capture) const;
	Capture cut_cornucopia(Capture const &capture, std::map<size_t, size_t> &cutindex_to_index, 
		std::map<size_t, std::vector<size_t>> &c_record, bool output_orig = false);

	void update_cut_record(int from, int to);
	std::unordered_set<int> get_cut_record(int ind) const;

	std::vector<std::unordered_set<int>> cut_record;
};