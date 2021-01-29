#pragma once

#include <gurobi_c++.h>
#include <vector>
#include <map>
#include <ostream>

#include "Cluster.h"

class StrokeOrientation {
public:
	StrokeOrientation(bool viz = false);
	void orient_strokes(const Input& input);
	void flip_strokes(Input* input);
	void orientation_debug(std::ostream& os, const Input& input);

private:
	GRBEnv grb;
	bool viz;
	std::map<int, std::vector<int>> orientations;

	struct DebugLine {
		glm::dvec2 from;
		glm::dvec2 to;
		std::string color;
	};
	std::vector<DebugLine> debug_lines;

	std::vector<int> orient_cluster_strokes(const Cluster& cluster);

	struct PairOrientation {
		int orientation;
		double weight;
	};
	PairOrientation orient_stroke_pair(const Cluster::Stroke& a, const Cluster::Stroke& b);

	struct PolicyResult {
		std::vector<std::pair<glm::dvec2, glm::dvec2>> violations;
		std::vector<double> connection_angles;
		std::vector<double> connection_dists;
		double shortest_cut;
	};
	PolicyResult evaluate_policy(
		std::vector<int> overlaps_a,
		std::vector<int> overlaps_b,
		const Cluster::Stroke& a,
		const Cluster::Stroke& b,
		int policy);

	double weight_for_angle(double angle);
	int closest_ortho_idx_on_curve(const Cluster::Stroke& from, size_t from_i, const Cluster::Stroke& to);
	size_t closest_idx_on_curve(const Cluster::Stroke& from, size_t from_i, const Cluster::Stroke& to);
	void fill_overlaps(const Cluster::Stroke& from, const Cluster::Stroke& to, std::vector<int>& overlaps);
};