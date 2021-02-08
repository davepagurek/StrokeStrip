#define _USE_MATH_DEFINES
#include <cmath>

#include "StrokeOrientation.h"
#include "Utils.h"
#include <vector>
#include <limits>
#include <glm/gtx/norm.hpp>
#include <algorithm>
#include <numeric>
#include <future>

#include "SvgUtils.h"

const int MAX_VIOLATIONS = 10;

StrokeOrientation::StrokeOrientation(const Context& context): context(context) {}

void StrokeOrientation::orientation_debug(std::ostream& os, const Input& input) {
	input.orientation_svg(os, [&](std::ostream& os) {
		for (auto& line : debug_lines) {
			line.from *= input.thickness;
			line.to *= input.thickness;
			SVG::line(os, line.from.x, line.from.y, line.to.x, line.to.y, 0.5, line.color);
		}
	});
}

void StrokeOrientation::add_debug_line(DebugLine line) {
	std::lock_guard<std::mutex> lock(debug_lock);
	debug_lines.push_back(line);
}

void StrokeOrientation::orient_strokes(Input* input) {
	orientations = map_clusters<std::vector<int>>(*input, [&](Cluster& c) { return orient_cluster_strokes(c);  });
}

void StrokeOrientation::flip_strokes(Input* input) {
	for (auto& kv : input->clusters) {
		auto& cluster = kv.second;
		auto& cluster_orientations = orientations[kv.first];
		for (int stroke = 0; stroke < cluster_orientations.size(); stroke++) {
			if (cluster_orientations[stroke] == -1) {
				std::reverse(cluster.strokes[stroke].points.begin(), cluster.strokes[stroke].points.end());
			}
		}
	}
}

std::vector<int> StrokeOrientation::orient_cluster_strokes(const Cluster& cluster) {
	GRBModel model(context.grb);
	std::vector<GRBVar> vars;
	vars.reserve(cluster.strokes.size());
	for (size_t i = 0; i < cluster.strokes.size(); ++i) {
		vars.push_back(model.addVar(0.0, 1.0, 1.0, GRB_BINARY));
	}

	GRBQuadExpr objective = 0;
	for (size_t i = 0; i < cluster.strokes.size(); ++i) {
		for (size_t j = i + 1; j < cluster.strokes.size(); ++j) {
			auto info = orient_stroke_pair(cluster.strokes[i], cluster.strokes[j]);
			objective += info.weight * info.orientation * (vars[i] - vars[j]) * (vars[i] - vars[j]);

			if (context.debug_viz) {
				std::string color = std::string(info.orientation == 1 ? "rgba(0,255,0," : "rgba(255,0,0,") + std::to_string(std::min(1.0, info.weight * 10)) + std::string(")");
				add_debug_line({ midpoint(cluster.strokes[i].points), midpoint(cluster.strokes[j].points), color });
			}
		}
	}

	// Fix a single term to have a single solution
	model.addConstr(vars[0] == 1);

	model.setObjective(objective, GRB_MINIMIZE);
	model.optimize();

	std::vector<int> cluster_orientations(cluster.strokes.size(), 1);
	for (size_t i = 0; i < cluster.strokes.size(); i++) {
		if (vars[i].get(GRB_DoubleAttr_X) == 0) {
			cluster_orientations[i] = -1;
		}
	}
	return cluster_orientations;
}

StrokeOrientation::PairOrientation StrokeOrientation::orient_stroke_pair(const Cluster::Stroke& a, const Cluster::Stroke& b) {
	std::vector<int> overlaps_a(a.points.size(), -1);
	std::vector<int> overlaps_b(b.points.size(), -1);

	StrokeOrientation::PairOrientation result = { 0, 0.0 };

	// Create orthogonal cross sections
	fill_overlaps(a, b, overlaps_a);
	fill_overlaps(b, a, overlaps_b);

	// Check if this is a direct continuation
	double endpoint_dist = std::numeric_limits<double>::infinity();
	int min_ti = 0; // \in { 0, 1 }
	int min_tj = 0; // \in { 0, 1 }
	for (int ti : { 0, 1 }) {
		for (int tj : { 0, 1 }) {
			double dist = glm::distance(a.points[ti * (a.points.size() - 1)], b.points[tj * (b.points.size() - 1)]);
			if (dist < endpoint_dist) {
				endpoint_dist = dist;
				min_ti = ti;
				min_tj = tj;
			}
		}
	}
	bool endpoint = endpoint_dist < 1e-1;

	bool has_overlap = std::any_of(overlaps_a.begin(), overlaps_a.end(), [](int i) { return i != -1; }) ||
		std::any_of(overlaps_b.begin(), overlaps_b.end(), [](int i) { return i != -1; });

	int policy = 0;
	StrokeOrientation::PolicyResult policy_result = { {}, {}, {}, std::numeric_limits<double>::infinity() };
	if (has_overlap) {
		StrokeOrientation::PolicyResult policy_result_fwd = evaluate_policy(overlaps_a, overlaps_b, a, b, 1);
		StrokeOrientation::PolicyResult policy_result_rev = evaluate_policy(overlaps_a, overlaps_b, a, b, -1);

		if (std::isinf(policy_result_fwd.shortest_cut)) {
			policy = 1;
		}
		else if (std::isinf(policy_result_rev.shortest_cut)) {
			policy = -1;
		}
		else if (policy_result_fwd.violations.size() < 2 != policy_result_rev.violations.size() < 2) {
			policy = policy_result_fwd.violations.size() < 2 ? 1 : -1;
		}
		else if (policy_result_fwd.violations.size() < MAX_VIOLATIONS && policy_result_rev.violations.size() < MAX_VIOLATIONS) {
			policy = policy_result_fwd.shortest_cut > policy_result_rev.shortest_cut ? 1 : -1;
		}
		else {
			policy = policy_result_fwd.violations.size() < policy_result_rev.violations.size() ? 1 : -1;
		}

		policy_result = policy == 1 ? policy_result_fwd : policy_result_rev;

		if (context.debug_viz) {
			for (auto& pair : policy_result.violations) {
				add_debug_line({ pair.first, pair.second, "rgba(255,0,255,0.25)" });
			}
		}
	}

	bool tiny_stroke = false;
	bool continuation = false;

	if (!has_overlap || policy_result.violations.size() > MAX_VIOLATIONS || endpoint) {
		double len_a = total_length(a.points);
		double len_b = total_length(b.points);
		tiny_stroke = std::min(len_a, len_b) < 2.5;
		if (tiny_stroke) {
			policy_result.violations = {};
			if (len_a < len_b) {
				size_t closest_b = closest_idx_on_curve(a, 0, b);
				double dot = glm::dot(tangent(a.points, 0), tangent(b.points, closest_b));
				result.orientation = dot > 0 ? 1 : -1;
				result.weight = 1e-1;
				policy_result.connection_angles = { std::acos(result.orientation * dot) };
				policy_result.connection_dists = { glm::distance(a.points[0], b.points[closest_b]) };
			}
			else {
				size_t closest_a = closest_idx_on_curve(b, 0, a);
				double dot = glm::dot(tangent(b.points, 0), tangent(a.points, closest_a));
				result.orientation = dot > 0 ? 1 : -1;
				result.weight = 1e-1;
				policy_result.connection_angles = { std::acos(result.orientation * dot) };
				policy_result.connection_dists = { glm::distance(b.points[0], a.points[closest_a]) };
			}
		}
		else {
			result.orientation = min_ti != min_tj ? 1 : -1;
			double endpoint_angle = std::acos(result.orientation * glm::dot(
				tangent(a.points, min_ti * (a.points.size() - 1)),
				tangent(b.points, min_tj * (b.points.size() - 1))
			));

			if (!has_overlap || policy_result.violations.size() > MAX_VIOLATIONS || policy_result.connection_angles.empty()) {
				continuation = true;
				policy_result.violations = {};
				policy_result.connection_angles = { endpoint_angle };
				policy_result.connection_dists = { endpoint_dist };
				result.weight = has_overlap ? 1e-2 : 1;
			}
			else {
				bool orig_angle_ok = endpoint_angle / M_PI * 180 < 180 - 40;
				bool orig_angle_good = endpoint_angle / M_PI * 180 <= 90;

				if (orig_angle_ok) {
					continuation = true;
					policy_result.violations = {};
					policy_result.connection_angles = { endpoint_angle };
					policy_result.connection_dists = { endpoint_dist };
					result.weight = orig_angle_good ? 0.8 : 0.1;
				}
			}
		}
	}

	if (!continuation) {
		if (!tiny_stroke) {
			result.orientation = policy;
		}
		if (policy_result.connection_angles.empty() ||
			(policy_result.connection_dists.size() <= 2 &&
				*std::min_element(policy_result.connection_dists.begin(), policy_result.connection_dists.end()) > 2)) {
			result.weight = 2e-2;
		}
		else {
			double angle = std::accumulate(
				policy_result.connection_angles.begin(),
				policy_result.connection_angles.end(),
				0) / policy_result.connection_angles.size();
			result.weight = weight_for_angle(angle) + 1e-1;
			double len_ratio = double(policy_result.connection_angles.size()) / (0.95 * std::min(a.points.size(), b.points.size()));
			len_ratio = std::min(1., std::max(0.5, len_ratio));
			result.weight *= len_ratio;
		}
	}

	double min_dist = 0;
	if (!policy_result.connection_dists.empty()) {
		size_t num_dists = policy_result.connection_dists.size();
		int off = 0.025 * num_dists;
		std::nth_element(policy_result.connection_dists.begin(), policy_result.connection_dists.begin() + off, policy_result.connection_dists.end());
		double min_dist = 3 * policy_result.connection_dists[off];
		result.weight /= 1. + min_dist * min_dist;
		//std::cout << "Dist: " << policy_result.connection_dists[off] << "; weight: " << result.weight << std::endl;
	}

	return result;
}

double StrokeOrientation::weight_for_angle(double angle) {
	double deg = angle / M_PI * 180;
	double sigma = 10. / 3.;
	double low = 10.;
	double high = 20.;
	if (deg < low) {
		return 1;
	}
	else if (deg < high) {
		return std::exp(-(deg - low)*(deg - low) / (2 * sigma * sigma));
	}
	else {
		return 0;
	}
}

StrokeOrientation::PolicyResult StrokeOrientation::evaluate_policy(
	std::vector<int> overlaps_a,
	std::vector<int> overlaps_b,
	const Cluster::Stroke& a,
	const Cluster::Stroke& b,
	int policy
) {
	StrokeOrientation::PolicyResult result{ {}, {}, {}, std::numeric_limits<double>::infinity() };

	// Remove overlaps that don't work for this orientation
	double shortest_cut = std::numeric_limits<double>::infinity();
	for (size_t i = 0; i < overlaps_a.size(); ++i) {
		if (overlaps_a[i] == -1) continue;
		if (policy * glm::dot(tangent(a.points, i), tangent(b.points, overlaps_a[i])) <= 0) {
			result.shortest_cut = std::min(result.shortest_cut, glm::distance(a.points[i], b.points[overlaps_a[i]]));
			overlaps_a[i] = -1;
		}
	}
	for (size_t j = 0; j < overlaps_b.size(); ++j) {
		if (overlaps_b[j] == -1) continue;
		if (policy * glm::dot(tangent(a.points, overlaps_b[j]), tangent(b.points, j)) <= 0) {
			result.shortest_cut = std::min(result.shortest_cut, glm::distance(a.points[overlaps_b[j]], b.points[j]));
			overlaps_b[j] = -1;
		}
	}

	// Check isoline violations
	struct Neighbour {
		int i;
		int j;
		int idx_dist;
		double xsec_len_sq;
	};
	auto check_violations = [&](const Cluster::Stroke& from, const Cluster::Stroke& to, const std::vector<int>& overlaps, const std::vector<int>& other_overlaps) {
		size_t n = overlaps.size();
		std::vector<std::pair<glm::dvec2, glm::dvec2>> curr_violations;
		for (size_t i = 0; i < n; ++i) {
			if (overlaps[i] != -1) {
				// If we have an orthogonal connection, add its angle to the list
				result.connection_angles.push_back(std::acos(policy * glm::dot(tangent(from.points, i), tangent(to.points, overlaps[i]))));
				result.connection_dists.push_back(glm::distance(from.points[i], to.points[overlaps[i]]));
				if (curr_violations.size() > result.violations.size()) result.violations = curr_violations;
				curr_violations.clear();
			}
			else {
				// See if we're connected from the other side
				auto it = std::find(other_overlaps.begin(), other_overlaps.end(), i);
				if (it != other_overlaps.end()) {
					size_t j = it - other_overlaps.begin();
					result.connection_angles.push_back(std::acos(policy * glm::dot(tangent(from.points, i), tangent(to.points, j))));
					result.connection_dists.push_back(glm::distance(from.points[i], to.points[j]));
					if (curr_violations.size() > result.violations.size()) result.violations = curr_violations;
					curr_violations.clear();
					continue;
				}

				// Otherwise, find the isoline this point is on

				std::vector<Neighbour> neighbours;
				// Neighbour to our left
				for (int other_i = int(i) - 1; other_i >= 0; --other_i) {
					if (overlaps[other_i] != -1) {
						neighbours.push_back({
							other_i,
							overlaps[other_i],
							std::abs(int(i) - other_i),
							glm::distance2(from.points[other_i], to.points[overlaps[other_i]])
						});
						break;
					}
				}

				// Neighbour to our right
				for (int other_i = i + 1; other_i < overlaps.size(); ++other_i) {
					if (overlaps[other_i] != -1) {
						neighbours.push_back({
							other_i,
							overlaps[other_i],
							std::abs(int(i) - other_i),
							glm::distance2(from.points[other_i], to.points[overlaps[other_i]])
						});
						break;
					}
				}

				// Neighbours on the other side
				int closest_j_left = -1;
				int closest_j_right = -1;
				for (int j = 0; j < other_overlaps.size(); ++j) {
					if (other_overlaps[j] == -1) continue;
					if (other_overlaps[j] < i && (closest_j_left == -1 || other_overlaps[j] > other_overlaps[closest_j_left])) {
						closest_j_left = j;
					}
					if (other_overlaps[j] > i && (closest_j_right == -1 || other_overlaps[j] < other_overlaps[closest_j_right])) {
						closest_j_right = j;
					}
				}
				if (closest_j_left != -1) {
					neighbours.push_back({
						other_overlaps[closest_j_left],
						closest_j_left,
						std::abs(int(i) - other_overlaps[closest_j_left]),
						glm::distance2(from.points[other_overlaps[closest_j_left]], to.points[closest_j_left])
					});
				}
				if (closest_j_right != -1) {
					neighbours.push_back({
						other_overlaps[closest_j_right],
						closest_j_right,
						std::abs(int(i) - other_overlaps[closest_j_right]),
						glm::distance2(from.points[other_overlaps[closest_j_right]], to.points[closest_j_right])
					});
				}

				if (neighbours.empty()) continue;
				auto neighbour = *std::min_element(
					neighbours.begin(),
					neighbours.end(),
					[](const Neighbour& a, const Neighbour& b) {
						if (a.idx_dist != b.idx_dist) return a.idx_dist < b.idx_dist;
						return a.xsec_len_sq < b.xsec_len_sq;
					});
				int dist = i - neighbour.i;
				int j = neighbour.j + policy * dist;

				// No isoline because we ran out of other curve
				if (j < 0 || j >= to.points.size()) continue;

				// No isoline because it's on its own xsec
				if (other_overlaps[j] != -1) continue;

				double dot = policy * glm::dot(tangent(from.points, i), tangent(to.points, j));
				if (dot < 0) {
					curr_violations.push_back({ from.points[i], to.points[j] });
				}
				else {
					result.connection_angles.push_back(std::acos(dot));
					result.connection_dists.push_back(glm::distance(from.points[i], to.points[j]));
				}
			}
		}
		if (curr_violations.size() > result.violations.size()) result.violations = curr_violations;
		curr_violations.clear();
	};
	check_violations(a, b, overlaps_a, overlaps_b);
	check_violations(b, a, overlaps_b, overlaps_a);

	return result;
}

size_t StrokeOrientation::closest_idx_on_curve(const Cluster::Stroke& from, size_t from_i, const Cluster::Stroke& to) {
	double min_dist = std::numeric_limits<double>::infinity();
	size_t min_idx = 0;
	for (size_t j = 0; j < to.points.size(); ++j) {
		double dist = glm::distance2(from.points[from_i], to.points[j]);
		if (dist < min_dist) {
			min_dist = dist;
			min_idx = j;
		}
	}
	return min_idx;
}

int StrokeOrientation::closest_ortho_idx_on_curve(const Cluster::Stroke& from, size_t from_i, const Cluster::Stroke& to) {
	auto origin = from.points[from_i];
	auto ortho = normal(tangent(from.points, from_i));
	auto p1 = origin - 100. * ortho;
	auto p2 = origin + 100. * ortho;
	auto ints = intersections(to.points, p1, p2);

	if (ints.empty()) {
		return -1;
	}

	double min_dist = std::numeric_limits<float>::infinity();
	int min_idx = -1;
	for (int i = 0; i < ints.size(); ++i) {
		auto& intersection = ints[i];
		double dist = glm::distance2(origin, intersection.pt);
		if (dist < min_dist) {
			min_dist = dist;
			min_idx = i;
		}
	}

	return int(std::round(ints[min_idx].i));
};

void StrokeOrientation::fill_overlaps(const Cluster::Stroke& from, const Cluster::Stroke& to, std::vector<int>& overlaps) {
	for (size_t i = 0; i < from.points.size(); ++i) {
		int j = closest_ortho_idx_on_curve(from, i, to);
		if (j == -1) continue;

		auto vec1 = to.points[j] - from.points[i];
		double dist1 = glm::length(vec1);
		vec1 = glm::normalize(vec1);

		bool vec1_ok = true;
		if (dist1 > 5e-1) {
			auto line1 = from.points[i] + 1e-1*vec1;
			auto line2 = to.points[j] - 1e-1*vec1;
			if (!intersections(from.points, line1, line2).empty()) {
				vec1_ok = false;
			}
		}

		int k = closest_ortho_idx_on_curve(to, j, from);
		bool vec2_ok = false;
		double dist2 = std::numeric_limits<double>::infinity();
		if (k != -1) {
			vec2_ok = true;
			auto vec2 = to.points[j] - from.points[k];
			dist2 = glm::length(vec2);
			vec2 = glm::normalize(vec2);

			if (dist2 > 5e-1) {
				auto line1 = from.points[k] + 1e-1*vec1;
				auto line2 = to.points[j] - 1e-1*vec1;
				if (!intersections(to.points, line1, line2).empty()) {
					vec2_ok = false;
				}
			}
		}

		if (!vec1_ok && !vec2_ok) {
			continue;
		}

		double dist = dist1;
		int used_i = i;
		if (!vec1_ok || (vec2_ok && dist2 < dist1)) {
			used_i = k;
			dist = dist2;
		}

		if (overlaps[used_i] == -1 || dist < glm::distance(from.points[used_i], to.points[overlaps[used_i]])) {
			overlaps[used_i] = j;
		}
	}

	if (false && context.debug_viz) {
		for (size_t i = 0; i < overlaps.size(); i++) {
			if (overlaps[i] != -1) {
				add_debug_line({ from.points[i], to.points[overlaps[i]], "rgba(100,100,255,0.5)" });
			}
		}
	}
}