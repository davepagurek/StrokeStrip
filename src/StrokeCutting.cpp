#define _USE_MATH_DEFINES
#include <cmath>

#include "StrokeCutting.h"

#include "glm/glm.hpp"

#include "Cornucopia/Cornucopia.h"
#include "Cornucopia/PathFinder.h"
#include "Cornucopia/GraphConstructor.h"
#include "Cornucopia/Resampler.h"
#include "Cornucopia/PrimitiveFitter.h"
#include "Cornucopia/Preprocessing.h"
#include "Cornucopia/CurvePrimitive.h"

// TODO: If I feel fancy. I may want to combine these two functions into one cut function
// which takes a higher order function to do the cut.

double short_seg_len = 15;//2;
double short_ratio = 0.15;
double max_curv_threshold = 1./2;

double curvature_cut_min_length = 5.;
double max_curvature_ratio = 3;

double min_reparam_step_length = 5.;
size_t min_num_sample_per_stroke = 30;
size_t max_num_sample_per_stroke = 100;
double input_median_length;

double epsilon_small;

void preprocess_cluster(int cut_opt, int width, int height, Capture* inout_capture) {
	Capture& capture = *inout_capture;

	epsilon_small = capture.thickness;

	StrokeCut cut;
	std::map<size_t, size_t> cutindex_to_index;
	std::map<size_t, std::vector<size_t>> c_record;
	std::map<size_t, size_t> reindex_to_index;
	std::map<size_t, size_t> index_to_reindex;

	// cut
	//if (!precomputed) {
	//capture = cut.cut_hooks_simple(capture);
	//capture = cut.cut_RDP_sharp_turns(capture);

	cut.prepare_cornucopia(capture);
	capture = cut.cut_cornucopia(capture, cutindex_to_index, c_record, cut_opt != 0);
	//capture = cut.cut_hooks_simple(capture);

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

	capture.thickness = epsilon_small;
}

glm::dvec2 to_glm_vec(const SketchUI::Point2D& pt) {
	return glm::dvec2(pt.x, pt.y);
}

Capture StrokeCut::cut_sharp_turns(Capture const &capture) const {
	double sharp_angle = M_PI / 2;
	double sharp_threshold = glm::length(glm::dvec2(std::cos(sharp_angle), std::sin(sharp_angle)) - glm::dvec2(1, 0));

	int max_id = 0;
	for (auto const &s : capture.sketchedPolylines) {
		max_id = std::max(max_id, s.stroke_ind);
	}

	size_t stroke_ind = capture.getSketchedPolylineSize();
	auto cut_stroke = [&stroke_ind, &max_id, &sharp_angle](Sketch const &s)->std::vector<Sketch> {
		std::vector<double> local_turning;
		local_turning.emplace_back(0);
		for (size_t i = 1; i + 1 < s.points.size(); i++) {
			size_t prev = i - 1;
			size_t next = i + 1;

			glm::dvec2 prev_p(s.points[prev].first.x, s.points[prev].first.y);
			glm::dvec2 cur_p(s.points[i].first.x, s.points[i].first.y);
			glm::dvec2 next_p(s.points[next].first.x, s.points[next].first.y);

			glm::dvec2 v1 = cur_p - prev_p;
			glm::dvec2 v2 = next_p - cur_p;

			if (glm::length(v1) > 0)
				v1 = glm::normalize(v1);
			if (glm::length(v2) > 0)
				v2 = glm::normalize(v2);

			glm::dvec2 norm_v1(-v1.y, v1.x);

			// Counter-clockwise: +
			double sign = (glm::dot(norm_v1, v2) > 0) ? 1 : -1;

			local_turning.emplace_back(sign * std::acos(glm::dot(v1, v2)));
		}
		local_turning.emplace_back(0);

		// Sum turning angles locally
		std::vector<double> large_sum_turning;
		large_sum_turning.resize(local_turning.size(), 0);
		std::vector<double> small_sum_turning;
		small_sum_turning.resize(local_turning.size(), 0);
		int half_width = 5;
		int small_half_width = 0;
		{
			for (size_t i = 0; i < local_turning.size(); i++) {
				for (int j = i - half_width; j <= i + half_width; j++) {
					double angle = 0;
					if (j >= 0 && j < local_turning.size())
						angle = local_turning[j];
					large_sum_turning[i] += angle;
				}
				for (int j = i - small_half_width; j <= i + small_half_width; j++) {
					double angle = 0;
					if (j >= 0 && j < local_turning.size())
						angle = local_turning[j];
					small_sum_turning[i] += angle;
				}
			}
		}

		// Check to see if there's any zero crossing for curvature changing rate
		std::vector<bool> cut_pos;
		cut_pos.resize(s.points.size(), false);
		{
			for (size_t i = 1; i + 1 < large_sum_turning.size(); i++) {
				double change_from = std::abs(small_sum_turning[i]) - std::abs(small_sum_turning[i - 1]);
				double change_to = std::abs(small_sum_turning[i + 1]) - std::abs(small_sum_turning[i]);
				if (change_from * change_to < 0 &&
					std::abs(small_sum_turning[i]) >= sharp_angle) {
					double angle = small_sum_turning[i] / M_PI * 180;
					cut_pos[i] = true;
					std::cout << s.stroke_ind << ": " << i << " - " << angle << std::endl;
				}
			}
		}

		std::vector<Sketch> cut;
		cut_pos[s.points.size() - 1] = true;
		Sketch cur_s;
		cur_s.stroke_ind = s.stroke_ind;

		for (size_t i = 0; i < s.points.size(); i++) {
			cur_s.points.emplace_back(s.points[i]);

			if (cut_pos[i]) {
				cur_s.group_ind = s.group_ind;
				if (cut.size() == 0)
					cur_s.stroke_ind = s.stroke_ind;
				else
					cur_s.stroke_ind = stroke_ind++;

				cut.emplace_back(cur_s);
				cur_s.points.clear();
				cur_s.stroke_ind = ++max_id;
				cur_s.points.emplace_back(s.points[i]);
			}
		}

		return cut;
	};
	Capture cut_capture;

	for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
		std::vector<Sketch> cut_strokes = cut_stroke(capture.sketchedPolylines[i]);
		for (auto const &s : cut_strokes) {
			cut_capture.sketchedPolylines.emplace_back(s);
		}
	}

	return cut_capture;
}

Capture StrokeCut::cut_spirals(Capture const &capture, std::map<size_t, size_t> &cutindex_to_index, std::map<size_t, std::vector<size_t>> &c_record) const {
	int max_id = 0;
	for (auto const &s : capture.sketchedPolylines) {
		max_id = std::max(max_id, s.stroke_ind);
	}
	auto cut_stroke = [&max_id](Sketch const &s)->std::vector<Sketch> {
		// I'm lazy, so I'll draw a line perpendicular to the starting segment
		// and use that to cut the stroke.

		std::vector<bool> cut_pos;
		std::vector<double> dist;
		double prev_dist = 0.;

		cut_pos.resize(s.points.size(), false);

		float turning_number = 0.f;
		size_t prev_ind = 0;
		for (size_t i = 0; i + 2 < s.points.size(); i++) {
			size_t j = i + 1;
			size_t k = i + 2;

			glm::dvec2 cur_p(s.points[i].first.x, s.points[i].first.y);
			glm::dvec2 next_p(s.points[j].first.x, s.points[j].first.y);
			glm::dvec2 next2_p(s.points[k].first.x, s.points[k].first.y);

			glm::dvec2 v1 = next_p - cur_p;
			glm::dvec2 v2 = next2_p - next_p;
			
			float angle2 = std::atan2(v2.y, v2.x);
			float angle1 = std::atan2(v1.y, v1.x);
			//angle1 = (angle1 < 0) ? angle1 + 2 * M_PI : angle1;
			//angle2 = (angle2 < 0) ? angle2 + 2 * M_PI : angle2;
			float angle = angle2 - angle1;
			angle = (angle > M_PI) ? 2 * M_PI - angle : angle;
			angle = (angle < -M_PI) ? 2 * M_PI + angle : angle;

			turning_number += angle / (2 * M_PI);

			// Use a threshold near 1
			// If it's nearly a round, cut it into 4 pieces
			if (std::fabsf(turning_number) > 0.8) {
				size_t cut_step = (prev_ind + j) / 4;
				for (size_t cut = 1; cut < 4 && prev_ind + cut * cut_step < j; cut++) {
					cut_pos[prev_ind + cut * cut_step] = true;
				}
				cut_pos[j] = true;
				prev_ind = j;
				turning_number = 0.f;
			}
		}

		std::vector<Sketch> cut;
		cut_pos[s.points.size() - 1] = true;
		Sketch cur_s;
		cur_s.stroke_ind = s.stroke_ind;

		for (size_t i = 0; i < s.points.size(); i++) {
			cur_s.points.emplace_back(s.points[i]);

			if (cut_pos[i]) {
				cur_s.group_ind = s.group_ind;
				cut.emplace_back(cur_s);
				cur_s.points.clear();
				cur_s.stroke_ind = ++max_id;
				cur_s.points.emplace_back(s.points[i]);
			}
		}

		return cut;
	};
	Capture cut_capture;

	for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
		std::vector<Sketch> cut_strokes = cut_stroke(capture.sketchedPolylines[i]);
		c_record[capture.sketchedPolylines[i].stroke_ind] = std::vector<size_t>();
		for (auto const &s : cut_strokes) {
			c_record[capture.sketchedPolylines[i].stroke_ind].emplace_back(s.stroke_ind);
			cut_capture.sketchedPolylines.emplace_back(s);

			if (cut_strokes.size() > 1)
				cutindex_to_index[s.stroke_ind] = capture.sketchedPolylines[i].stroke_ind;
		}
	}

	return cut_capture;
}

std::vector<Sketch> StrokeCut::cut_evenly(std::vector<Sketch> const &strokes, double cut_distance) const {
	int max_id = 0;
	for (auto const &s : strokes) {
		max_id = std::max(max_id, s.stroke_ind);
	}

	auto cut_stroke = [&](Sketch const &s)->std::vector<Sketch> {
		// I'm lazy, so I'll draw a line perpendicular to the starting segment
		// and use that to cut the stroke.

		std::vector<bool> cut_pos;
		std::vector<double> dist;
		cut_pos.resize(s.points.size(), false);
		s.computeSegLen(dist);

		double sum_dist = 0.;
		for (size_t i = 0; i < dist.size(); i++) {
			sum_dist += dist[i];
			if (sum_dist > cut_distance) {
				cut_pos[i + 1] = true;
				sum_dist -= cut_distance;
			}
		}

		std::vector<Sketch> cut;
		cut_pos[s.points.size() - 1] = true;
		Sketch cur_s;
		cur_s.stroke_ind = s.stroke_ind;

		for (size_t i = 0; i < s.points.size(); i++) {
			cur_s.points.emplace_back(s.points[i]);

			if (cut_pos[i]) {
				cur_s.group_ind = s.group_ind;
				cut.emplace_back(cur_s);
				cur_s.points.clear();
				cur_s.stroke_ind = ++max_id;
				cur_s.points.emplace_back(s.points[i]);
			}
		}

		return cut;
	};
	std::vector<Sketch> results;

	for (size_t i = 0; i < strokes.size(); i++) {
		std::vector<Sketch> cut_strokes = cut_stroke(strokes[i]);
		results.insert(results.end(), cut_strokes.begin(), cut_strokes.end());
	}

	return results;
}

void RDP(std::vector<glm::dvec2> const &points, std::vector<size_t> &indices, double epsilon) {
	// Find the point furtherest away from the line connecting the first and last points
	size_t far_index = 0;
	size_t far_x_index = 0;
	bool folding = false;
	double max_distance = -1;
	double max_x_distance = -1;
	{
		glm::dvec2 simp_vector = points[indices.back()] - points[indices.front()];
		glm::dvec2 normed = (glm::length(simp_vector) > std::numeric_limits<double>::epsilon()) ? glm::normalize(simp_vector) : glm::dvec2(1, 0);

		glm::dvec2 tan(0.0, 0.0);
		glm::dvec2 pre_tan(0, 0);

		for (size_t i = 0; i < indices.size(); i++) {
			glm::dvec2 p = points[indices[i]];

			glm::dvec2 dist_vector = p - points[indices.front()];
			double dist = glm::length(dist_vector - glm::dot(dist_vector, normed) * normed);
			if (dist > max_distance) {
				max_distance = dist;
				far_index = i;
			}

			// Extension: zig-zag stroke
			if (i + 1 < indices.size()) {
				tan = points[indices[i + 1]] - points[indices[i]];
			}
			if (i >= 1) {
				pre_tan = points[indices[i]] - points[indices[i - 1]];
			}

			// A sharp turn
			dist = std::abs(glm::dot(dist_vector, normed));
			glm::dvec2 dist_vector2 = p - points[indices.back()];
			dist = std::min(dist, std::abs(glm::dot(dist_vector2, normed)));
			if ((glm::dot(tan, simp_vector) < 0 ||
				glm::dot(pre_tan, simp_vector) < 0) &&
				dist > max_x_distance) {
				max_x_distance = dist;
				far_x_index = i;
			}
		}
	}

	if (max_x_distance > max_distance) {
		max_distance = max_x_distance;
		far_index = far_x_index;
		folding = true;
	}

	// Call simplification recursively
	if (max_distance > epsilon || folding && max_distance > 5) {
		std::vector<size_t> sub_left_indices, sub_right_indices;
		sub_left_indices.resize(far_index + 1);
		sub_right_indices.resize(indices.size() - far_index);

		std::copy(indices.begin(), indices.begin() + far_index + 1, sub_left_indices.begin());
		std::copy(indices.begin() + far_index, indices.end(), sub_right_indices.begin());

		RDP(points, sub_left_indices, epsilon);
		RDP(points, sub_right_indices, epsilon);

		indices.clear();
		indices.resize(sub_left_indices.size() + sub_right_indices.size() - 1);
		std::copy(sub_left_indices.begin(), sub_left_indices.end(), indices.begin());
		std::copy(sub_right_indices.begin() + 1, sub_right_indices.end(), indices.begin() + sub_left_indices.size());
	}
	else { // Or terminate
		std::vector<size_t> remaining;
		remaining.emplace_back(size_t(indices.front()));
		remaining.emplace_back(size_t(indices.back()));

		indices = remaining;
	}
}

void RDP(Sketch &s, double epsilon) {
	std::vector<glm::dvec2> points;
	std::vector<size_t> indices;
	for (size_t i = 0; i < s.points.size(); i++) {
		points.emplace_back(to_glm_vec(s.points[i].first));
		indices.emplace_back(i);
	}
	RDP(points, indices, epsilon);

	s.points.clear();

	for (size_t i : indices) {
		s.points.emplace_back(SketchPoint(points[i].x, points[i].y), 0);
	}
}

// https://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm
Capture StrokeCut::cut_RDP(Capture const &capture) const {
	double RDP_epsilon = 100;//35.;

	int max_id = 0;
	for (auto const &s : capture.sketchedPolylines) {
		max_id = std::max(max_id, s.stroke_ind);
	}

	auto cut_stroke = [&max_id, &RDP_epsilon](Sketch const &s)->std::vector<Sketch> {
		// Simplify using RDP
		std::vector<glm::dvec2> points;
		std::vector<size_t> indices;
		for (size_t i = 0; i < s.points.size(); i++) {
			points.emplace_back(to_glm_vec(s.points[i].first));
			indices.emplace_back(i);
		}
		RDP(points, indices, RDP_epsilon);

		// Check to see if there's any zero crossing for curvature changing rate
		std::vector<bool> cut_pos;
		cut_pos.resize(s.points.size(), false);
		{
			for (size_t i : indices) {
				cut_pos[i] = true;
			}
			cut_pos.front() = false;
			cut_pos.back() = false;
		}

		std::vector<Sketch> cut;
		cut_pos[s.points.size() - 1] = true;
		Sketch cur_s;
		cur_s.stroke_ind = s.stroke_ind;

		for (size_t i = 0; i < s.points.size(); i++) {
			cur_s.points.emplace_back(s.points[i]);

			if (cut_pos[i]) {
				cur_s.group_ind = s.group_ind;
				if (cut.size() == 0)
					cur_s.stroke_ind = s.stroke_ind;
				else
					cur_s.stroke_ind = ++max_id;

				cut.emplace_back(cur_s);
				cur_s.points.clear();
				cur_s.points.emplace_back(s.points[i]);
			}
		}

		return cut;
	};
	Capture cut_capture;

	for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
		std::vector<Sketch> cut_strokes = cut_stroke(capture.sketchedPolylines[i]);
		for (auto const &s : cut_strokes) {
			cut_capture.sketchedPolylines.emplace_back(s);
		}
	}

	return cut_capture;
}

Capture StrokeCut::cut_RDP_sharp_turns(Capture const &capture) const {
	double RDP_epsilon = 10.;

	auto compute_cut_positions_RDP = [&RDP_epsilon](Sketch const &s, std::vector<bool> &cut_pos) {
		// Simplify using RDP
		std::vector<glm::dvec2> points;
		std::vector<size_t> indices;
		for (size_t i = 0; i < s.points.size(); i++) {
			points.emplace_back(to_glm_vec(s.points[i].first));
			indices.emplace_back(i);
		}
		double cur_epsilon = s.totalLen() / 20;
		cur_epsilon = std::max(cur_epsilon, RDP_epsilon);
		//std::cout << "RDP: " << cur_epsilon << std::endl;
		RDP(points, indices, cur_epsilon);

		cut_pos.resize(s.points.size(), false);
		{
			for (size_t i : indices) {
				cut_pos[i] = true;
			}
			cut_pos.front() = false;
			cut_pos.back() = false;
		}

		cut_pos[s.points.size() - 1] = true;
	};

	double sharp_angle = M_PI / 3;
	double sharp_threshold = glm::length(glm::dvec2(std::cos(sharp_angle), std::sin(sharp_angle)) - glm::dvec2(1, 0));

	auto compute_cut_positions_sharp_turn = [&sharp_angle](Sketch const &s, std::vector<bool> &cut_pos) {
		std::vector<double> local_turning;
		local_turning.emplace_back(0);
		for (size_t i = 1; i + 1 < s.points.size(); i++) {
			size_t prev = i - 1;
			size_t next = i + 1;

			glm::dvec2 prev_p(s.points[prev].first.x, s.points[prev].first.y);
			glm::dvec2 cur_p(s.points[i].first.x, s.points[i].first.y);
			glm::dvec2 next_p(s.points[next].first.x, s.points[next].first.y);

			glm::dvec2 v1 = cur_p - prev_p;
			glm::dvec2 v2 = next_p - cur_p;

			if (glm::length(v1) > 0)
				v1 = glm::normalize(v1);
			if (glm::length(v2) > 0)
				v2 = glm::normalize(v2);

			glm::dvec2 norm_v1(-v1.y, v1.x);

			// Counter-clockwise: +
			double sign = (glm::dot(norm_v1, v2) > 0) ? 1 : -1;

			local_turning.emplace_back(sign * std::acos(glm::dot(v1, v2)));
		}
		local_turning.emplace_back(0);

		// Sum turning angles locally
		std::vector<double> large_sum_turning;
		large_sum_turning.resize(local_turning.size(), 0);
		std::vector<double> small_sum_turning;
		small_sum_turning.resize(local_turning.size(), 0);
		int half_width = 5;
		int small_half_width = 0;
		{
			for (size_t i = 0; i < local_turning.size(); i++) {
				for (int j = i - half_width; j <= i + half_width; j++) {
					double angle = 0;
					if (j >= 0 && j < local_turning.size())
						angle = local_turning[j];
					large_sum_turning[i] += angle;
				}
				for (int j = i - small_half_width; j <= i + small_half_width; j++) {
					double angle = 0;
					if (j >= 0 && j < local_turning.size())
						angle = local_turning[j];
					small_sum_turning[i] += angle;
				}
			}
		}

		// Check to see if there's any zero crossing for curvature changing rate
		cut_pos.resize(s.points.size(), false);
		{
			for (size_t i = 1; i + 1 < large_sum_turning.size(); i++) {
				double change_from = std::abs(small_sum_turning[i]) - std::abs(small_sum_turning[i - 1]);
				double change_to = std::abs(small_sum_turning[i + 1]) - std::abs(small_sum_turning[i]);
				if (change_from * change_to < 0 &&
					std::abs(small_sum_turning[i]) >= sharp_angle) {
					double angle = small_sum_turning[i] / M_PI * 180;
					cut_pos[i] = true;
					/*cut_pos[i - 1] = true;
					cut_pos[i + 1] = true;*/
#ifdef DEBUG_PRINTOUT
					std::cout << s.stroke_ind << ": " << i << " - " << angle << std::endl;
#endif
				}
			}
		}

		cut_pos[s.points.size() - 1] = true;
	};

	int max_id = 0;
	for (auto const &s : capture.sketchedPolylines) {
		max_id = std::max(max_id, s.stroke_ind);
	}

	auto cut_stroke = [&max_id](Sketch const &s, std::vector<bool> const &cut_pos)->std::vector<Sketch> {
		std::vector<Sketch> cut;
		Sketch cur_s;
		cur_s.stroke_ind = s.stroke_ind;

		for (size_t i = 0; i < s.points.size(); i++) {
			cur_s.points.emplace_back(s.points[i]);

			if (cut_pos[i]) {
				cur_s.group_ind = s.group_ind;
				if (cut.size() == 0)
					cur_s.stroke_ind = s.stroke_ind;
				else
					cur_s.stroke_ind = ++max_id;

				cut.emplace_back(cur_s);
				cur_s.points.clear();
				cur_s.points.emplace_back(s.points[i]);
			}
		}

		return cut;
	};
	Capture cut_capture;

	for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
		std::vector<bool> cut_pos_RDP;
		compute_cut_positions_RDP(capture.sketchedPolylines[i], cut_pos_RDP);
		
		std::vector<glm::dvec2> points;
		std::vector<size_t> indices;
		for (size_t j = 0; j < capture.sketchedPolylines[i].points.size(); j++) {
			points.emplace_back(to_glm_vec(capture.sketchedPolylines[i].points[j].first));
			indices.emplace_back(j);
		}
		double epsilon = 3;
		RDP(points, indices, epsilon);

		Sketch sharp_sketch;
		for (size_t ind : indices) {
			sharp_sketch.points.emplace_back(capture.sketchedPolylines[i].points[ind]);
		}

		std::vector<bool> cut_pos_simple;
		compute_cut_positions_sharp_turn(sharp_sketch, cut_pos_simple);

		std::vector<bool> cut_pos_sharp;
		cut_pos_sharp.resize(capture.sketchedPolylines[i].points.size(), false);
		for (size_t j = 0; j < cut_pos_simple.size(); j++) {
			cut_pos_sharp[indices[j]] = cut_pos_simple[j];
		}

		std::vector<bool> cut_pos;
		for (size_t j = 0; j < cut_pos_RDP.size(); j++) {
			cut_pos.push_back(cut_pos_RDP[j] && cut_pos_sharp[j]);
		}

		std::vector<Sketch> cut_strokes = cut_stroke(capture.sketchedPolylines[i], cut_pos);
		for (auto const &s : cut_strokes) {
			cut_capture.sketchedPolylines.emplace_back(s);
		}
	}

	return cut_capture;
}

Capture StrokeCut::cut_hooks_simple(Capture const &capture) const {
	int max_id = 0;
	for (auto const &s : capture.sketchedPolylines) {
		max_id = std::max(max_id, s.stroke_ind);
	}

	auto cut_stroke = [&max_id](Sketch const &s)->std::vector<Sketch> {
		std::vector<bool> cut_pos;
		cut_pos.resize(s.points.size(), false);
		if (s.points.size() > 2) {
			auto is_a_turn = [](SketchPoint p1, SketchPoint p2, SketchPoint p3)->bool {
				return glm::dot((to_glm_vec(p1) - to_glm_vec(p2)), (to_glm_vec(p2) - to_glm_vec(p3))) < 0;
			};
			cut_pos[1] = is_a_turn(s.points[0].first, s.points[1].first, s.points[2].first);
			cut_pos[cut_pos.size() - 2] = is_a_turn(s.points[s.points.size() - 1].first, 
				s.points[s.points.size() - 2].first, s.points[s.points.size() - 3].first);

			for (size_t i = 1; i + 1< s.points.size(); i++) {
				cut_pos[i] = is_a_turn(s.points[i-1].first, s.points[i].first, s.points[i+1].first);
			}
		}

		std::vector<Sketch> cut;
		cut_pos[s.points.size() - 1] = true;
		Sketch cur_s;
		cur_s.stroke_ind = s.stroke_ind;

		for (size_t i = 0; i < s.points.size(); i++) {
			cur_s.points.emplace_back(s.points[i]);

			if (cut_pos[i]) {
				cur_s.group_ind = s.group_ind;
				if (cut.size() == 0)
					cur_s.stroke_ind = s.stroke_ind;
				else
					cur_s.stroke_ind = ++max_id;

				cut.emplace_back(cur_s);
				cur_s.points.clear();
				cur_s.points.emplace_back(s.points[i]);
			}
		}

		return cut;
	};
	Capture cut_capture;

	for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
		std::vector<Sketch> cut_strokes = cut_stroke(capture.sketchedPolylines[i]);

		std::cout << "Cut to " << cut_strokes.size() << std::endl;

		// Keep the largest one
		double max_length = 0;
		Sketch max_s = capture.sketchedPolylines[i];
		for (auto const &s : cut_strokes) {
			if (s.totalLen() > max_length) {
				max_length = s.totalLen();
				max_s = s;
				max_s.stroke_ind = capture.sketchedPolylines[i].stroke_ind;
			}
		}
		cut_capture.sketchedPolylines.emplace_back(max_s);
	}

	return cut_capture;
}

// Assume evenly sampled
void smooth_stroke(Sketch &s) {
	s.reindex();

	// Need to have at least min_num_sample_per_stroke (> 7)
	if (s.points.size() < min_num_sample_per_stroke/2) {
		s.reparameterize(s.totalLen() / (min_num_sample_per_stroke/2));
	}

	assert(s.points.size() >= 7);

	// Handle the two ends by mimicing laplacian
	// http://graphics.stanford.edu/courses/cs468-12-spring/LectureSlides/06_smoothing.pdf
	auto get_weighted_laplacian = [](glm::dvec2 p1, glm::dvec2 p2, glm::dvec2 p3)->glm::dvec2 {
		double len12 = glm::distance(p1, p2);
		double len23 = glm::distance(p2, p3);

		return (p1 / len12 + p3 / len23) / (1/len12 + 1/len23) - p2;
	};

	auto handle_end = [&get_weighted_laplacian](Sketch &s, size_t end_ind, int to_neighbor) {
		int n_ind = end_ind + to_neighbor;
		assert(n_ind >= 0 && n_ind < s.points.size());

		glm::dvec2 n_lap = get_weighted_laplacian(to_glm_vec(s.points[n_ind + to_neighbor].first),
			to_glm_vec(s.points[n_ind].first),
			to_glm_vec(s.points[end_ind].first));
		glm::dvec2 end_lap = get_weighted_laplacian(to_glm_vec(s.points[n_ind].first),
			to_glm_vec(s.points[end_ind].first),
			to_glm_vec(s.points[end_ind - to_neighbor].first));

		double n_lap_len = glm::distance(n_lap, glm::dvec2(0));

		// Normalize the end laplacian
		// If degenerated, use the norm
		glm::dvec2 tan = to_glm_vec(s.points[end_ind - to_neighbor].first) - to_glm_vec(s.points[n_ind].first);
		glm::dvec2 norm_end_lap = (glm::distance(end_lap, glm::dvec2(0)) > std::numeric_limits<double>::epsilon()) ? 
			glm::normalize(end_lap) : glm::normalize(glm::dvec2(-tan.y, tan.x));

		// Pick the direction aligned with the neighbor laplacian
		if (glm::dot(norm_end_lap, n_lap) < 0) {
			norm_end_lap = -norm_end_lap;
		}

		glm::dvec2 end_pos = to_glm_vec(s.points[end_ind].first);
		glm::dvec2 end_new_pos = end_lap + end_pos - n_lap_len * norm_end_lap;

		s.points[end_ind].first.x = end_new_pos.x;
		s.points[end_ind].first.y = end_new_pos.y;

		/*{
		glm::dvec2 n_pos = to_glm_vec(s.points[n_ind].first);
		draw_log.record_to(std::unique_ptr<DrawSegment>(new DrawSegment(n_pos, n_pos + 5. * n_lap,
		nvgRGBA(150, 50, 12, 255))), FilterAttribute(s.group_ind, s.stroke_ind), 0);
		draw_log.record_to(std::unique_ptr<DrawSegment>(new DrawSegment(end_pos, end_pos + 5. * end_lap,
		nvgRGBA(150, 50, 12, 255))), FilterAttribute(s.group_ind, s.stroke_ind), 0);
		}*/
	};
	handle_end(s, 1, 1);
	handle_end(s, s.points.size() - 2, -1);

	// Smooth
	{
		auto smooth_range_2 = [](Sketch &s, size_t start, size_t end) {
			// We don't smooth the endpoints
			assert(start > 1 && start + 2 < s.points.size());
			assert(end > 1 && end + 2 < s.points.size());

			std::vector<SketchPoint> smoothed_points;
			smoothed_points.reserve(end - start + 1);
			for (size_t i = start; i <= end; i++) {
				SketchPoint m = (s.points[i - 1].first + s.points[i + 1].first);
				m /= 2;
				SketchPoint m_prim = (s.points[i - 2].first + s.points[i + 2].first);
				m_prim /= 2;

				SketchPoint offset = m - m_prim;
				offset /= 4;

				smoothed_points.emplace_back(m + offset);
			}

			for (size_t i = start; i <= end; i++) {
				s.points[i].first = smoothed_points[i - start];
			}
		};

		auto smooth_range = smooth_range_2;
		size_t overall_start = 2;
		size_t overall_end = s.points.size() - 3;

		if (s.totalLen() < input_median_length * 0.2) {
			smooth_range(s, overall_start, overall_end);
		}
		else { // Use short_ratio for the hook region
			// Left
			size_t left_end;
			double len = 0;
			double total_len = s.totalLen();
			for (size_t i = 0; i < s.points.size(); i++) {
				if (i > 0) {
					len += (s.points[i].first - s.points[i - 1].first).SqLength();
				}
				double t = len / total_len;

				if (t >= short_ratio) {
					left_end = i;
					break;
				}
			}

			// Right
			left_end = std::max(left_end, overall_start);
			size_t left_start = overall_start;
			size_t right_start = std::max(overall_start, s.points.size() - 1 - left_end);
			size_t right_end = overall_end;

			if (left_end <= right_end) {
				smooth_range(s, left_start, left_end);
				smooth_range(s, right_start, right_end);
			}
		}
	}

	handle_end(s, 1, 1);
	handle_end(s, s.points.size() - 2, -1);
}

Capture StrokeCut::cut_cornucopia(Capture const &capture, std::map<size_t, size_t> &cutindex_to_index, 
								  std::map<size_t, std::vector<size_t>> &c_record, bool output_orig) {
	int max_id = 0;
	for (auto const &s : capture.sketchedPolylines) {
		max_id = std::max(max_id, s.stroke_ind);
	}

	auto is_too_short = [](Sketch const &orig_s, Sketch const &s)->bool {
		double seg_len = s.totalLen();
		double short_threshold = std::min(short_seg_len, orig_s.totalLen() * short_ratio);
		if (seg_len < short_threshold || s.points.size() == 2) {
#ifdef DEBUG_PRINTOUT
			std::cout << "\t Cut: " << seg_len << " from " << orig_s.totalLen() << std::endl;
#endif
		}
		return (seg_len < short_threshold);
		//return s.totalLen() < short_seg_len;
	};

	auto cut_stroke = [&](Sketch const &s, bool &is_closed)->std::vector<Sketch> {
		std::vector<Sketch> cut;

		{
			using namespace Cornu;
			Fitter fitter;
			Parameters params; 
			params.set(Cornu::Parameters::MAX_SAMPLING_INTERVAL, s.sampling_interval);
			params.set(Parameters::ERROR_COST, 5);
			params.set(Parameters::MAX_RESCALE, std::numeric_limits<double>::infinity());
			params.set(Parameters::CLOSEDNESS_THRESHOLD, std::numeric_limits<double>::epsilon()); // only close the curve if it is closed
			//params.set(Parameters::SHORTNESS_COST, 5);
			//params.setAlgorithm(AlgorithmStage::PRIMITIVE_FITTING, 1);
			params.setAlgorithm(AlgorithmStage::RESAMPLING, 2); // Use my version of default resampler which use MAX_SAMPLING_INTERVAL in the correct scaled way
			params.setAlgorithm(AlgorithmStage::CURVE_CLOSING, 1);
			fitter.setParams(params);

			VectorC<Eigen::Vector2d> pts(s.points.size(), NOT_CIRCULAR);

			for (size_t i = 0; i < s.points.size(); i++) {
				pts[i] = Eigen::Vector2d(s.points[i].first.x, s.points[i].first.y);
			}

			//pass it to the fitter and process it
			fitter.setOriginalSketch(new Polyline(pts));
			fitter.run(); //to see why this prints debugging output, look at DebuggingTestImpl in Test.cpp

			//process the output -- count the number of primitives of each type
			PrimitiveSequenceConstPtr output = fitter.finalOutput();

			// From combiner.cpp
			smart_ptr<const AlgorithmOutput<GRAPH_CONSTRUCTION> > graph = fitter.output<GRAPH_CONSTRUCTION>();
			const std::vector<int> &path = fitter.output<PATH_FINDING>()->path;
			std::vector<int> _primIdcs;
			std::vector<int> _continuities; //continuity[i] is between curves i and i + 1
			std::vector<FitPrimitive> _primitives(fitter.output<PRIMITIVE_FITTING>()->primitives);
			is_closed = fitter.output<CURVE_CLOSING>()->closed;

			if (path.size() == 0) {
				cut.emplace_back(s);
				return cut;
			}

			for(int i = 0; i < (int)path.size(); ++i)
			{
				_primIdcs.push_back(graph->edges[path[i]].startVtx);
				_continuities.push_back(graph->edges[path[i]].continuity);
			}
			if (!is_closed && _primIdcs.back() != graph->edges[path.back()].endVtx)
				_primIdcs.push_back(graph->edges[path.back()].endVtx);

			///////////////// Original sample output ////////////////////
			VectorC<CurvePrimitivePtr> _curves = VectorC<CurvePrimitivePtr>((int)_primIdcs.size(), is_closed ? CIRCULAR : NOT_CIRCULAR);
			VectorC<std::pair<int, int> > _curveRanges = VectorC<std::pair<int, int> >((int)_primIdcs.size(), _curves.circular());

			for(int i = 0; i < (int)_primIdcs.size(); ++i) {
				_curveRanges[i] = std::make_pair(_primitives[_primIdcs[i]].startIdx, _primitives[_primIdcs[i]].endIdx);
				_curves[i] = _primitives[_primIdcs[i]].curve->clone();
			}

			/*if (is_closed)
				_curveRanges.back().second = _curveRanges.front().first;*/

			// trim curves and ranges
			const int sampledPts = fitter.output<RESAMPLING>()->output->pts().size();
			for(int i = 0; i < (int)_continuities.size(); ++i)
			{
				if(_continuities[i] <= 0)
					continue;

				//trim
				Eigen::Vector2d trimPt = 0.5 * (_curves[i]->endPos() + _curves[i + 1]->startPos());
				if(!_primitives[_primIdcs[i]].isFixed())
					_curves[i]->trim(0, _curves[i]->project(trimPt));
				if(_curves.circular() || !_primitives[_primIdcs[i + 1]].isFixed())
					_curves[i + 1]->trim(_curves[i + 1]->project(trimPt), _curves[i + 1]->length());

				//the ranges over which error is computed should not overlap
				if(_continuities[i] == 1)
					std::swap(_curveRanges[i].second, _curveRanges[i + 1].first); //they overlap by 1 originally
				if(_continuities[i] == 2) //they overlap by 2
				{
					_curveRanges[i + 1].first = (_curveRanges[i + 1].first + 2) % sampledPts;
					_curveRanges[i].second = (_curveRanges[i].second + sampledPts - 2) % sampledPts;
				}
			}

			VectorC<Eigen::Vector2d> _pts(fitter.output<RESAMPLING>()->output->pts());
			auto add_to_stroke = [&](int from, int to, int curv_ind, Sketch &cur_stroke) {
				if(from < 0 || to >= (int)_pts.size())
					return ;

				for(VectorC<Eigen::Vector2d>::Circulator circ = _pts.circulator(from); ; ++circ)
				{
					int idx = circ.index();
					bool last = (idx == to);
					const Eigen::Vector2d &pt = _pts.flatAt(idx);
					double s = _curves[curv_ind]->project(pt);
					Eigen::Vector2d tangent = _curves[curv_ind]->der2(s);
					cur_stroke.points.emplace_back(SketchPoint(pt.x(), pt.y(), tangent.x(), tangent.y(), _curves[curv_ind]->curvature(s)), 0);
					if(last)
						break;
				}
			};

			auto find_max_curv = [&](int from, int to, int curv_ind, int &max_ind)->double {
				double max_curv = -1;
				for(VectorC<Eigen::Vector2d>::Circulator circ = _pts.circulator(from); ; ++circ)
				{
					int idx = circ.index();
					bool last = (idx == to);

					const Eigen::Vector2d &pt = _pts.flatAt(idx);

					double s = _curves[curv_ind]->project(pt);
					if (std::abs(_curves[curv_ind]->curvature(s)) > max_curv) {
						max_curv = std::abs(_curves[curv_ind]->curvature(s));
						max_ind = idx;
					}

					if(last)
						break;
				}
				return max_curv;
			};

			auto add_break = [&](int from, int to, int curv_ind, Sketch &cur_stroke) {
				int max_ind;
				double max_curv = find_max_curv(from, to, curv_ind, max_ind);
				if (max_curv > max_curv_threshold) {
					add_to_stroke(from, max_ind, curv_ind, cur_stroke);

					// Indicate it's a curvature cut
					cur_stroke.points.back().second = 2;

					cut.emplace_back(Sketch());
					cut.back().group_ind = s.group_ind;
					cut.back().stroke_ind = ++max_id;
					cutindex_to_index[cut.back().stroke_ind] = s.stroke_ind;
					add_to_stroke(max_ind, to, curv_ind, cut.back());

					// Indicate it's a curvature cut
					cut.back().points.front().second = 2;
				}
				else {
					add_to_stroke(from, to, curv_ind, cur_stroke);
				}
			};

			///////////////// Cornucopia output for smoothing ////////////////////
			auto add_to_stroke_prim = [&](CurvePrimitiveConstPtr const primitive, double _s0, double _s1,
				Sketch &cur_stroke) {
				const Eigen::Vector2d &pt0 = primitive->pos(_s0);
				Eigen::Vector2d tangent = primitive->der2(_s0);
				tangent.normalize();
				cur_stroke.points.emplace_back(SketchPoint(pt0.x(), pt0.y(), tangent.x(), tangent.y(), primitive->curvature(_s0)), 0);

				for(double _s = _s0 + s.sampling_interval; _s + s.sampling_interval < _s1; _s += s.sampling_interval) //skip pixels 
				{
					const Eigen::Vector2d &pt = primitive->pos(_s);
					Eigen::Vector2d tangent = primitive->der2(_s);
					tangent.normalize();

					cur_stroke.points.emplace_back(SketchPoint(pt.x(), pt.y(), tangent.x(), tangent.y(), primitive->curvature(_s)), 0);
				}

				const Eigen::Vector2d &pt1 = primitive->pos(_s1);
				tangent = primitive->der2(_s1);
				tangent.normalize();
				cur_stroke.points.emplace_back(SketchPoint(pt1.x(), pt1.y(), tangent.x(), tangent.y(), primitive->curvature(_s1)), 0);

				cur_stroke.sampling_interval = s.sampling_interval;
			};

			auto find_max_curv_prim = [&](CurvePrimitiveConstPtr const primitive, double _s0, double _s1,
				double &max_s)->double {
				double max_curv = -1;

				for(double _s = _s0; _s < _s1; _s += s.sampling_interval) //skip pixels 
				{
					if (std::abs(primitive->curvature(_s)) > max_curv) {
						max_curv = std::abs(primitive->curvature(_s));
						max_s = _s;
					}
				}

				return max_curv;
			};

			auto add_break_prim = [&](CurvePrimitiveConstPtr const primitive, double _s0, double _s1, Sketch &cur_stroke) {
				double max_s;
				double max_curv = find_max_curv_prim(primitive, _s0, _s1, max_s);

				if (max_curv > max_curv_threshold) {
					add_to_stroke_prim(primitive, _s0, max_s, cur_stroke);
					
					// Indicate it's a curvature cut
					cur_stroke.points.back().second = 2;

					cut.emplace_back(Sketch());
					cut.back().group_ind = s.group_ind;
					cut.back().stroke_ind = ++max_id;//s.stroke_ind;//
					cutindex_to_index[cut.back().stroke_ind] = s.stroke_ind;
					add_to_stroke_prim(primitive, max_s, _s1, cut.back());

					// Indicate it's a curvature cut
					cut.back().points.front().second = 2;
				}
				else {
					add_to_stroke_prim(primitive, _s0, _s1, cur_stroke);
				}
			};
			//////////////////////////////////////

			// Use the fitting from cornucopia
			if (!output_orig) {
				PrimitiveSequenceConstPtr prim_output = fitter.finalOutput();
				if (prim_output->primitives().size() > 0) {
					cut.emplace_back(Sketch());
					cut.back().group_ind = s.group_ind;
					cut.back().stroke_ind = s.stroke_ind;

					add_break_prim(prim_output->primitives()[0], 0, prim_output->primitives()[0]->length(), cut.back());

					for (int i = 0; i + 1 < (int)prim_output->primitives().size(); ++i) {
						if (_continuities[i] == 0) {
							cut.emplace_back(Sketch());
							cut.back().group_ind = s.group_ind;
							cut.back().stroke_ind = ++max_id;//s.stroke_ind;//
							cutindex_to_index[cut.back().stroke_ind] = s.stroke_ind;
						}

						add_break_prim(prim_output->primitives()[i + 1], 0, prim_output->primitives()[i + 1]->length(), cut.back());
					}
				}
			}
			else {
				if (_curveRanges.size() > 0) {
					cut.emplace_back(Sketch());
					cut.back().group_ind = s.group_ind;
					cut.back().stroke_ind = s.stroke_ind;
					add_break(_curveRanges[0].first, _curveRanges[0].second, 0, cut.back());
					for (int i = 0; i + 1 < (int)_curveRanges.size(); ++i) {
						if (_continuities[i] == 0) {
							cut.emplace_back(Sketch());
							cut.back().group_ind = s.group_ind;
							cut.back().stroke_ind = ++max_id;
							cutindex_to_index[cut.back().stroke_ind] = s.stroke_ind;
						}
						add_break(_curveRanges[i + 1].first, _curveRanges[i + 1].second, i + 1, cut.back());
					}
				}
			}
		}

		std::vector<Sketch> dezero;
		dezero.reserve(cut.size());
		for (auto const &s : cut) {
			if (!s.points.empty() &&
				s.totalLen() > std::numeric_limits<double>::epsilon()) {
				dezero.emplace_back(s);
			}
		}

		// Find if there's any pair of curvature cut points that are too close to each other.
		// If yes, split it to its neighbors.
		// Note that this kind of segments can only in the middle of an input stroke, 
		// since the start and end of the entire stroke are not curvature cut points.
		{
			for (size_t i = 1; i + 1 < dezero.size(); i++) {
				double length = dezero[i].totalLen();
				if (/*dezero[i].points.front().second == 2 &&
					dezero[i].points.back().second == 2 &&*/
					length < curvature_cut_min_length) {
					dezero[i].points.front().second = 0;
					dezero[i].points.back().second = 0;

					int prev_ind = i - 1;
					while (dezero[prev_ind].points.size() == 0 && prev_ind >= 0)
						prev_ind--;

					size_t num_p = dezero[i].points.size();

					if (num_p > 2) {
						dezero[prev_ind].points.insert(dezero[prev_ind].points.end(), 
							dezero[i].points.begin() + 1, dezero[i].points.begin() + (num_p >> 1));
						dezero[i + 1].points.insert(dezero[i + 1].points.begin(), 
							dezero[i].points.begin() + (num_p >> 1), dezero[i].points.end() - 1);
					}
					else {
						// Interpolate the mid point
						Sketch seg;
						seg.points.emplace_back(dezero[i].points.front());
						seg.points.emplace_back(dezero[i].points.back());
						seg.reparameterize(seg.totalLen(), 3);

						dezero[prev_ind].points.insert(dezero[prev_ind].points.end(), seg.points[1]);
						dezero[i + 1].points.insert(dezero[i + 1].points.begin(), seg.points[1]);
					}

					// Keep the curvature cut flag
					dezero[prev_ind].points.back().second = 2;//dezero[i].points.front().second;
					dezero[i + 1].points.front().second = 2;//dezero[i].points.back().second;

					dezero[i].points.clear();
				}
			}
		}

		std::vector<Sketch> declose;
		for (auto const &s : dezero) {
			if (!s.points.empty() &&
				s.totalLen() > std::numeric_limits<double>::epsilon()) {
				declose.emplace_back(s);
			}
		}

		// Find if there's any curvature cut point which is significantly larger than its neighbors.
		{
			auto compute_median_curvature = [](std::vector<double> &curv)->double {
				std::nth_element(curv.begin(), curv.begin() + curv.size()/2, curv.end());
				return curv[curv.size()/2];
			};
			auto compute_median_absolute_deviation = [](std::vector<double> const &curv, double median)->double {
				std::vector<double> abs_dev = curv;
				for (double &v : abs_dev) v = std::abs(v - median);

				std::nth_element(abs_dev.begin(), abs_dev.begin() + abs_dev.size()/2, abs_dev.end());
				return abs_dev[abs_dev.size()/2];
			};

			auto compute_mean = [](std::vector<double> &curv)->double {
				double sum = 0;
				for (double &c : curv)
					sum += c;
				return (curv.empty()) ? 0 : sum / curv.size();
			};
			auto compute_standard_deviation = [&](std::vector<double> const &curv, double mean)->double {
				std::vector<double> x2 = curv;
				for (double &v : x2) v *= v;

				double e_x2 = compute_mean(x2);
				return (curv.empty()) ? 0 : std::sqrt(e_x2 - mean * mean);
			};

			std::vector<double> curv;
			for (auto const &seg : declose) {
				for (auto const &p: seg.points) {
					curv.emplace_back(std::abs(p.first.curvature));
				}
			}

			for (size_t i = 0; i + 1 < declose.size(); i++) {
				if (!declose[i].points.empty() &&
					declose[i].points.back().second == 2) {
					std::vector<double> local_curv;
					for (auto const &p: declose[i].points) {
						local_curv.emplace_back(std::abs(p.first.curvature));
					}

					double median = compute_median_curvature(curv);

					// Debug
					if (median < std::numeric_limits<double>::epsilon())
						median = 1e-6;
					//

					//double mad = compute_median_absolute_deviation(curv, median);

					//double mean = compute_mean(curv);
					//double stddev = compute_standard_deviation(curv, mean);

					double max_curv = *std::max_element(local_curv.begin(), local_curv.end());

					double ratio = std::abs(max_curv - median) / median;

#ifdef DEBUG_PRINTOUT
					std::cout << "\t curv ratio: " << ratio << std::endl;
#endif

					// Not a significantly large peak
					if (ratio < max_curvature_ratio) {
						declose[i].points.back().second = 0;
						declose[i + 1].points.insert(declose[i + 1].points.begin(), declose[i].points.begin(), declose[i].points.end());
						declose[i].points.clear();
					}
				}
			}
		}

		cut.clear();
		for (auto const &s : declose) {
			if (!s.points.empty() &&
				s.totalLen() > std::numeric_limits<double>::epsilon()) {
				cut.emplace_back(s);
			}
		}

		if (is_closed && 
			(cut.front().points.front().first - cut.back().points.back().first).Length() > std::numeric_limits<double>::epsilon()) {
			SketchPoint p = cut.front().points.front().first;
			cut.back().points.emplace_back(p, 0);
		}

		return cut;
	};
	Capture cut_capture;

	for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
		// Degenerated case: single point
		Sketch s = capture.sketchedPolylines[i];
		s.reindex();
		if (s.points.size() == 1)
			continue;

		bool is_closed;

		std::vector<Sketch> cut_strokes = cut_stroke(capture.sketchedPolylines[i], is_closed);

#ifdef DEBUG_PRINTOUT
		std::cout << "\t" << capture.sketchedPolylines[i].stroke_ind << ": " << is_closed << std::endl;
#endif

		std::vector<bool> to_keep;
		to_keep.resize(cut_strokes.size(), true);

		if (cut_strokes.size() == 0) continue;

		if (!is_closed) {
			to_keep.front() = !is_too_short(capture.sketchedPolylines[i], cut_strokes.front());
			to_keep.back() = !is_too_short(capture.sketchedPolylines[i], cut_strokes.back());
		}

		c_record[s.stroke_ind] = std::vector<size_t>();
		std::unordered_set<int> cut_stroke_indices;
		for (size_t j = 0; j < cut_strokes.size(); j++) {
			if (to_keep[j]) {
				cut_capture.sketchedPolylines.emplace_back(cut_strokes[j]);
				//cut_capture.sketchedPolylines.emplace_back(cut_strokes[j]);

				// Avoid deleting updated indices later
				cut_stroke_indices.insert(-1 * (int)cut_strokes[j].stroke_ind);

				c_record[s.stroke_ind].emplace_back(cut_strokes[j].stroke_ind);
			}
		}

		if (!cut_stroke_indices.empty())
			cut_record.emplace_back(cut_stroke_indices);
	}

	// Recompute sampling rate
	for (size_t i = 0; i < cut_capture.sketchedPolylines.size(); i++) {
		cut_capture.sketchedPolylines[i].sampling_interval = 0;

		for (size_t j = 0; j + 1 < cut_capture.sketchedPolylines[i].points.size(); j++) {
			cut_capture.sketchedPolylines[i].sampling_interval += 
				(cut_capture.sketchedPolylines[i].points[j + 1].first - cut_capture.sketchedPolylines[i].points[j].first).Length();
		}
		cut_capture.sketchedPolylines[i].sampling_interval = (cut_capture.sketchedPolylines[i].points.size() >= 2) ? 
			cut_capture.sketchedPolylines[i].sampling_interval/cut_capture.sketchedPolylines[i].points.size() : 0;
		cut_capture.sketchedPolylines[i].reparameterize(cut_capture.sketchedPolylines[i].sampling_interval);

		//if (i % 2 == 0)
		smooth_stroke(cut_capture.sketchedPolylines[i]);

		/*cut_capture.sketchedPolylines[i].reparameterize(cut_capture.sketchedPolylines[i].sampling_interval);
		smooth_stroke(cut_capture.sketchedPolylines[i]);*/
	}

	return cut_capture;
}

void StrokeCut::prepare_cornucopia(Capture &capture) const {
	// Compute sampling rate for cornucopia
	for (size_t i = 0; i < capture.sketchedPolylines.size(); i++) {
		double reparam_step_length = 5.;//std::max(5., capture.sketchedPolylines[i].sampling_interval);
		double step_length = (capture.getSketchedPolyline(i).totalLen() / min_num_sample_per_stroke);
		step_length = std::min(step_length, reparam_step_length);
		step_length = std::max(step_length, capture.sketchedPolylines[i].sampling_interval);
		//capture.getSketchedPolyline(i).reparameterize(step_length);

		step_length = std::max(step_length, capture.getSketchedPolyline(i).totalLen() / max_num_sample_per_stroke);

		size_t num_steps = capture.getSketchedPolyline(i).totalLen() / step_length;

		capture.getSketchedPolyline(i).sampling_interval = step_length; // temporarily set this value for cut_cornucopia func
	}
}

void StrokeCut::update_cut_record(int from, int to) {
	for (auto &s : cut_record) {
		if (s.count(from) == 0) continue;
		s.erase(from);
		s.emplace(to);
		break;
	}
}

std::unordered_set<int> StrokeCut::get_cut_record(int ind) const {
	std::unordered_set<int> record;
	for (auto const &s : cut_record) {
		if (s.count(ind) == 0) continue;
		
		record = s;
		break;
	}

	return record;
}
