#define _USE_MATH_DEFINES
#include <cmath>

#include <future>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <deque>
#include <unordered_set>
#include <unordered_map>
#include <glm/gtx/norm.hpp>

#include "Parameterization.h"
#include "Utils.h"
#include "SvgUtils.h"

Parameterization::Parameterization(const Context& context): context(context) {}

void Parameterization::parameterize(Input* input) {
	map_clusters<int>(*input, [&](Cluster& c) { parameterize_cluster(&c); return 0; });
}

void Parameterization::isolines_svg(std::ostream& os, const Input& input) {
	input.cluster_svg(os, [&](std::ostream& os) {
		for (auto& kv : input.clusters) {
			for (auto& xsec : kv.second.xsecs) {
				for (auto& connection : xsec.connections) {
					if (std::abs(int(connection.a_idx) - int(connection.b_idx)) == 1) {
						std::stringstream ss;
						ss << "#";
						if (connection.weight > 0.5) {
							ss << std::setfill('0') << std::setw(2) << std::hex << int(map(connection.weight, 0.5, 1.0, 255.0, 0.0)); // red
							ss << std::setfill('0') << std::setw(2) << std::hex << int(map(connection.weight, 0.5, 1.0, 255.0 / 2.0, 255.0)); // green
							ss << "00"; // blue
						}
						else {
							ss << "FF", // red
							ss << std::setfill('0') << std::setw(2) << std::hex << int(map(connection.weight, 0.0, 0.5, 0.0, 255.0 / 2.0)); // green
							ss << "00"; // blue
						}
						SVG::line(
							os,
							xsec.points[connection.a_idx].point.x * input.thickness,
							xsec.points[connection.a_idx].point.y * input.thickness,
							xsec.points[connection.b_idx].point.x * input.thickness,
							xsec.points[connection.b_idx].point.y * input.thickness,
							0.5,
							ss.str()
						);
					}
				}
			}
		}
	});
}

void Parameterization::debug_svg(std::ostream& os, const Input& input) {
	input.cluster_svg(os, [&](std::ostream& os) {
		for (auto& line : debug_lines) {
			line.from *= input.thickness;
			line.to *= input.thickness;
			SVG::line(os, line.from.x, line.from.y, line.to.x, line.to.y, 1.0, line.color);
		}
	});
}

void Parameterization::add_debug_line(Parameterization::DebugLine line) {
	std::lock_guard<std::mutex> lock(viz_lock);
	debug_lines.push_back(line);
}

void Parameterization::parameterize_cluster(Cluster* cluster) {
	cluster->xsecs = orthogonal_xsecs(*cluster);
	std::vector<std::vector<double>> prev_u;

	auto record_current_u = [&]() -> void {
		prev_u.clear();
		for (auto& stroke : cluster->strokes) {
			prev_u.push_back(stroke.u);
		}
	};

	auto converged = [&]() -> bool {
		bool ok = true;
		for (size_t stroke = 0; ok && stroke < prev_u.size(); ++stroke) {
			for (size_t i = 0; ok && i < prev_u[stroke].size(); ++i) {
				double diff = std::abs(prev_u[stroke][i] - cluster->strokes[stroke].u[i]);
				if (i == prev_u[stroke].size() - 1) {
					ok = diff * diff < 0.5 * 0.5 * glm::distance2(cluster->strokes[stroke].points[i - 1], cluster->strokes[stroke].points[i]);
				}
				else {
					ok = diff * diff < 0.5 * 0.5 * glm::distance2(cluster->strokes[stroke].points[i], cluster->strokes[stroke].points[i + 1]);
				}
			}
		}
		return ok;
	};

	for (int it = 0; it < 5; ++it) {
		record_current_u();
		ensure_connected(cluster);
		params_from_xsecs(cluster, it == 0, nullptr);

		if (converged()) {
			break;
		}
		else {
			// Adjust connection weights
			double total_len = cluster->max_u();
			for (auto& xsec : cluster->xsecs) {
				if (xsec.connector) continue;
				auto get_u = [&](size_t idx) -> double {
					auto& pt = xsec.points[idx];
					double mix = pt.i - std::floor(pt.i);
					return (1. - mix) * cluster->strokes[pt.stroke_idx].u[std::floor(pt.i)] +
						mix * cluster->strokes[pt.stroke_idx].u[std::ceil(pt.i)];
				};
				for (auto& connection : xsec.connections) {
					double diff = get_u(connection.a_idx) - get_u(connection.b_idx);
					connection.weight = gaussian(diff, total_len / 30., 0.);
				}
			}
		}
	}

	for (int it = 0; it < 5; ++it) {
		record_current_u();
		cluster->xsecs = xsecs_from_params(*cluster);
		ensure_connected(cluster);
		params_from_xsecs(cluster, false, nullptr);

		if (converged()) {
			break;
		}
	}

	record_current_u();
	check_periodic(cluster);
	if (cluster->periodic) {
		auto old_u = [&](size_t stroke, double i) {
			double mix = i - std::floor(i);
			return (1. - mix) * prev_u[stroke][std::floor(i)] + mix * prev_u[stroke][std::ceil(i)];
		};
		auto is_cut = [&](const Cluster::XSec& xsec) {
			for (auto& conn : xsec.connections) {
				if (std::abs(
					old_u(xsec.points[conn.a_idx].stroke_idx, xsec.points[conn.a_idx].i) -
					old_u(xsec.points[conn.b_idx].stroke_idx, xsec.points[conn.b_idx].i))
				> 5.) {
					return true;
				}
			}
			return false;
		};
		auto cut_it = std::find_if(cluster->xsecs.begin(), cluster->xsecs.end(), is_cut);
		if (cut_it != cluster->xsecs.end()) {
			auto cut = *cut_it;
			ensure_connected(cluster, &cut);
			params_from_xsecs(cluster, false, &cut);
		}
	}
	cluster->xsecs = xsecs_from_params(*cluster);
}

void Parameterization::params_from_xsecs(Cluster* cluster, bool initial, Cluster::XSec* cut) {
	GRBModel model(context.grb);

	std::vector<std::vector<GRBVar>> param_vars;
	for (auto& stroke : cluster->strokes) {
		param_vars.emplace_back();
		auto& vars = param_vars.back();
		for (size_t i = 0; i < stroke.points.size(); ++i) {
			std::string name = std::string("u_{") + std::to_string(param_vars.size()-1) + "," + std::to_string(i) + "}";
			vars.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS, name));
		}
	}

	auto points_cross_cut = [&](size_t stroke, int a, int b) {
		for (auto& pt : cut->points) {
			//if (pt.stroke_idx == stroke && std::ceil(pt.i) >= a && std::floor(pt.i) <= b) {
			if (pt.stroke_idx == stroke && pt.i >= a && pt.i <= b) {
				return true;
			}
		}
		return false;
	};

	std::vector<GRBLinExpr> velocity_terms;
	std::vector<GRBLinExpr> alignment_terms;
	for (auto& xsec : cluster->xsecs) {
		glm::dvec2 tangent = xsec.avg_tangent();

		auto connection_crosses_cut = [&](const Cluster::XSecConnection conn) {
			for (size_t idx : { conn.a_idx, conn.b_idx }) {
				size_t stroke = xsec.points[idx].stroke_idx;
				double i = xsec.points[idx].i;

				for (auto& pt : cut->points) {
					if (pt.stroke_idx == stroke && std::ceil(pt.i) >= std::ceil(i) && std::floor(pt.i) <= std::floor(i)) {
						return true;
					}
				}
			}
			return false;
		};

		// 1. Add velocity term
		{
			velocity_terms.emplace_back();
			auto& velocity_term = velocity_terms.back();

			double total_weight = 0.0;
			for (size_t i = 0; i < xsec.points.size(); ++i) {
				auto& point = xsec.points[i];
				double weight = xsec.distance_weight(i) + 1; // Regularize to avoid zero weights on very close strokes
				double coefficient = weight * glm::dot(point.tangent, tangent) / glm::length(point.to_next);
				if (std::isnan(coefficient)) {
					throw "coefficient NaN";
				}
				if (std::isinf(coefficient)) {
					//throw "coefficient inf";
					continue;
				}

				int a = std::floor(point.i);
				int b = std::ceil(point.i);
				double mix = point.i - a;
				if (a == b) {
					if (b == cluster->strokes[point.stroke_idx].points.size() - 1) {
						--a;
					}
					else {
						++b;
					}
				}

				if (cut && points_cross_cut(point.stroke_idx, a, b)) continue;

				total_weight += weight;
				velocity_term += coefficient * (1 - mix) *
					(param_vars[point.stroke_idx][b] - param_vars[point.stroke_idx][a]);

			}

			if (total_weight > 0) {
				velocity_term /= total_weight;
				velocity_term -= 1.0;
			}
			else {
				velocity_terms.pop_back();
			}
		}

		// 2. Add alignment terms
		if (!cut || !std::any_of(xsec.connections.begin(), xsec.connections.end(), connection_crosses_cut)) {
			for (auto& connection : xsec.connections) {
				auto& pt_a = xsec.points[connection.a_idx];
				auto& pt_b = xsec.points[connection.b_idx];

				double mix_a = pt_a.i - std::floor(pt_a.i);
				GRBLinExpr u_a = (1 - mix_a) * param_vars[pt_a.stroke_idx][std::floor(pt_a.i)] +
					mix_a * param_vars[pt_a.stroke_idx][std::ceil(pt_a.i)];

				double mix_b = pt_b.i - std::floor(pt_b.i);
				GRBLinExpr u_b = (1 - mix_b) * param_vars[pt_b.stroke_idx][std::floor(pt_b.i)] +
					mix_b * param_vars[pt_b.stroke_idx][std::ceil(pt_b.i)];

				double proj_dist = glm::dot(pt_a.point - pt_b.point, tangent);
				alignment_terms.push_back(connection.weight / (xsec.connections.size()) * ((u_a - u_b) - proj_dist));
			}
		}
	}

	// 3. Enforce monotonicity
	for (size_t stroke_idx = 0; stroke_idx < cluster->strokes.size(); ++stroke_idx) {
		auto& stroke = cluster->strokes[stroke_idx];
		for (size_t i = 0; i < stroke.points.size() - 1; ++i) {
			if (cut && points_cross_cut(stroke_idx, i, i + 1)) continue;

			model.addConstr(
				param_vars[stroke_idx][i + 1] - param_vars[stroke_idx][i] >=
					0.5 * glm::distance(stroke.points[i + 1], stroke.points[i]));
		}
	}

	// 4. Boundary
	if (cut) {
		for (auto& pt : cut->points) {
			double dist_from_cut = glm::distance(pt.point, point(cluster->strokes[pt.stroke_idx].points, std::ceil(pt.i)));
			model.addConstr(param_vars[pt.stroke_idx][std::ceil(pt.i)] == dist_from_cut);
		}
	} else {
		model.addConstr(param_vars[0][0] == 0.0);
	}

	GRBQuadExpr objective;
	if (initial) {
		double scale = 1e-5;
		if (alignment_terms.size() > 6 * velocity_terms.size()) {
			scale = 1e-6;
		}
		objective = l2_norm_sq(&model, velocity_terms) + scale * l1_norm(&model, alignment_terms);
	}
	else {
		objective = l2_norm_sq(&model, velocity_terms) + l2_norm_sq(&model, alignment_terms);
	}

	try {
		model.setObjective(objective, GRB_MINIMIZE);
	}
	catch (GRBException e) {
		std::cout << "Error code = " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;

		model.write("model.lp");

		throw e;
	}

	context.optimize_model(&model);

	double min_u = std::numeric_limits<double>::infinity();
	for (size_t stroke_idx = 0; stroke_idx < cluster->strokes.size(); ++stroke_idx) {
		auto& stroke = cluster->strokes[stroke_idx];
		for (size_t i = 0; i < stroke.points.size(); ++i) {
			stroke.u[i] = param_vars[stroke_idx][i].get(GRB_DoubleAttr_X);
			min_u = std::min(min_u, stroke.u[i]);
		}
	}
	for (auto& stroke : cluster->strokes) {
		for (size_t i = 0; i < stroke.u.size(); ++i) {
			stroke.u[i] -= min_u;
		}
	}
}

std::vector<Cluster::XSec> Parameterization::orthogonal_xsecs(const Cluster& cluster, double angle_tolerance) {
	const double TARGET_ANGLE = 20. / 180. * M_PI;
	const double MIN_GAP_CUTOFF = 5.;
	const double MAX_GAP_CUTOFF = 20.;

	std::vector<Cluster::XSec> xsecs;

	for (size_t stroke = 0; stroke < cluster.strokes.size(); ++stroke) {
		for (size_t i = 0; i < cluster.strokes[stroke].points.size(); ++i) {
			xsecs.push_back(orthogonal_xsec_at(cluster, stroke, i));
		}
	}

	// Local filtering

	// 1. Angle filtering
	{
		struct AngleSample {
			glm::dvec2 mid;
			double angle;
		};
		std::vector<AngleSample> samples;
		for (auto& xsec : xsecs) {
			for (size_t i = 1; i < xsec.points.size(); ++i) {
				samples.push_back({
					(xsec.points[i - 1].point + xsec.points[i].point) * 0.5,
					std::acos(glm::dot(xsec.points[i - 1].tangent, xsec.points[i].tangent)),
				});
			}
		}
		for (auto& xsec : xsecs) {
			std::deque<AngleSample> neighbourhood;
			double r_sq = std::pow(std::max(30.0, glm::distance(xsec.points.front().point, xsec.points.back().point)), 2);

			auto sample_within_radius = [&](const AngleSample& sample) {
				for (auto& point : xsec.points) {
					if (glm::distance2(point.point, sample.mid) <= r_sq) {
						return true;
					}
				}
				return false;
			};

			// TODO k-d tree lookup?
			for (auto& sample : samples) {
				if (sample_within_radius(sample)) {
					neighbourhood.push_back(sample);
				}
			}

			double max_angle = TARGET_ANGLE;
			if (!neighbourhood.empty()) {

				// Sort by angle descending
				std::sort(neighbourhood.begin(), neighbourhood.end(), [](const AngleSample& a, const AngleSample& b) { return a.angle > b.angle; });

				// Remove samples until the `angle_tolerence` percentile is below the target value
				if (neighbourhood.back().angle <= TARGET_ANGLE) {
					max_angle = neighbourhood.front().angle;
					while (!neighbourhood.empty() && neighbourhood[int(angle_tolerance * (neighbourhood.size() - 1))].angle > TARGET_ANGLE) {
						neighbourhood.pop_front();
						max_angle = neighbourhood.front().angle;
					}
				}
			}

			// Remove connections above the threshold to the right of the center
			for (int i = xsec.center_idx + 1; i < xsec.points.size(); ++i) {
				if (std::acos(glm::dot(xsec.points[i - 1].tangent, xsec.points[i].tangent)) > max_angle) {
					xsec.points = std::vector<Cluster::XSecPoint>(xsec.points.begin(), xsec.points.begin() + i);
					break;
				}
			}
			// Remove connections above the threshold to the left of the center
			for (int i = xsec.center_idx - 1; i >= 0; --i) {
				if (std::acos(glm::dot(xsec.points[i + 1].tangent, xsec.points[i].tangent)) > max_angle) {
					xsec.points = std::vector<Cluster::XSecPoint>(xsec.points.begin() + i + 1, xsec.points.end());
					xsec.center_idx -= i + 1;
					break;
				}
			}
		}
	}

	// 2. Gap filtering
	{
		struct GapSample {
			glm::dvec2 left;
			glm::dvec2 right;
			double gap_sq;
		};
		std::vector<GapSample> samples;
		for (auto& xsec : xsecs) {
			for (size_t i = 1; i < xsec.points.size(); ++i) {
				samples.push_back({
					xsec.points[i - 1].point,
					xsec.points[i].point,
					glm::distance2(xsec.points.front().point, xsec.points.back().point),
				});
			}
		}
		for (auto& xsec : xsecs) {
			if (xsec.points.size() == 1) continue;

			std::unordered_set<double> neighbourhood;
			double r_sq = std::max(200.0*200.0, glm::distance2(xsec.points.front().point, xsec.points.back().point));

			auto sample_within_radius = [&](const GapSample& sample) {
				for (auto& point : xsec.points) {
					for (auto& other : { sample.left, sample.right }) {
						if (glm::distance2(point.point, other) <= r_sq) {
							return true;
						}
					}
				}
				return false;
			};

			// TODO k-d tree lookup?
			for (auto& sample : samples) {
				if (sample_within_radius(sample)) {
					neighbourhood.insert(sample.gap_sq);
				}
			}

			if (neighbourhood.empty()) continue;

			// Get median gap
			std::vector<double> gaps(neighbourhood.begin(), neighbourhood.end());
			size_t median_offset = gaps.size() * 0.75;
			std::nth_element(
				gaps.begin(),
				gaps.begin() + median_offset,
				gaps.end());
			double median_gap_sq = gaps[median_offset];
			double gap_cutoff = std::min(MAX_GAP_CUTOFF*MAX_GAP_CUTOFF, std::max(MIN_GAP_CUTOFF*MIN_GAP_CUTOFF, 1.2*1.2*median_gap_sq));

			// Remove connections above the threshold to the right of the center
			for (int i = xsec.center_idx + 1; i < xsec.points.size(); ++i) {
				if (glm::distance2(xsec.points[i - 1].point, xsec.points[i].point) > gap_cutoff) {
					xsec.points = std::vector<Cluster::XSecPoint>(xsec.points.begin(), xsec.points.begin() + i);
					break;
				}
			}
			// Remove connections above the threshold to the left of the center
			for (int i = xsec.center_idx - 1; i >= 0; --i) {
				if (glm::distance2(xsec.points[i + 1].point, xsec.points[i].point) > gap_cutoff) {
					xsec.points = std::vector<Cluster::XSecPoint>(xsec.points.begin() + i + 1, xsec.points.end());
					xsec.center_idx -= i + 1;
					break;
				}
			}
		}
	}

	for (auto& xsec : xsecs) {
		for (size_t i = 0; i < xsec.points.size(); ++i) {
			for (size_t j = i + 1; j < xsec.points.size(); ++j) {
				xsec.connections.push_back({ i, j, 1.0 });
			}
		}
	}

	return xsecs;
}

Cluster::XSec Parameterization::orthogonal_xsec_at(const Cluster& cluster, size_t stroke, double i) {
	const double SCALE = 100.;

	auto& stroke_points = cluster.strokes[stroke].points;

	size_t a = std::floor(i);
	double mix = i - double(a);
	glm::dvec2 origin = point(stroke_points, i);
	glm::dvec2 tan = tangent(stroke_points, i);
	glm::dvec2 ortho = normal(tan);
	glm::dvec2 to_next(0.0, 0.0);
	if (a < stroke_points.size() - 1) {
		to_next = (1 - mix) * (stroke_points[a + 1] - stroke_points[a]);
	}
	else {
		to_next = stroke_points[a] - stroke_points[a - 1];
	}

	glm::dvec2 pt1 = origin - SCALE * ortho;
	glm::dvec2 pt2 = origin + SCALE * ortho;

	struct PotentialInt {
		Cluster::XSecPoint point;
		double signed_dist;
	};
	std::vector<PotentialInt> potential_ints = {
		{
			{
				stroke,
				i,
				origin,
				to_next,
				tan
			},
			0.0
		}
	};

	for (size_t j = 0; j < cluster.strokes.size(); ++j) {
		auto ints = intersections(cluster.strokes[j].points, pt1, pt2);
		for (auto& intersection : ints) {

			// Ignore origin point
			if (stroke == j && std::abs(intersection.i - i) < 1e-1) continue;

			glm::dvec2 to_next(0.0, 0.0);
			double other_i = intersection.i;
			size_t floor_other_i = std::floor(intersection.i);
			if (floor_other_i < cluster.strokes[j].points.size() - 1) {
				to_next = (1 - (other_i - double(floor_other_i))) * (cluster.strokes[j].points[floor_other_i + 1] - cluster.strokes[j].points[floor_other_i]);
			}
			else if (floor_other_i == other_i) {
				to_next = cluster.strokes[j].points[floor_other_i] - cluster.strokes[j].points[floor_other_i - 1];
			}

			potential_ints.push_back({
				{
					j,
					intersection.i,
					intersection.pt,
					to_next,
					tangent(cluster.strokes[j].points, intersection.i)
				},
				glm::dot(intersection.pt - origin, ortho)
			});
			if (glm::any(glm::isnan(potential_ints.back().point.tangent))) {
				potential_ints.pop_back();
			}
		}
	}

	std::sort(potential_ints.begin(), potential_ints.end(), [](const PotentialInt& a, const PotentialInt& b) { return a.signed_dist < b.signed_dist; });

	if (potential_ints.size() > 1) {
		// Cut off intersections after any bad tangent angle
		int first_bad_tangent = -1;
		for (int n = potential_ints.size() - 1; n >= 0; --n) {
			if (potential_ints[n].signed_dist < 0 && glm::dot(potential_ints[n].point.tangent, tan) < 0) {
				first_bad_tangent = n;
				break;
			}
		}
		int last_bad_tangent = potential_ints.size();
		for (int n = 0; n < potential_ints.size(); ++n) {
			if (potential_ints[n].signed_dist > 0 && glm::dot(potential_ints[n].point.tangent, tan) < 0) {
				last_bad_tangent = n;
				break;
			}
		}

		potential_ints = std::vector<PotentialInt>(potential_ints.begin() + (1 + first_bad_tangent), potential_ints.begin() + last_bad_tangent);
	}

	Cluster::XSec xsec = { {}, {}, 0, 0.0, false };
	for (auto& potential_int : potential_ints) {
		xsec.points.push_back(potential_int.point);
		if (xsec.points.back().stroke_idx == stroke && xsec.points.back().i == i) {
			xsec.center_idx = xsec.points.size() - 1;
		}
	}

	return xsec;
}

std::vector<Cluster::XSec> Parameterization::xsecs_from_params(const Cluster& cluster, double sample_rate, bool nonlinear) {
	double max_u = cluster.max_u();
	sample_rate = std::max(sample_rate, cluster.max_u() / 1e4);
	size_t num_xsecs = std::ceil(max_u / sample_rate);
	if (nonlinear) num_xsecs *= 3;

	std::vector<Cluster::XSec> result;
	result.reserve(num_xsecs);

	std::vector<std::vector<bool>> sampled;
	for (auto& stroke : cluster.strokes) {
		sampled.emplace_back(stroke.points.size(), false);
	}

	for (size_t n = 0; n < num_xsecs; ++n) {
		double t = double(n) / double(num_xsecs - 1);
		double u;
		if (nonlinear) {
			u = poly_in_out(t, 2.25) * max_u;
		}
		else {
			u = t * max_u;
		}
		auto xsec = xsec_at_u(cluster, u);
		if (xsec.points.size() > 0) {
			result.push_back(xsec);

			for (auto& pt : xsec.points) {
				sampled[pt.stroke_idx][std::floor(pt.i)] = true;
				sampled[pt.stroke_idx][std::ceil(pt.i)] = true;
			}
		}
	}

	// Make sure all sample points have a cross section referencing them
	for (size_t stroke_idx = 0; stroke_idx < cluster.strokes.size(); ++stroke_idx) {
		for (size_t i = 0; i < cluster.strokes[stroke_idx].u.size(); ++i) {
			if (!sampled[stroke_idx][i]) {
				auto xsec = xsec_at_u(cluster, cluster.strokes[stroke_idx].u[i]);
				if (xsec.points.empty()) {
					glm::dvec2 to_next(0.0, 0.0);
					if (i < cluster.strokes[stroke_idx].points.size() - 1) {
						to_next = cluster.strokes[stroke_idx].points[i + 1] - cluster.strokes[stroke_idx].points[i];
					}
					else {
						to_next = cluster.strokes[stroke_idx].points[i] - cluster.strokes[stroke_idx].points[i - 1];
					}
					xsec = Cluster::XSec{
						{{ stroke_idx, double(i), cluster.strokes[stroke_idx].points[i], to_next, tangent(cluster.strokes[stroke_idx].points, double(i)) }},
						{},
						0,
						cluster.strokes[stroke_idx].u[i],
						false
					};
				}
				/*for (auto& pt : xsec.points) {
					sampled[pt.stroke_idx][std::floor(pt.i)] = true;
					sampled[pt.stroke_idx][std::ceil(pt.i)] = true;
				}*/
				result.push_back(xsec);
			}
		}
	}

	return result;
}

Cluster::XSec Parameterization::xsec_at_u(const Cluster& cluster, double u) {
	Cluster::XSec xsec = { {}, {}, 0, u, false };
	double period = cluster.max_u() + 1e-1;
	for (size_t stroke_idx = 0; stroke_idx < cluster.strokes.size(); ++stroke_idx) {
		auto& stroke = cluster.strokes[stroke_idx];
		auto from = stroke.u.begin();

		// Periodic parameterizations are no longer monotinic; instead, they are composed of
		// monotonic chunks. We want to iterate over each chunk individually.
		auto find_to = [](std::vector<double>::const_iterator it, std::vector<double>::const_iterator end) {
			do {
				++it;
			} while (it != end && *it > *(it - 1));
			return it;
		};

		for (auto to = from; from != stroke.u.end(); from = to) {
			to = find_to(from, stroke.u.end());

			auto it_gt = std::lower_bound(from, to, u);
			auto it_le = it_gt;
			if (it_gt == from) {
				// This will only be a problem if we're at the beginning of a stroke.
				// Periodic strokes will be special-cased later.
				if (it_gt == stroke.u.begin()) {
					continue;
				}
				else {
					--it_le;
				}
			}
			else if (it_gt != to && *it_gt - u < 1e-6) {
				++it_gt;
			}
			else {
				--it_le;
			}
			if (it_gt == to) {
				// This is only a problem if we're at the end of a stroke. Periodic
				// strokes will be special-cased later.
				if (to == stroke.u.end()) {
					if (*it_le == u) {
						it_gt = it_le;
					}
					else {
						continue;
					}
				}
			}

			size_t a = (it_le - stroke.u.begin());
			size_t b = (it_gt - stroke.u.begin());

			double u_le = *it_le;
			double u_gt = *it_gt;
			if (u_gt < u_le) {
				if (u_le > u) {
					u_le -= period;
				}
				else {
					u_gt += period;
				}
			}

			double mix = 0;
			if (u_gt != u_le) {
				mix = (u - u_le) / (u_gt - u_le);
			}
			if (mix < 0 || mix > 1) {
				continue;
			}
			double i = a + mix;
			glm::dvec2 point = (1. - mix) * stroke.points[a] + mix * stroke.points[b];

			glm::dvec2 to_next(0.0, 0.0);
			if (it_gt != it_le) {
				to_next = stroke.points[b] - point;
			}
			else {
				to_next = stroke.points[b] - stroke.points[b - 1];
			}

			xsec.points.push_back({
				stroke_idx,
				i,
				point,
				to_next,
				tangent(stroke.points, i)
			});
		}
	}

	for (size_t a = 0; a < xsec.points.size(); ++a) {
		for (size_t b = a + 1; b < xsec.points.size(); ++b) {
			xsec.connections.push_back({ a, b, 1. });
		}
	}

	return xsec;
}

void Parameterization::ensure_connected(Cluster* cluster, Cluster::XSec* cut) {
	struct StrokeSegment {
		size_t stroke_idx;
		double begin;
		double end;
		bool small;
		bool cut;
	};
	std::unordered_map<size_t, std::vector<StrokeSegment>> segments;

	auto connection_crosses_cut = [&](const Cluster::XSec& xsec) {
		for (auto& conn : xsec.connections) {
			for (size_t idx : { conn.a_idx, conn.b_idx }) {
				size_t stroke = xsec.points[idx].stroke_idx;
				double i = xsec.points[idx].i;

				for (auto& pt : cut->points) {
					if (pt.stroke_idx == stroke && std::ceil(pt.i) >= std::ceil(i) && std::floor(pt.i) <= std::floor(i)) {
						return true;
					}
				}
			}
		}
		return false;
	};

	for (size_t stroke_idx = 0; stroke_idx < cluster->strokes.size(); ++stroke_idx) {
		std::vector<double> cuts;
		if (cut) {
			for (auto& pt : cut->points) {
				if (pt.stroke_idx == stroke_idx) {
					cuts.push_back(pt.i);
				}
			}
		}
		std::sort(cuts.begin(), cuts.end());
		StrokeSegment segment{ stroke_idx, -1., 0.0, false, !cuts.empty() };
		for (double cut_i : cuts) {
			segment.end = cut_i;
			segment.small = std::min<int>(cluster->strokes[stroke_idx].points.size() - 1, std::floor(segment.end)) -
				std::max<int>(0, std::ceil(segment.begin)) <= 2;
			//if (std::min<int>(cluster->strokes[stroke_idx].points.size() - 1, std::floor(segment.end)) -
			//	std::max<int>(0, std::ceil(segment.begin)) > 2) {
				segments[stroke_idx].push_back(segment);
			//}
			segment.begin = cut_i;
		}
		segment.end = cluster->strokes[stroke_idx].points.size();
		segment.small = std::min<int>(cluster->strokes[stroke_idx].points.size() - 1, std::floor(segment.end)) -
			std::max<int>(0, std::ceil(segment.begin)) <= 2;
		//if (std::min<int>(cluster->strokes[stroke_idx].points.size() - 1, std::floor(segment.end)) -
		//	std::max<int>(0, std::ceil(segment.begin)) > 2) {
			segments[stroke_idx].push_back(segment);
		//}
	}

	StrokeSegment* first_non_cut_segment = &segments[0][0];
	for (size_t stroke_idx = 0; stroke_idx < cluster->strokes.size(); ++stroke_idx) {
		if (segments[stroke_idx].size() == 1) {
			first_non_cut_segment = &segments[stroke_idx][0];
			break;
		}
	}

	bool ok = false;
	while (!ok) {
		std::unordered_map<StrokeSegment*, bool> connected;
		for (auto& kv : segments) {
			for (auto& seg : kv.second) {
				connected[&seg] = false;
			}
		}

		std::function<void(StrokeSegment*)> visit = [&](StrokeSegment* seg) -> void {
			if (connected[seg]) return;
			connected[seg] = true;

			for (auto& xsec : cluster->xsecs) {
				for (auto& c : xsec.connections) {
					if (c.weight < 1e-4) continue;
					bool connected_a = xsec.points[c.a_idx].stroke_idx == seg->stroke_idx && xsec.points[c.a_idx].i > seg->begin && xsec.points[c.a_idx].i < seg->end;
					bool connected_b = xsec.points[c.b_idx].stroke_idx == seg->stroke_idx && xsec.points[c.b_idx].i > seg->begin && xsec.points[c.b_idx].i < seg->end;

					if (connected_a || connected_b) {
						size_t other = connected_a ? c.b_idx : c.a_idx;
						size_t other_stroke = xsec.points[other].stroke_idx;
						double other_i = xsec.points[other].i;

						auto it = std::find_if(
							segments[other_stroke].begin(),
							segments[other_stroke].end(),
							[&](const StrokeSegment& s) { return s.begin <= other_i && s.end > other_i;  }
						);
						if (it != segments[other_stroke].end()) {
							auto& other_seg = *it;
							visit(&other_seg);
						}
						else {
							std::cout << "Can't find it" << std::endl;
						}
					}
				}
			}
		};
		visit(first_non_cut_segment);

		double closest_dist = std::numeric_limits<double>::infinity();
		Cluster::XSec connector = { {}, {{0, 1, 1.0}}, 0, 0.0, true };
		ok = true;
		for (auto& disconnected : connected) {
			if (disconnected.second) continue;
			if (disconnected.first->cut && disconnected.first->small) continue;
			ok = false;
			size_t stroke_a = disconnected.first->stroke_idx;
			double min_i = disconnected.first->begin;
			double max_i = disconnected.first->end;

			for (auto& kv : connected) {
				if (!kv.second) continue;
				if (kv.first->cut && kv.first->small) continue;
				size_t stroke_b = kv.first->stroke_idx;
				double min_j = kv.first->begin;
				double max_j = kv.first->end;

				for (int ia : { std::max<int>(0, std::ceil(min_i)), std::min<int>(cluster->strokes[stroke_a].points.size() - 1, std::floor(max_i)) }) {
					for (int ib : { std::max<int>(0, std::ceil(min_j)), std::min<int>(cluster->strokes[stroke_b].points.size() - 1, std::floor(max_j)) }) {

						double dist = glm::distance2(cluster->strokes[stroke_a].points[ia], cluster->strokes[stroke_b].points[ib]);
						if (dist < closest_dist) {
							glm::dvec2 to_next_a(0.0, 0.0);
							if (ia == 0) to_next_a = cluster->strokes[stroke_a].points[ia + 1] - cluster->strokes[stroke_a].points[ia];
							else to_next_a = cluster->strokes[stroke_a].points[ia] - cluster->strokes[stroke_a].points[ia - 1];

							glm::dvec2 to_next_b(0.0, 0.0);
							if (ib == 0) to_next_b = cluster->strokes[stroke_b].points[ib + 1] - cluster->strokes[stroke_b].points[ib];
							else to_next_b = cluster->strokes[stroke_b].points[ib] - cluster->strokes[stroke_b].points[ib - 1];

							auto new_connector = connector;
							new_connector.points = {
								{ stroke_a, double(ia), cluster->strokes[stroke_a].points[ia], to_next_a, tangent(cluster->strokes[stroke_a].points, double(ia)) },
								{ stroke_b, double(ib), cluster->strokes[stroke_b].points[ib], to_next_b, tangent(cluster->strokes[stroke_b].points, double(ib)) },
							};
							new_connector.u = 0.5 * (cluster->strokes[stroke_a].u[ia] + cluster->strokes[stroke_b].u[ib]);

							if (!cut || !connection_crosses_cut(new_connector)) {
								closest_dist = dist;
								connector = new_connector;
							}
						}
					}
				}
			}
		}

		if (!ok) {
			if (connector.points.size() == 0) {
				break;
			}
			cluster->xsecs.push_back(connector);
		}
	}
}

void Parameterization::check_periodic(Cluster* cluster) {
	std::vector<Cluster::XSec> ortho_xsecs = orthogonal_xsecs(*cluster, 0.05);

	struct Node;

	struct Edge {
		Node* to;
		int xsec_idx;
	};

	struct Node {
		size_t stroke_idx;
		int i;
		std::vector<Edge> edges;
	};

	std::vector<std::vector<Node>> nodes;
	nodes.reserve(cluster->strokes.size());
	for (size_t stroke_idx = 0; stroke_idx < cluster->strokes.size(); ++stroke_idx) {
		nodes.emplace_back(cluster->strokes[stroke_idx].u.size(), Node{ stroke_idx, 0, {} });
		for (int i = 0; i < nodes[stroke_idx].size(); ++i) {
			nodes[stroke_idx][i].i = i;

			// Connections due to monotonicity
			if (i < nodes[stroke_idx].size() - 1) {
				nodes[stroke_idx][i].edges.push_back({ &nodes[stroke_idx][i + 1], -1 });
			}
		}
	}
	// Known isolines
	for (int xsec_idx = 0; xsec_idx < cluster->xsecs.size(); xsec_idx += 3) {
		auto& xsec = cluster->xsecs[xsec_idx];
		for (auto& conn : xsec.connections) {
			size_t sa = xsec.points[conn.a_idx].stroke_idx;
			size_t sb = xsec.points[conn.b_idx].stroke_idx;
			double raw_ia = xsec.points[conn.a_idx].i;
			double raw_ib = xsec.points[conn.b_idx].i;
			int ia = std::round(raw_ia);
			int ib = std::round(raw_ib);
			nodes[sa][ia].edges.push_back({ &nodes[sb][ib], -1 });
			nodes[sb][ib].edges.push_back({ &nodes[sa][ia], -1 });
		}
	}
	// Orthogonal jumps
	for (int xsec_idx = 0; xsec_idx < ortho_xsecs.size(); ++xsec_idx) {
		auto& xsec = ortho_xsecs[xsec_idx];
		for (auto& conn : xsec.connections) {
			size_t sa = xsec.points[conn.a_idx].stroke_idx;
			size_t sb = xsec.points[conn.b_idx].stroke_idx;
			double raw_ia = xsec.points[conn.a_idx].i;
			double raw_ib = xsec.points[conn.b_idx].i;
			int ia = std::round(raw_ia);
			int ib = std::round(raw_ib);
			nodes[sa][ia].edges.push_back({ &nodes[sb][ib], xsec_idx });
			nodes[sb][ib].edges.push_back({ &nodes[sa][ia], xsec_idx });
		}
	}

	auto is_close = [&](const Node& a, const Cluster::XSec& b) {
		for (auto& p_b : b.points) {
			if (a.stroke_idx == p_b.stroke_idx &&
				std::abs(std::double_t(a.i) - p_b.i) < 5.0) {
				return true;
			}
		}
		return false;
	};

	auto find_path = [&](const Node& from, const Cluster::XSec& to) -> std::vector<int> {
		int allowed_jumps = 2;
		std::unordered_map<int, std::vector<std::vector<bool>>> visited;
		for (int jumped = 0; jumped <= allowed_jumps; ++jumped) {
			visited[jumped].reserve(nodes.size());
			for (auto& stroke_nodes : nodes) {
				visited[jumped].emplace_back(stroke_nodes.size(), false);
			}
		}

		std::function<std::vector<int>(const Node*, int)> visit = [&](const Node* n, int jumps) -> std::vector<int> {
			if (visited[jumps][n->stroke_idx][n->i]) return {};
			visited[jumps][n->stroke_idx][n->i] = true;

			for (auto& edge : n->edges) {
				if (jumps == 0 && edge.xsec_idx != -1) continue;

				if (is_close(*edge.to, to)) {
					if (context.debug_viz) {
						add_debug_line({
							cluster->strokes[n->stroke_idx].points[n->i],
							cluster->strokes[edge.to->stroke_idx].points[edge.to->i],
							edge.xsec_idx == -1 ? "#000" : "#F0F"
							});
					}
					return { edge.xsec_idx };
				}
				else {
					auto res = visit(edge.to, jumps - (edge.xsec_idx != -1 ? 1 : 0));
					if (!res.empty()) {
						if (context.debug_viz) {
							add_debug_line({
								cluster->strokes[n->stroke_idx].points[n->i],
								cluster->strokes[edge.to->stroke_idx].points[edge.to->i],
								edge.xsec_idx == -1 ? "#000" : "#F0F"
								});
						}
						res.push_back(edge.xsec_idx);
						return res;
					}
				}
			}

			return {};
		};
		
		return visit(&from, allowed_jumps);
	};

	auto pair = std::minmax_element(cluster->xsecs.begin(), cluster->xsecs.end(), [](const Cluster::XSec& a, const Cluster::XSec& b) { return a.u < b.u; });
	Cluster::XSec& from = *pair.second;
	Cluster::XSec& to = *pair.first;

	auto get_u = [&](const Cluster::XSecPoint& pt) {
		double mix = pt.i - std::floor(pt.i);
		return (1. - mix) * cluster->strokes[pt.stroke_idx].u[std::floor(pt.i)] + mix * cluster->strokes[pt.stroke_idx].u[std::ceil(pt.i)];
	};

	std::vector<double> period_samples;
	for (auto& stroke_nodes : nodes) {
		for (auto& node : stroke_nodes) {
			if (is_close(node, from)) {
				auto path = find_path(node, to);
				if (!path.empty()) {
					for (int xsec_idx : path) {
						if (xsec_idx == -1) continue;
						auto& step = ortho_xsecs[xsec_idx];
						for (auto& conn : step.connections) {
							double u_a = get_u(step.points[conn.a_idx]);
							double u_b = get_u(step.points[conn.b_idx]);
							double diff = std::abs(u_a - u_b);
							if (diff > 10.) {
								period_samples.push_back(diff);
							}
						}
					}
				}
			}
		}
	}

	// Not periodic
	if (period_samples.empty()) return;

	cluster->periodic = true;

	int off = 0.2 * (period_samples.size() - 1);
	std::nth_element(period_samples.begin(), period_samples.begin() + off, period_samples.end());
	double period = period_samples[off];

	// Make all u values be \in [0, period)
	for (auto& stroke : cluster->strokes) {
		for (auto& u : stroke.u) {
			while (u >= period) {
				u -= period;
			}
		}
	}
	cluster->xsecs = xsecs_from_params(*cluster);
}
