#include "Parameterization.h"

#include <future>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <deque>
#include <unordered_set>
#include <glm/gtx/norm.hpp>

#include "Utils.h"
#include "SvgUtils.h"

Parameterization::Parameterization(bool viz) : viz(viz), grb(true) {
	//grb.set(GRB_IntParam_LogToConsole, 0);
	grb.start();
}

void Parameterization::parameterize(Input* input) {
	std::map<int, std::future<void>> futures;
	for (auto& kv : input->clusters) {
		futures[kv.first] = std::async(std::launch::async, [&]() -> void {
			return parameterize_cluster(&kv.second);
		});
	}
	for (auto& kv : input->clusters) {
		futures[kv.first].get();
	}
}

void Parameterization::isolines_svg(std::ostream& os, const Input& input) {
	input.cluster_svg(os, [&](std::ostream& os) {
		for (auto& kv : input.clusters) {
			for (auto& xsec : kv.second.xsecs) {
				for (auto& connection : xsec.connections) {
					if (std::abs(int(connection.a_idx) - int(connection.b_idx)) == 1) {
						std::stringstream ss;
						ss << "#";
						ss << std::setfill('0') << std::setw(2);
						if (connection.weight > 0.5) {
							ss << std::hex << int(map(connection.weight, 0.5, 1.0, 255.0, 0.0)); // red
							ss << std::hex << int(map(connection.weight, 0.5, 1.0, 255.0 / 2.0, 255.0)); // green
							ss << "00"; // blue
						}
						else {
							ss << "FF", // red
							ss << std::hex << int(map(connection.weight, 0.0, 0.5, 0.0, 255.0 / 2.0)); // green
							ss << "00"; // blue
						}
						SVG::line(
							os,
							xsec.points[connection.a_idx].point.x,
							xsec.points[connection.a_idx].point.y,
							xsec.points[connection.b_idx].point.x,
							xsec.points[connection.b_idx].point.y,
							0.5,
							ss.str()
						);
					}
				}
			}
		}
	});
}

void Parameterization::parameterize_cluster(Cluster* cluster) {
	cluster->xsecs = orthogonal_xsecs(*cluster);
	params_from_xsecs(cluster, true, nullptr);
}

GRBLinExpr l1_norm(GRBModel* model, const std::vector<GRBLinExpr>& x) {
	// min ||x||_1
	//
	// ...is equivalent to:
	//
	// min t
	// s.t.
	// x_i <= y_i,
	// -x_i <= y_i,
	// \sum_i y_i = t

	GRBLinExpr sum_y = 0.0;
	auto t = model->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
	std::vector<GRBVar> y;
	y.reserve(x.size());
	for (auto& term : x) {
		y.push_back(model->addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS));
		sum_y += y.back();
		model->addConstr(term <= y.back());
		model->addConstr(-term <= y.back());
	}
	model->addConstr(sum_y == t);

	return t;
}

GRBQuadExpr l2_norm_sq(GRBModel* model, const std::vector<GRBLinExpr>& x) {
	GRBQuadExpr result = 0.0;
	for (auto& term : x) {
		result += term * term;
	}

	return result;
}

void Parameterization::params_from_xsecs(Cluster* cluster, bool initial, Cluster::XSec* cut) {
	GRBModel model(grb);

	std::vector<std::vector<GRBVar>> param_vars;
	for (auto& stroke : cluster->strokes) {
		param_vars.emplace_back();
		auto& vars = param_vars.back();
		for (size_t i = 0; i < stroke.points.size(); ++i) {
			//std::string name = std::to_string(param_vars.size()-1) + "." + std::to_string(i);
			vars.push_back(model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS));
		}
	}

	std::vector<GRBLinExpr> velocity_terms;
	std::vector<GRBLinExpr> alignment_terms;
	for (auto& xsec : cluster->xsecs) {
		glm::dvec2 tangent = xsec.avg_tangent();

		// 1. Add velocity term
		{
			velocity_terms.emplace_back();
			auto& velocity_term = velocity_terms.back();

			double total_weight = 0.0;
			for (size_t i = 0; i < xsec.points.size(); ++i) {
				auto& point = xsec.points[i];
				double weight = xsec.distance_weight(i);
				total_weight += weight;
				double coefficient = weight * glm::dot(point.tangent, tangent) / glm::length(point.to_next);
				if (std::isnan(coefficient)) {
					throw "coefficient NaN";
				}

				int a = std::floor(point.i);
				int b = std::ceil(point.i);
				double mix = point.i - a;
				if (a == b) {
					if (b == cluster->strokes[point.stroke_idx].points.size() - 1) {
						continue;
					}
					else {
						++b;
					}
				}
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
		{
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
				alignment_terms.push_back(connection.weight * ((u_a - u_b) - proj_dist));
			}
		}
	}

	// 3. Enforce monotonicity
	for (size_t stroke_idx = 0; stroke_idx < cluster->strokes.size(); ++stroke_idx) {
		auto& stroke = cluster->strokes[stroke_idx];
		for (size_t i = 0; i < stroke.points.size() - 1; ++i) {
			model.addConstr(
				param_vars[stroke_idx][i + 1] - param_vars[stroke_idx][i] >=
					0.5 * glm::distance(stroke.points[i + 1], stroke.points[i]));
		}
	}

	// 4. Boundary
	//GRBQuadExpr boundary = param_vars[0][0] * param_vars[0][0];
	model.addConstr(param_vars[0][0] == 0.0);

	GRBQuadExpr objective;
	if (initial) {
		objective = l2_norm_sq(&model, velocity_terms) + 1e-5 * l1_norm(&model, alignment_terms);
	}
	else {
		objective = l2_norm_sq(&model, velocity_terms) + l2_norm_sq(&model, alignment_terms);
	}

	try {
		model.set(GRB_DoubleParam_FeasibilityTol, 1e-9);
		model.set(GRB_IntParam_DualReductions, 0);
		model.setObjective(objective, GRB_MINIMIZE);
		model.optimize();

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
	catch (GRBException e) {
		model.computeIIS();
		std::cout << "\nThe following constraint(s) "
			<< "cannot be satisfied:" << std::endl;
		GRBConstr* c = model.getConstrs();
		for (int i = 0; i < model.get(GRB_IntAttr_NumConstrs); ++i)
		{
			if (c[i].get(GRB_IntAttr_IISConstr) == 1)
			{
				std::cout << c[i].get(GRB_StringAttr_ConstrName) << std::endl;
			}
		}
		std::cout << "Done" << std::endl;

		model.write("D:\\model.lp");

		std::cout << "Error code = " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
		throw e;
	}
}

std::vector<Cluster::XSec> Parameterization::orthogonal_xsecs(const Cluster& cluster, double angle_tolerance) {
	const double TARGET_ANGLE = 20. / 180. * M_PI;
	const double MIN_GAP_CUTOFF = 15.;
	const double MAX_GAP_CUTOFF = 30.;

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

			if (neighbourhood.empty()) continue;

			// Sort by angle descending
			std::sort(neighbourhood.begin(), neighbourhood.end(), [](const AngleSample& a, const AngleSample& b) { return a.angle > b.angle; });

			// Remove samples until the `angle_tolerence` percentile is below the target value
			double max_angle = neighbourhood.front().angle;
			if (neighbourhood.back().angle <= TARGET_ANGLE) {
				while (!neighbourhood.empty() && neighbourhood[int(angle_tolerance * (neighbourhood.size() - 1))].angle > TARGET_ANGLE) {
					neighbourhood.pop_front();
					max_angle = neighbourhood.front().angle;
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
			double r_sq = std::max(100.0*100.0, glm::distance2(xsec.points.front().point, xsec.points.back().point));

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
	glm::dvec2 to_next;
	if (a < stroke_points.size() - 1) {
		to_next = (1 - mix) * (stroke_points[a + 1] - stroke_points[a]);
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

			glm::dvec2 to_next;
			size_t floor_j = std::floor(j);
			if (floor_j < cluster.strokes[j].points.size() - 1) {
				to_next = (1 - (j - double(floor_j))) * (cluster.strokes[j].points[floor_j + 1] - cluster.strokes[j].points[floor_j]);
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