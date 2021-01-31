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
	size_t b = std::ceil(i);
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
			if (i == j && std::abs(intersection.i - i) > 1e-1) continue;

			glm::dvec2 to_next;
			size_t floor_j = std::floor(j);
			if (floor_j < cluster.strokes[j].points.size() - 1) {
				to_next = (j - double(floor_j)) * (cluster.strokes[j].points[floor_j + 1] - cluster.strokes[j].points[floor_j]);
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