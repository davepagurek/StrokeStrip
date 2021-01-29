#include "Parameterization.h"

#include <future>
#include <sstream>
#include <iomanip>
#include <cmath>

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

std::vector<Cluster::XSec> Parameterization::orthogonal_xsecs(const Cluster& cluster) {
	std::vector<Cluster::XSec> xsecs;

	for (size_t stroke = 0; stroke < cluster.strokes.size(); ++stroke) {
		for (size_t i = 0; i < cluster.strokes[stroke].points.size(); ++i) {
			xsecs.push_back(orthogonal_xsec_at(cluster, stroke, i));
		}
	}

	return xsecs;
}

Cluster::XSec Parameterization::orthogonal_xsec_at(const Cluster& cluster, size_t stroke, double i) {
	const double SCALE = 50.;

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
					tangent(cluster.strokes[j].points, j)
				},
				glm::dot(intersection.pt - origin, ortho)
			});
		}
	}

	std::sort(potential_ints.begin(), potential_ints.end(), [](const PotentialInt& a, const PotentialInt& b) { return a.signed_dist < b.signed_dist; });

	if (!potential_ints.size() > 1) {
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

		potential_ints = std::vector<PotentialInt>(potential_ints.begin() + (1 + first_bad_tangent), potential_ints.begin() + (last_bad_tangent - 1));
	}

	Cluster::XSec xsec = { {}, {}, 0.0, false };
	for (auto& potential_int : potential_ints) {
		xsec.points.push_back(potential_int.point);
	}
	for (size_t a = 0; a < xsec.points.size(); ++a) {
		for (size_t b = a + 1; b < xsec.points.size(); ++b) {
			xsec.connections.push_back({ a, b, 1.0 });
		}
	}

	return xsec;
}