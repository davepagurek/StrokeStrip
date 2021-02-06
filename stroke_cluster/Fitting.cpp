#include <future>

#include "Fitting.h"
#include "Parameterization.h"
#include "Utils.h"
#include "SvgUtils.h"

Fitting::Fitting(const Context& context) : context(context) {}

void Fitting::fit_svg(std::ostream& os, const Input& input, const std::map<int, std::vector<glm::dvec2>>& fits) {
	double padding = input.thickness;
	double w = input.width * input.thickness + 2 * padding;
	double h = input.height * input.thickness + 2 * padding;
	double x = -w / 2;
	double y = -h / 2;

	SVG::begin(os, x, y, w, h);
	for (auto& kv : fits) {
		std::vector<glm::dvec2> path = kv.second;
		for (auto& pt : path) {
			pt *= input.thickness;
		}
		SVG::polyline(os, path, input.thickness);
	}
	SVG::end(os);
}

std::map<int, std::vector<glm::dvec2>> Fitting::fit(const Input& input) {
	std::map<int, std::vector<glm::dvec2>> fits;
	std::map<int, std::future<std::vector<glm::dvec2>>> futures;
	for (auto& kv : input.clusters) {
		futures[kv.first] = std::async(std::launch::async, [&]() -> std::vector<glm::dvec2> {
			return fit_cluster(kv.second);
			});
	}
	for (auto& kv : input.clusters) {
		fits[kv.first] = futures[kv.first].get();
	}
	return fits;
	/*std::map<int, std::vector<glm::dvec2>> fits;
	for (auto& kv : input.clusters) {
		fits[kv.first] = fit_cluster(kv.second);
	}
	return fits;*/
}

std::vector<glm::dvec2> Fitting::fit_cluster(Cluster cluster) {
	// 1. Get xsec samples more densely sampled near endpoints
	{
		Parameterization parameterization(context);
		double rate = cluster.periodic ? 0.025 : 0.1;
		cluster.xsecs = parameterization.xsecs_from_params(cluster, rate, !cluster.periodic);
	}

	// 2. Make one fitting sample per xsec
	std::vector<Fitting::Sample> samples;
	for (size_t i = 0; i < cluster.xsecs.size(); ++i) {
		auto res = samples_from_xsec(cluster, i);
		samples.insert(samples.end(), res.begin(), res.end());
	}

	// 3. Sovle for tangents
	auto tangents = fit_tangents(samples, cluster.periodic);

	// 4. Solve for positions
	return fit_positions(samples, tangents, cluster.periodic);
}

std::vector<glm::dvec2> Fitting::fit_positions(const std::vector<Sample>& samples, const std::vector<glm::dvec2>& tangents, bool periodic) {
	std::vector<glm::dvec2> positions;

	GRBModel model(context.grb);
	struct VecVar {
		GRBVar x;
		GRBVar y;
	};
	std::vector<VecVar> vars;
	vars.reserve(samples.size());
	for (size_t i = 0; i < samples.size(); ++i) {
		GRBVar x = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
		GRBVar y = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
		vars.push_back({ x, y });
	}

	std::vector<GRBLinExpr> pos_terms;
	std::vector<GRBLinExpr> tan_terms;
	for (int i = 0; i < samples.size(); ++i) {
		pos_terms.emplace_back(vars[i].x - samples[i].point.x);
		pos_terms.emplace_back(vars[i].y - samples[i].point.y);

		if (periodic || i < samples.size() - 1) {
			int j = (i + 1) % samples.size();

			double target_len = std::max(
				MIN_DIST,
				glm::dot(samples[j].point - samples[i].point, tangents[i]));
			tan_terms.emplace_back((vars[j].x - vars[i].x) / target_len - tangents[i].x);
			tan_terms.emplace_back((vars[j].y - vars[i].y) / target_len - tangents[i].y);
		}
	}

	model.setObjective(l2_norm_sq(&model, tan_terms) + POS_WEIGHT * l2_norm_sq(&model, pos_terms));

	context.optimize_model(&model);

	positions.reserve(vars.size());
	for (auto& vec : vars) {
		positions.emplace_back(vec.x.get(GRB_DoubleAttr_X), vec.y.get(GRB_DoubleAttr_X));
	}
	if (periodic) {
		positions.push_back(positions.front());
	}

	return positions;
}

std::vector<glm::dvec2> Fitting::fit_tangents(const std::vector<Sample>& samples, bool periodic) {
	std::vector<glm::dvec2> tangents;

	GRBModel model(context.grb);
	struct VecVar {
		GRBVar x;
		GRBVar y;
	};
	std::vector<VecVar> vars;
	vars.reserve(samples.size());
	for (size_t i = 0; i < samples.size(); ++i) {
		GRBVar x = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
		GRBVar y = model.addVar(-GRB_INFINITY, GRB_INFINITY, 1.0, GRB_CONTINUOUS);
		vars.push_back({ x, y });
	}

	std::vector<GRBLinExpr> tan_terms;
	std::vector<GRBLinExpr> k_terms;
	for (int i = 0; i < samples.size(); ++i) {
		tan_terms.emplace_back(vars[i].x - samples[i].tangent.x);
		tan_terms.emplace_back(vars[i].y - samples[i].tangent.y);

		if (periodic || i < samples.size() - 1) {
			int j = (i + 1) % samples.size();
			glm::dvec2 avg_tangent = glm::normalize(samples[i].tangent + samples[j].tangent);
			glm::dvec2 ortho = normal(avg_tangent);

			double scale = samples[i].no_k ? 0.1 : 1.;

			k_terms.emplace_back(scale * (vars[j].x - vars[i].x - samples[i].k * ortho.x));
			k_terms.emplace_back(scale * (vars[j].y - vars[i].y - samples[i].k * ortho.y));
		}
	}

	model.setObjective(l2_norm_sq(&model, tan_terms) + K_WEIGHT * l2_norm_sq(&model, k_terms));
	
	context.optimize_model(&model);

	tangents.reserve(vars.size());
	for (auto& vec : vars) {
		tangents.push_back(glm::normalize(glm::dvec2{ vec.x.get(GRB_DoubleAttr_X), vec.y.get(GRB_DoubleAttr_X) }));
	}

	return tangents;
}

std::vector<Fitting::Sample> Fitting::samples_from_xsec(const Cluster& cluster, size_t xsec_idx) {
	std::vector<Fitting::Sample> result;
	auto& xsec = cluster.xsecs[xsec_idx];

	if (xsec.connector && glm::distance(xsec.points[0].point, xsec.points[1].point) > 1e-4) {
		glm::dvec2 pt1 = xsec.points[0].point;
		glm::dvec2 pt2 = xsec.points[1].point;
		glm::dvec2 tan1 = xsec.points[0].tangent;
		glm::dvec2 tan2 = xsec.points[1].tangent;

		glm::dvec2 tan = xsec.avg_tangent();
		glm::dvec2 ortho = normal(tan);

		// Find out which point leads into the next
		size_t stroke1 = xsec.points[0].stroke_idx;
		if (xsec_idx > 0) {
			auto& prev = cluster.xsecs[xsec_idx - 1];
			if (!std::any_of(prev.points.begin(), prev.points.end(), [=](const Cluster::XSecPoint& p) { return p.stroke_idx == stroke1; })) {
				std::swap(pt1, pt2);
				std::swap(tan1, tan2);
			}
		}
		else {
			auto& next = cluster.xsecs[xsec_idx + 1];
			if (std::any_of(next.points.begin(), next.points.end(), [=](const Cluster::XSecPoint& p) { return p.stroke_idx == stroke1; })) {
				std::swap(pt1, pt2);
				std::swap(tan1, tan2);
			}
		}

		double dist = glm::distance(pt1, pt2);
		size_t num_samples = std::max(1., std::ceil(dist / 0.05));
		double sign = glm::dot(tan2 - tan1, ortho) > 0 ? 1 : -1;
		double k = glm::length(tan2 - tan1) / double(num_samples) * sign;
		for (size_t i = 0; i < num_samples; ++i) {
			double t = double(i) / double(num_samples - 1);
			result.push_back({
				(1. - t) * pt1 + t * pt2,
				tan,
				k,
				false,
				true
			});
		}
	}
	else {
		result.push_back({ xsec.avg_point(), xsec.avg_tangent(), 0., false, false });
		auto& sample = result.back();

		double change_mags = 0.;
		int num_changes = 0;
		for (auto& p : xsec.points) {
			double di = 0.1;

			// Don't add tangent change at endpoints
			if (p.i < di || p.i >= cluster.strokes[p.stroke_idx].points.size() - 1 - di) continue;

			glm::dvec2 ortho = normal(p.tangent);

			glm::dvec2 prev_pt(point(cluster.strokes[p.stroke_idx].points, p.i - di));
			glm::dvec2 pt(p.point);
			glm::dvec2 next_pt(point(cluster.strokes[p.stroke_idx].points, p.i + di));

			glm::dvec2 prev_tan(tangent(cluster.strokes[p.stroke_idx].points, p.i - di));
			glm::dvec2 tan(p.tangent);
			glm::dvec2 next_tan(tangent(cluster.strokes[p.stroke_idx].points, p.i + di));

			double k1 = glm::dot(ortho, tan - prev_tan) / std::max(1e-4, glm::distance(pt, prev_pt));
			double k2 = glm::dot(ortho, next_tan - tan) / std::max(1e-4, glm::distance(pt, next_pt));
			change_mags += (k1 + k2) / 2.;
			++num_changes;
		}

		if (num_changes > 0) {
			change_mags /= double(num_changes);
		}
		else {
			sample.no_k = true;
		}
		double dist;
		if (xsec_idx < cluster.xsecs.size() - 1) {
			dist = glm::dot(normal(sample.tangent), cluster.xsecs[xsec_idx + 1].avg_point() - sample.point);
		}
		else {
			dist = glm::dot(normal(sample.tangent), sample.point - cluster.xsecs[xsec_idx - 1].avg_point());
		}
		sample.k = std::max(MIN_DIST, dist) * change_mags;
	}

	return result;
}