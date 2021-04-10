#pragma once

#include <map>

#include "Cluster.h"
#include "Context.h"

struct FittedCurve {
	std::vector<glm::dvec2> centerline;
	std::vector<double> widths;
};

struct Fitting {
	Fitting(const Context& context);
	std::map<int, FittedCurve> fit(Input* input);
	FittedCurve fit_cluster(Cluster cluster);

	void fit_svg(std::ostream& os, const Input& input, const std::map<int, FittedCurve>& fits);

private:
	const Context& context;

	const double K_WEIGHT = 1e3;
	const double POS_WEIGHT = 5e-5;
	const double MIN_DIST = 1e-4;

	struct Sample {
		glm::dvec2 point;
		glm::dvec2 tangent;
		double k;
		bool no_k;
		bool gap;
		double width;
	};
	std::vector<Sample> samples_from_xsec(const Cluster& cluster, size_t xsec_idx);

	std::vector<glm::dvec2> fit_tangents(const std::vector<Sample>& samples, bool periodic);
	std::vector<glm::dvec2> fit_positions(const std::vector<Sample>& samples, const std::vector<glm::dvec2>& tangents, bool periodic);
	std::vector<double> fit_widths(const std::vector<Sample>& samples, bool periodic);
};
