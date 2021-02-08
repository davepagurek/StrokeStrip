#pragma once
#define _USE_MATH_DEFINES
#include <ostream>
#include <mutex>

#include "Cluster.h"
#include "Context.h"

class Parameterization {
public:
	Parameterization(const Context& context);
	void parameterize(Input* input);
	void isolines_svg(std::ostream& os, const Input& input);
	void debug_svg(std::ostream& os, const Input& input);

	std::vector<Cluster::XSec> xsecs_from_params(const Cluster& cluster, double sample_rate = 1., bool nonlinear = false);

private:
	bool viz;

	const Context& context;
	std::mutex viz_lock;

	struct DebugLine {
		glm::dvec2 from;
		glm::dvec2 to;
		std::string color;
	};
	void add_debug_line(DebugLine line);
	std::vector<DebugLine> debug_lines;

	void parameterize_cluster(Cluster* cluster);
	void params_from_xsecs(Cluster* cluster, bool initial = false, Cluster::XSec* cut = nullptr);
	void ensure_connected(Cluster* cluster, Cluster::XSec* cut = nullptr);
	void check_periodic(Cluster* cluster);

	std::vector<Cluster::XSec> orthogonal_xsecs(const Cluster& cluster, double angle_tolerance = 0.0);
	Cluster::XSec xsec_at_u(const Cluster& cluster, double u);
	Cluster::XSec orthogonal_xsec_at(const Cluster& cluster, size_t stroke, double i);
};