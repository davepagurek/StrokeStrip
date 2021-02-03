#pragma once
#include <gurobi_c++.h>
#include <ostream>

#include "Cluster.h"

class Parameterization {
public:
	Parameterization(bool viz);
	void parameterize(Input* input);
	void isolines_svg(std::ostream& os, const Input& input);

private:
	bool viz;

	GRBEnv grb;
	void parameterize_cluster(Cluster* cluster);
	void params_from_xsecs(Cluster* cluster, bool initial = false, Cluster::XSec* cut = nullptr);
	void ensure_connected(Cluster* cluster);

	std::vector<Cluster::XSec> orthogonal_xsecs(const Cluster& cluster, double angle_tolerance = 0.0);
	std::vector<Cluster::XSec> xsecs_from_params(const Cluster& cluster, bool nonlinear = false);
	Cluster::XSec xsec_at_u(const Cluster& cluster, double u);
	Cluster::XSec orthogonal_xsec_at(const Cluster& cluster, size_t stroke, double i);
};