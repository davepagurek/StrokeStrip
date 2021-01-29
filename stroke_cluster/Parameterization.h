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

	std::vector<Cluster::XSec> orthogonal_xsecs(const Cluster& cluster);
	Cluster::XSec orthogonal_xsec_at(const Cluster& cluster, size_t stroke, double i);
};