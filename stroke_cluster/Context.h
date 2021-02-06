#pragma once

#include <gurobi_c++.h>

struct Context {
	Context(bool grb_log = false);
	void optimize_model(GRBModel* model) const;

	GRBEnv grb;
	bool debug_viz = false;
};