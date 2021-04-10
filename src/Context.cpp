#define _USE_MATH_DEFINES
#include <cmath>

#include "Context.h"

Context::Context(bool grb_log) : grb(true) {
	if (!grb_log) {
		grb.set(GRB_IntParam_LogToConsole, 0);
	}
	grb.start();
}

void Context::optimize_model(GRBModel* model) const {
	try {
		//model.set(GRB_DoubleParam_FeasibilityTol, 1e-2);
		model->optimize();
	}
	catch (GRBException e) {
		try {
			//model.set(GRB_IntParam_DualReductions, 0);
			model->set(GRB_IntParam_BarHomogeneous, 1);
			model->set(GRB_DoubleParam_PSDTol, 1e-4);
			//model->set(GRB_IntParam_NonConvex, 2);
			model->optimize();
		}
		catch (GRBException e2) {
			std::cout << "Error code = " << e2.getErrorCode() << std::endl;
			std::cout << e2.getMessage() << std::endl;
			throw e2;
		}
	}
	if (model->get(GRB_IntAttr_Status) == GRB_NUMERIC) {
		//model.set(GRB_IntParam_DualReductions, 0);
		model->set(GRB_IntParam_BarHomogeneous, 1);
		model->optimize();
	}
	if (model->get(GRB_IntAttr_Status) == GRB_NUMERIC) {
		//model->set(GRB_IntParam_BarHomogeneous, -1);
		model->set(GRB_IntParam_ScaleFlag, 2);
		model->set(GRB_DoubleParam_ObjScale, -0.5);
		model->optimize();
	}
}