#include "RegularizedTangentFitting.h"

#include <iostream>

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/LU>

double compute_edge_alignment(PointSample source, PointSample target) {
	glm::dvec2 s_tan = source.tangent;
	glm::dvec2 t_tan = target.tangent;

	glm::dvec2 e_dir = glm::normalize(target.position - source.position);

	double alignment = glm::dot(s_tan, e_dir);
	double inv_alignment = glm::dot(t_tan, e_dir);

	return std::max(alignment, inv_alignment);
}

std::vector<glm::dvec2> ls_solve_open(std::vector<glm::dvec2> const &samples, std::vector<glm::dvec2> const &tangents, 
									  std::vector<double> const &weights, std::vector<double> const &pos_weights, 
									  double gamma) {
	// Basically a first-order ODE with a regularization term instead of boundary conditions
	assert(tangents.size() >= 1);
	assert(tangents.size() == samples.size() - 1);

	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2 * samples.size() - 1, samples.size());
	std::vector<Eigen::VectorXd> b;
	b.resize(tangents.front().length(), Eigen::VectorXd(2 * samples.size() - 1));

	// Regularization with the curvature term
	// Computed by the 1D Laplacian
	Eigen::MatrixXd L = Eigen::MatrixXd::Zero(samples.size() - 2, samples.size());

	double endpoint_weight = std::sqrt(2.);
	double sqrt_gamma = std::sqrt(gamma);
	std::vector<double> W_diag;
	for (size_t j = 0; j + 1 < samples.size(); j++) {
		W_diag.emplace_back(glm::distance(samples[j], samples[j + 1]));
	}

	// A
	{
		// Forward difference matrix D
		for (size_t i = 0; i < tangents.size(); i++) {
			A.coeffRef(i, i) = -1 / W_diag[i];
			A.coeffRef(i, i + 1) = 1 / W_diag[i];
		}

		// Position term
		for (size_t i = 0; i < samples.size(); i++) {
			A.coeffRef(tangents.size() + i, i) = sqrt_gamma * pos_weights[i];

			if (i == 0 || i == samples.size() - 1)
				A.coeffRef(tangents.size() + i, i) *= endpoint_weight;
		}
	}

	// b
	{
		// tangent
		for (size_t i = 0; i < tangents.size(); i++) {
			for (size_t j = 0; j < tangents.front().length(); j++) {
				b[j][i] = tangents[i][j];
			}
		}

		// Initial sample positions
		for (size_t i = 0; i < samples.size(); i++) {
			for (size_t j = 0; j < tangents.front().length(); j++) {
				b[j][tangents.size() + i] = sqrt_gamma * pos_weights[i] * samples[i][j];

				if (i == 0 || i == samples.size() - 1)
					b[j][tangents.size() + i] *= endpoint_weight;
			}
		}
	}

	// L
	{
		// Laplacian matrix
		// TODO: weights?
		for (size_t i = 0; i + 2 < samples.size(); i++) {
			L.coeffRef(i, i) = 1 * weights[i + 1];
			L.coeffRef(i, i + 1) = -2 * weights[i + 1];
			L.coeffRef(i, i + 2) = 1 * weights[i + 1];
		}
	}

	// Solve for the normal equation
	// Hmm... seems that if the polyline size is fix and we move the weighting to the RHS,
	// we can LU the same difference matrix only once
	// Though speed is not a problem (at least for now)

	//Eigen::PartialPivLU<Eigen::MatrixXd> solver;
	Eigen::LDLT<Eigen::MatrixXd> solver;
	solver.compute(A.transpose() * A + L.transpose() * L);

	std::vector<Eigen::VectorXd> x;
	x.resize(tangents.front().length());
	for (size_t j = 0; j < tangents.front().length(); j++) {
		x[j] = solver.solve(A.transpose() * b[j]);
	}

	std::vector<glm::dvec2> results;
	results.resize(samples.size());

	for (size_t i = 0; i < samples.size(); i++) {
		for (size_t j = 0; j < tangents.front().length(); j++) {
			results[i][j] = x[j][i];
		}
	}

	return results;
}

std::vector<glm::dvec2> ls_solve_closed(std::vector<glm::dvec2> const &samples, std::vector<glm::dvec2> const &tangents, 
										std::vector<double> const &weights, std::vector<double> const &pos_weights, 
										double gamma) {
	// Basically a first-order ODE with a regularization term instead of boundary conditions
	assert(tangents.size() >= 1);
	assert(tangents.size() == samples.size());

	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2 * samples.size(), samples.size());
	std::vector<Eigen::VectorXd> b;
	b.resize(tangents.front().length(), Eigen::VectorXd(2 * samples.size()));

	// Regularization with the curvature term
	// Computed by the 1D Laplacian
	Eigen::MatrixXd L = Eigen::MatrixXd::Zero(samples.size(), samples.size());

	double endpoint_weight = std::sqrt(2.);
	double sqrt_gamma = std::sqrt(gamma);
	std::vector<double> W_diag;
	for (size_t j = 0; j < samples.size(); j++) {
		if (j + 1 < samples.size())
			W_diag.emplace_back(glm::distance(samples[j], samples[j + 1]));
		else
			W_diag.emplace_back(glm::distance(samples[j], samples[0]));
	}

	// A
	{
		// Circulant difference matrix D
		// divided by W
		for (size_t i = 0; i < samples.size(); i++) {
			if (i + 1 != samples.size()) {
				A.coeffRef(i, i) = -1 / W_diag[i];
				A.coeffRef(i, i + 1) = 1 / W_diag[i];
			}
			else {
				A.coeffRef(i, i) = -1 / W_diag[i];
				A.coeffRef(i, 0) = 1 / W_diag[i];
			}
		}

		// Position term
		for (size_t i = 0; i < samples.size(); i++) {
			A.coeffRef(samples.size() + i, i) = sqrt_gamma * pos_weights[i];

			if (i == 0 || i == samples.size() - 1)
				A.coeffRef(samples.size() + i, i) *= endpoint_weight;
		}
	}

	// b
	{
		// tangent
		for (size_t i = 0; i < tangents.size(); i++) {
			for (size_t j = 0; j < tangents.front().length(); j++) {
				b[j][i] = tangents[i][j];
			}
		}

		// Initial sample positions
		for (size_t i = 0; i < samples.size(); i++) {
			for (size_t j = 0; j < tangents.front().length(); j++) {
				b[j][tangents.size() + i] = sqrt_gamma * pos_weights[i] * samples[i][j];

				if (i == 0 || i == samples.size() - 1)
					b[j][tangents.size() + i] *= endpoint_weight;
			}
		}
	}

	// L
	{
		// Laplacian matrix
		// TODO: weights?
		for (size_t i = 0; i < samples.size(); i++) {
			size_t prev = (i + samples.size() - 1) % samples.size();
			size_t next = (i + 1) % samples.size();

			L.coeffRef(i, i) = -2 * weights[i];
			L.coeffRef(i, next) = 1 * weights[i];
			L.coeffRef(i, prev) = 1 * weights[i];
		}
	}

	// Solve for the normal equation
	// Hmm... seems that if the polyline size is fix and we move the weighting to the RHS,
	// we can LU the same difference matrix only once
	// Though speed is not a problem (at least for now)
	Eigen::LDLT<Eigen::MatrixXd> solver;
	solver.compute(A.transpose() * A + L.transpose() * L);

	std::vector<Eigen::VectorXd> x;
	x.resize(tangents.front().length());
	for (size_t j = 0; j < tangents.front().length(); j++) {
		x[j] = solver.solve(A.transpose() * b[j]);
	}

	std::vector<glm::dvec2> results;
	results.resize(samples.size());

	for (size_t i = 0; i < samples.size(); i++) {
		for (size_t j = 0; j < tangents.front().length(); j++) {
			results[i][j] = x[j][i];
		}
	}

	return results;
}

// samples: n; tangents: n - 1 (on edge), weighted based on the length of the corresponding segment
std::vector<glm::dvec2> ls_solve(std::vector<glm::dvec2> const &samples, std::vector<glm::dvec2> const &tangents, 
								 std::vector<double> const &weights, std::vector<double> const &pos_weights, 
								 double gamma, bool is_closed) {
	assert(tangents.size() >= 1);

	// Let f: [0, L] -> R^2, [0, L] in R, L is the length of the spline
	// Using [0, L] instead of [0, 1] to normalize the tangent

	if (!is_closed)
		return ls_solve_open(samples, tangents, weights, pos_weights, gamma);
	else
		return ls_solve_closed(samples, tangents, weights, pos_weights, gamma);
}