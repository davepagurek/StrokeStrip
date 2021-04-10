#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <vector>
#include <future>
#include <mutex>
#include <deque>
#include <glm/glm.hpp>
#include <gurobi_c++.h>
#include "Cluster.h"

struct Intersection {
	double i;
	glm::dvec2 pt;
};
std::vector<Intersection> intersections(const std::vector<glm::dvec2>& polyline, glm::dvec2 from, glm::dvec2 to);

template<typename T> T map(T in, T in_from, T in_to, T out_from, T out_to) {
	T mixed = (in - in_from) / (in_to - in_from);
	T res = mixed * (out_to - out_from) + out_from;
	if (res < std::min(out_from, out_to)) res = std::min(out_from, out_to);
	if (res > std::max(out_from, out_to)) res = std::max(out_from, out_to);
	return res;
}

double gaussian(double x, double sigma, double mu);

double poly_in_out(double t, double k);

GRBLinExpr l1_norm(GRBModel* model, const std::vector<GRBLinExpr>& x);
GRBQuadExpr l2_norm_sq(GRBModel* model, const std::vector<GRBLinExpr>& x);

template<typename RetType>
std::map<int, RetType> map_clusters(Input& input, std::function<RetType(Cluster&)> process) {
	const int NUM_TASKS = 4;
	if (NUM_TASKS == 1) {
		std::map<int, RetType> vals;
		for (auto& kv : input.clusters) {
			vals[kv.first] = process(kv.second);
		}
		return vals;
	}
	else {
		std::deque<int> keys;
		std::vector<std::future<void>> futures;
		std::map<int, RetType> vals;
		std::mutex queue_lock;
		std::mutex vals_lock;
		for (auto& kv : input.clusters) {
			keys.push_back(kv.first);
			vals[kv.first] = {};
		}

		for (size_t task = 0; task < NUM_TASKS; ++task) {
			futures.push_back(std::async(std::launch::async, [&]() -> void {
				while (true) {
					int id = -1;
					{
						std::lock_guard<std::mutex> lock(queue_lock);
						if (keys.empty()) break;
						id = keys.front();
						keys.pop_front();
					}

					auto val = process(input.clusters[id]);
					{
						std::lock_guard<std::mutex> lock(vals_lock);
						vals[id] = val;
					}
				}
				}));
		}

		for (auto& future : futures) {
			future.get();
		}

		return vals;
	}
}