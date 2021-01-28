#pragma once

#include <vector>

#include <glm/glm.hpp>

#include "andres/graph/graph.hxx"

#include "SketchInfo.h"

std::vector<glm::dvec2> build_proximity_graph(std::vector<Sketch> const &strokes, double threshold, 
						   andres::graph::Graph<> &sample_graph, std::vector<size_t> &labels, std::vector<double> &weights);

bool order_strokes(std::vector<Sketch> const &strokes, double epsilon, std::vector<int> &directions, double step_size = -1);

void order_strokes_simple(std::vector<Sketch> const &strokes, std::vector<int> &directions, size_t end_offset = 0);

bool order_strokes_h0(std::vector<Sketch> const &strokes, double epsilon, std::vector<int> &directions);