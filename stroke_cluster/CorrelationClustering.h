#pragma once

#include <iostream>
#include <vector>
#include <map>

#include "andres/graph/graph.hxx"
#include "SketchInfo.h"

#define EDGE_WEIGHTER OverlappingOnlyEdgeWeighter
//#define EDGE_WEIGHTER OverlappingDropoffEdgeWeighter
//#define EDGE_WEIGHTER OverlappingEdgeWeighter

class EdgeDeterminer {
public:
	bool edge_exists(Sketch const &s1, Sketch const &s2) const;
};

class EdgeWeighter {
public:
	virtual void pre_compute(Sketch const &s1, Sketch const &s2) { }
	virtual bool is_similar(std::vector<Sketch> const &strokes, size_t ind1, size_t ind2) const;
	virtual double weight_edge(Sketch const &s1, Sketch const &s2, bool is_similar = true) const;
};

class OverlappingEdgeWeighter : public EdgeWeighter {
public:
	OverlappingEdgeWeighter() : length_short(-1), length_tiny(-1) { 
		tangent_diff.resize(2, 0);
	}
	OverlappingEdgeWeighter(double _short, double _tiny) : length_short(_short), length_tiny(_tiny) { 
		tangent_diff.resize(2, 0);
	}

	void pre_compute(Sketch const &s1, Sketch const &s2) override;
	bool is_similar(std::vector<Sketch> const &strokes, size_t ind1, size_t ind2) const override;
	double weight_edge(Sketch const &s1, Sketch const &s2, bool is_similar = true) const override;

	std::vector<double> tangent_diff;

	double length_short;
	double length_tiny;
};

class OverlappingDropoffEdgeWeighter : public EdgeWeighter {
public:
	OverlappingDropoffEdgeWeighter() : length_short(-1), length_tiny(-1) { 
		tangent_diff.resize(2, 0);
	}
	OverlappingDropoffEdgeWeighter(double _short, double _tiny) : length_short(_short), length_tiny(_tiny) { 
		tangent_diff.resize(2, 0);
	}

	void pre_compute(Sketch const &s1, Sketch const &s2) override;
	bool is_similar(std::vector<Sketch> const &strokes, size_t ind1, size_t ind2) const override;
	double weight_edge(Sketch const &s1, Sketch const &s2, bool is_similar = true) const override;

	std::vector<double> tangent_diff;

	double length_short;
	double length_tiny;
};

class OverlappingOnlyEdgeWeighter : public EdgeWeighter {
public:
	OverlappingOnlyEdgeWeighter() : length_short(-1), length_tiny(-1) { 
		tangent_diff.resize(3, 0);
	}
	OverlappingOnlyEdgeWeighter(double _short, double _tiny) : length_short(_short), length_tiny(_tiny) { 
		tangent_diff.resize(3, 0);
	}

	void pre_compute(Sketch const &s1, Sketch const &s2) override;
	bool is_similar(std::vector<Sketch> const &strokes, size_t ind1, size_t ind2) const override;
	double weight_edge(Sketch const &s1, Sketch const &s2, bool is_similar = true) const override;

	std::vector<double> tangent_diff;
	Sketch merged;

	double length_short;
	double length_tiny;
};


class CorrelationClustering {
public:
	CorrelationClustering() {}
	~CorrelationClustering() {}

	static std::set<std::string> pairnames;

	void init(std::vector<Sketch> const &strokes, EdgeDeterminer const &determiner, EDGE_WEIGHTER &weighter);
	void cluster(std::vector<uint16_t> &vertex_labels) const;
	
	void clear();

	size_t get_num_of_edges() const {
		return stroke_graph.numberOfEdges();
	}

	double get_weight(int i, int j) const;

private:
	andres::graph::Graph<> stroke_graph;
	std::vector<double> weights;
};
