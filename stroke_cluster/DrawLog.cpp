#include "DrawLog.h"
#include "Util.h"

#include <iostream>

void DrawLog::record_to(std::unique_ptr<DrawPrimitive> prim_ptr, FilterAttribute attribute,
						int layer_num) {
//#ifdef TO_PARALLEL
	return;
//#endif
	if (layers.count(layer_num) == 0) {
		layers[layer_num] = Layer(layer_num);
	}

	//if (layer_num == 30)
	//	std::cout << "### " << layer_num << " <- " << attribute.group_ind << std::endl;

	layers[layer_num].primitives.push_back(std::move(prim_ptr));
	layers[layer_num].attributes.emplace_back(attribute);
}

void DrawLog::set_layer_visibility(int layer_num, bool on) {
	if (layers.count(layer_num) == 0)
		return;
	layers[layer_num].to_draw = on;
	//std::cout << "Layer " << layer_num << " : " << layers[layer_num].to_draw << std::endl;
}

bool DrawLog::set_layer_filter(int layer_num, FilterAttribute filter, bool insert) {
	if (layers.count(layer_num) == 0)
		return false;

	if (insert) {
		if (layers[layer_num].filters.count(filter) == 0) {
			layers[layer_num].filters.emplace(filter);
			//std::cout << "Add filter " << layer_num << " : " << filter.group_ind << " - " << filter.stroke_ind << std::endl;
			return true;
		}
	}
	else {
		if (layers[layer_num].filters.count(filter) > 0) {
			layers[layer_num].filters.erase(filter);
			//std::cout << "Delete filter " << layer_num << " : " << filter.group_ind << " - " << filter.stroke_ind << std::endl;
			return true;
		}
	}

	return false;
}

void DrawLog::set_all_layers_filter(FilterAttribute filter, bool insert) {
	// Combine filters if the stroke index is -1
	// First, delete all individual filters with the input group index
	if (filter.stroke_ind == -1) {
		for (auto &l : layers) {
			for (auto itr = l.second.filters.begin(); itr != l.second.filters.end(); ) {
				bool to_delete = false;
				auto delete_itr = itr;

				if (delete_itr->group_ind == filter.group_ind)
					to_delete = true;
				itr++;

				if (to_delete)
					l.second.filters.erase(delete_itr);
			}
		}
	}

	// Add/delete as (group, -1)
	for (auto &l : layers) {
		if (insert)
			l.second.filters.emplace(filter);
		else {
			l.second.filters.erase(filter);
		}
	}
}

void DrawLog::reset_layer_filter(int layer_num) {
	if (layers.count(layer_num) == 0)
		return;
	layers[layer_num].filters.clear();
}

std::vector<int> DrawLog::get_layer_numbers() const {
	std::vector<int> numbers;
	for (auto &l : layers) {
		numbers.emplace_back(l.first);
	}

	return numbers;
}

void DrawLog::draw() const {
	for (auto const &l : layers) {
		if (l.second.to_draw) {
			assert(l.second.primitives.size() == l.second.attributes.size());
			for (size_t i = 0; i < l.second.primitives.size(); i++) {

				// Apply filter within the layer
				bool to_draw = false;
				for (auto const &f : l.second.filters) {
					if (f.match(l.second.attributes[i])) {
						to_draw = true;
						break;
					}
				}

				if (to_draw)
					l.second.primitives[i]->draw();
			}
		}
	}
}

std::string DrawLog::to_string() const {
	std::string str = "";

	for (auto const &l : layers) {
		if (l.second.to_draw) {
			assert(l.second.primitives.size() == l.second.attributes.size());
			for (size_t i = 0; i < l.second.primitives.size(); i++) {

				// Apply filter within the layer
				bool to_draw = false;
				for (auto const &f : l.second.filters) {
					if (f.match(l.second.attributes[i])) {
						to_draw = true;
						break;
					}
				}

				if (to_draw)
					str += l.second.primitives[i]->to_string();
			}
		}
	}

	return str;
}