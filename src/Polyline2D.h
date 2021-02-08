#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <vector>
#include <cstdint>
#include <sstream>
#include <regex>
#include <assert.h>
#include <glm/glm.hpp>

extern double epsilon_small;

namespace SketchUI
{
	inline std::vector<std::string> split_str(std::string str, std::string pattern) {
		std::vector<std::string> strs;
		
		size_t pos = 0;
		std::string token;
		while ((pos = str.find(pattern)) != std::string::npos) {
			token = str.substr(0, pos);
			strs.emplace_back(token);
			str.erase(0, pos + pattern.length());
		}
		strs.emplace_back(str);

		return strs;
	}

	class Point2D
	{
	public: // data
		double x;
		double y;

		// the target tangent
		// if none, these two value should be zeros
		double tangent_x;
		double tangent_y;

		// curvature only used for preprocessing
		double curvature;

	public: // constructor
		Point2D(void) : tangent_x(0), tangent_y(0), curvature(0) {}
		Point2D(double p_x, double p_y) : x(p_x), y(p_y), tangent_x(0), tangent_y(0), curvature(0) {}
		Point2D(double p_x, double p_y, double t_x, double t_y) : x(p_x), y(p_y), tangent_x(t_x), tangent_y(t_y), curvature(0) {}
		Point2D(double p_x, double p_y, double t_x, double t_y, double curv) : x(p_x), y(p_y), tangent_x(t_x), tangent_y(t_y), curvature(curv) {}

	public: // operations
		const Point2D operator-(const Point2D& right) const
		{
			return Point2D(x - right.x, y - right.y);
		}

		const Point2D operator+(const Point2D& right) const
		{
			return Point2D(x + right.x, y + right.y);
		}
		Point2D& operator+=(const Point2D& right)
		{
			x += right.x;
			y += right.y;
			return *this;
		}

		Point2D& operator/=(double divider)
		{
			x /= divider;
			y /= divider;
			return *this;
		}

		Point2D& operator=(const Point2D& right)
		{
			x = right.x;
			y = right.y;
			tangent_x = right.tangent_x;
			tangent_y = right.tangent_y;

			return *this;
		}

		double SqLength(void) const { return x*x+y*y; }
		double Length(void) const { return std::sqrt(x*x+y*y); }

		void Normalize(void)
		{
			*this /= Length();
		}

		std::string to_string() const {
			std::ostringstream oss;
			oss << "\t" << x << "\t" << y;
			return oss.str();
		}

		void from_string(std::string point_str) {
			point_str = point_str.substr(point_str.find_first_not_of("\t"));
			std::vector<std::string> strs = SketchUI::split_str(point_str, "\t");

			assert(strs.size() >= 2);
			x = std::stod(strs[0]);
			y = std::stod(strs[1]);
		}
	};

	class Polyline2D
	{
	public: // data
		int stroke_ind;
		int group_ind;
		std::vector<std::pair<Point2D, std::int64_t>> points;
		double width;
		double sampling_interval;

	public: 
		Polyline2D() : sampling_interval(5), width(epsilon_small) { }

		// attributes
		// calculate total length of the polyline as the sum of each line segment's length
		double totalLen(void) const
		{
			if(points.size() == 0)
			{
				puts("Can't calculate polyline total length! There's no points.");
				return -1.0;
			}
			else if(points.size() == 1)
			{
				return 0.0;
			}
			else
			{
				double len = 0.0;
				for(size_t i = 0; i < points.size() - 1; i++)
				{
					len += (points[i+1].first - points[i].first).Length();
				}
				return len;
			}
		}

		void computeSegLen(std::vector<double>& segLen) const
		{
			if(points.size() < 2)
			{
				puts("Can't calculate segment lengths! There's no segments.");
				return;
			}

			segLen.resize(points.size() - 1);
			for(size_t i = 0; i < points.size() - 1; i++)
			{
				segLen[i] = (points[i].first - points[i+1].first).Length();
			}
		}

		void computeSegLen2(std::vector<double>& segLen) const
		{
			if(points.empty())
			{
				puts("Can't calculate segment lengths! There's no segments.");
				return;
			}

			segLen.resize(points.size());
			segLen[0] = 0;
			for(size_t i = 1; i < points.size(); i++)
			{
				segLen[i] = (points[i].first - points[i - 1].first).Length();
			}
		}

		void reparameterize(double sample_step_length, size_t min_samples = 3) {
			double total = totalLen();
			size_t num_samples = std::max(min_samples, (size_t)ceil(total / sample_step_length));
			sample_step_length = total / num_samples;
			std::vector<double> segLen;
			computeSegLen(segLen);

			sampling_interval = sample_step_length;

			if (points.empty()) return;

			std::vector<std::pair<Point2D, std::int64_t>> reparam_points;
			reparam_points.reserve(num_samples + 2);
			reparam_points.emplace_back(points[0]);

			double acc = 0.;
			double local_t;
			for (size_t i = 0; i < segLen.size(); i++) {
				acc += segLen[i];
				double new_acc = acc;
				size_t j = 1;

				// reparameterize the current step
				while (acc >= j * sample_step_length) {
					local_t = (segLen[i] - (acc - j * sample_step_length)) / segLen[i];
					Point2D reparam_p((1 - local_t) * points[i].first.x + local_t * points[i + 1].first.x,
						(1 - local_t) * points[i].first.y + local_t * points[i + 1].first.y);
					reparam_p.tangent_x = (1 - local_t) * points[i].first.tangent_x + local_t * points[i + 1].first.tangent_x;
					reparam_p.tangent_y = (1 - local_t) * points[i].first.tangent_y + local_t * points[i + 1].first.tangent_y;
					double length = std::sqrt(reparam_p.tangent_x * reparam_p.tangent_x + reparam_p.tangent_y * reparam_p.tangent_y);
					if (length > std::numeric_limits<double>::epsilon()) {
						reparam_p.tangent_x /= length;
						reparam_p.tangent_y /= length;
					}

					long double reparam_t = (1 - local_t) * points[i].second + local_t * points[i + 1].second;

					reparam_points.emplace_back(reparam_p, (std::int64_t)reparam_t);

					// left segment length
					new_acc = (acc - j * sample_step_length);

					j++;
				}

				acc = new_acc;
			}

			if ((points[points.size() - 1].first - reparam_points.back().first).Length() > 1.) {
				reparam_points.emplace_back(points[points.size() - 1]);
			}
			else {
				reparam_points.back() = points[points.size() - 1];
			}

			points = reparam_points;
		}

		// calculate normalized tangent directions at each point,
		// NOTE: there mustn't be duplicated points!!!
		void getTangentDirs(std::vector<Point2D>& tangentDirs) const
		{
			tangentDirs.resize(0);
			tangentDirs.reserve(points.size());

			if(points.size() < 2)
			{
				puts("Can't get polyline tangent vectors! Polyline points less than 2.");
				return;
			}

			for(size_t i = 0; i < points.size() - 1; i++)
			{
				tangentDirs.push_back(points[i+1].first-points[i].first);
				tangentDirs.back().Normalize();
			}
			tangentDirs.push_back(tangentDirs.back());
			for(size_t i = points.size() - 2; i > 0; i--)
			{
				tangentDirs[i] = tangentDirs[i-1] + tangentDirs[i];
				tangentDirs[i].Normalize();
			}
		}

		void reindex() {
			if (points.size() == 0) return;

			std::vector<std::pair<Point2D, std::int64_t>> dedup_points = points;
			points.clear();

			points.emplace_back(dedup_points.front());
			for (size_t i = 1; i < dedup_points.size(); i++) {
				if ((dedup_points[i].first - points.back().first).Length() > std::numeric_limits<double>::epsilon())
					points.emplace_back(dedup_points[i]);
			}
		}

		double get_t(size_t ind) const {
			assert(ind < points.size());

			std::vector<double> seg;
			computeSegLen2(seg);
			double all_seg = 0;
			double total = totalLen();
			for (auto &s : seg) {
				all_seg += s;
				s = all_seg / total;
			}

			return seg[ind];
		}

		bool is_roughly_circular() const {
			//std::cout << "\t In is_roughly_circular" << std::endl;
			if (points.empty()) return false;
			if ((points.front().first - points.back().first).Length() < std::numeric_limits<double>::epsilon())
				return true;

			auto to_glm_vec = [](Point2D p)->glm::dvec2 {
				return glm::dvec2(p.x, p.y);
			};

			double turning = 0;
			double min_angle_dist = 5 * M_PI;
			glm::dvec2 near_round_p = to_glm_vec(points.at(0).first);
			for (size_t i = 1; i + 1 < points.size(); i++) {
				glm::dvec2 t0 = to_glm_vec(points.at(i).first - points.at(i - 1).first);
				glm::dvec2 t1 = to_glm_vec(points.at(i + 1).first - points.at(i).first);

				t0 = glm::normalize(t0);
				t1 = glm::normalize(t1);

				glm::dvec2 up0(-t0.y, t0.x);
				double sign = (glm::dot(up0, t1) > 0) ? 1 : -1;

				turning += sign * std::acos(std::max(-1., std::min(glm::dot(t0, t1), 1.)));

				double angle_dist = std::abs(std::abs(turning) - 2 * M_PI);
				if (angle_dist < min_angle_dist) {
					near_round_p = to_glm_vec(points.at(i).first);
					min_angle_dist = angle_dist;
				}
			}

			bool angle_check = (std::abs(turning) > 250. / 180 * M_PI);
			bool dist_check = glm::distance(near_round_p, to_glm_vec(points.at(0).first)) < 20;

			//std::cout << "\t Turn: " << turning/M_PI * 180 << "; " << glm::distance(near_round_p, to_glm_vec(points.at(0).first)) << std::endl;

			return angle_check && dist_check;
		}
		
		void laplacian_smooth(size_t center_ind, double t_range) {
			auto to_glm_vec = [](Point2D p)->glm::dvec2 {
				return glm::dvec2(p.x, p.y);
			};

			std::vector<double> seg;
			computeSegLen2(seg);
			double all_seg = 0;
			double total = totalLen();
			for (auto &s : seg) {
				all_seg += s;
				s = all_seg / total;
			}

			std::vector<std::pair<Point2D, std::int64_t>> smooth_points = points;
			for (int i = 0; i < smooth_points.size(); i++) {
				if ((std::abs(seg[i] - seg[center_ind]) < t_range || std::abs(i - (int)center_ind) == 1)
					&& i != 0 && i + 1 != smooth_points.size()) {
					double wij = 1. / (points[i].first - points[i - 1].first).Length();//1
					double wik = 1. / (points[i].first - points[i + 1].first).Length();//1

					glm::dvec2 smoothed = (wij * to_glm_vec(points[i - 1].first) + wik * to_glm_vec(points[i + 1].first))
						/ (wij + wik);

					smooth_points[i].first = Point2D(smoothed.x, smoothed.y);
				}
			}

			points = smooth_points;
		}

		void laplacian_smooth(size_t center_ind, size_t ind_range) {
			auto to_glm_vec = [](Point2D p)->glm::dvec2 {
				return glm::dvec2(p.x, p.y);
			};

			std::vector<std::pair<Point2D, std::int64_t>> smooth_points = points;
			for (int i = 0; i < smooth_points.size(); i++) {
				if ((std::abs(i - (int)center_ind) <= ind_range)
					&& i != 0 && i + 1 != smooth_points.size()) {
					double wij = 1. / (points[i].first - points[i - 1].first).Length();//1
					double wik = 1. / (points[i].first - points[i + 1].first).Length();//1

					glm::dvec2 smoothed = (wij * to_glm_vec(points[i - 1].first) + wik * to_glm_vec(points[i + 1].first))
						/ (wij + wik);

					smooth_points[i].first = Point2D(smoothed.x, smoothed.y);
				}
			}

			points = smooth_points;
		}

		void laplacian_move(size_t center_ind, Point2D center_new, size_t ind_range, size_t num_itr = 1) {
			auto to_glm_vec = [](Point2D p)->glm::dvec2 {
				return glm::dvec2(p.x, p.y);
			};

			// First, compute laplacians before we move
			struct Laplacian {
				double lap_mag;
				double angle;

				Laplacian() : lap_mag(0), angle(0) { }
			};
			std::vector<Laplacian> laplacians;
			laplacians.resize(points.size());

			auto compute_laplacians = [&](std::vector<Laplacian> &laplacians) {
				for (size_t i = 1; i + 1 < points.size(); i++) {
					glm::dvec2 neighbor_vec = to_glm_vec(points[i + 1].first - points[i - 1].first);
					glm::dvec2 neighbor_dir = glm::normalize(neighbor_vec);
					glm::dvec2 lap = to_glm_vec(points[i].first) - 0.5 * to_glm_vec(points[i + 1].first + points[i - 1].first);
					glm::dvec2 lap_dir = glm::normalize(lap);

					glm::dvec2 neighbor_up(-neighbor_dir.y, neighbor_dir.x);
					double sign = (glm::dot(neighbor_up, lap_dir) > 0) ? -1 : 1;

					double dot_prod = std::max(-1., std::min(1., glm::dot(neighbor_dir, lap_dir)));

					laplacians[i].angle = sign * std::acos(dot_prod);
					laplacians[i].lap_mag = std::sqrt(glm::dot(lap, lap));
				}
			};
			compute_laplacians(laplacians);

			// Now we move
			points[center_ind].first = center_new;

			// Finally we smooth
			for (size_t k = 0; k < num_itr; k++) {
				std::vector<std::pair<Point2D, std::int64_t>> smooth_points = points;
				for (int i = 0; i < smooth_points.size(); i++) {
					if (i != 0 && i + 1 != smooth_points.size()) {
						glm::dvec2 avg_neighbor = 0.5 * to_glm_vec(points[i + 1].first + points[i - 1].first);
						glm::dvec2 neighbor_vec = to_glm_vec(points[i + 1].first - points[i - 1].first);
						glm::dvec2 neighbor_dir = glm::normalize(neighbor_vec);
						glm::dmat2x2 rotation(std::cos(laplacians[i].angle), -std::sin(laplacians[i].angle), std::sin(laplacians[i].angle), std::cos(laplacians[i].angle));
						glm::dvec2 lap = rotation * neighbor_dir;
						lap *= laplacians[i].lap_mag;

						glm::dvec2 smoothed = avg_neighbor + lap;

						if (!(std::abs(i - (int)center_ind) <= ind_range || std::abs(i - (int)center_ind) == 1)) {
							smoothed = 0.1 * smoothed + 0.9 * to_glm_vec(points[i].first);
						}

						smooth_points[i].first = Point2D(smoothed.x, smoothed.y);

					}
				}

				points = smooth_points;
			}
		}

		Point2D get_endpoint_inward_tangent(size_t end_ind) const
		{
			Point2D tangent(0, 0);
			if (end_ind == 0) {
				tangent = points[1].first-points[0].first;
				tangent.Normalize();
			}
			else if (end_ind == points.size() - 1) {
				tangent = points[points.size() - 2].first-points[points.size() - 1].first;
				tangent.Normalize();
			}

			return tangent;
		}

		// calculate center point of the polyline by equally weighted averaging
		Point2D center(void) const
		{
			if(points.size() == 0)
			{
				puts("Can't calculate polyline center! There's no points.");
				return Point2D(0.0, 0.0);
			}
			else if(points.size() == 1)
			{
				return points[0].first;
			}
			else
			{
				Point2D c(0.0, 0.0);
				for(auto vi = points.begin(); vi != points.end(); vi++)
				{
					c += vi->first;
				}
				c.x /= points.size();
				c.y /= points.size();
				return c;
			}
		}

		// compute normalized principal direction of the polyline by PCA
		Point2D principalDir(void) const
		{
			if(points.size() < 2)
			{
				puts("Can't compute principal direction! Polyline points less than 2.");
				return Point2D(0.0, 0.0);
			}
			else if(points.size() == 2)
			{
				return points[1].first - points[0].first;
			}
			else
			{
				Point2D centerPoint = this->center();
				double a = 0.0, b = 0.0, d = 0.0;
				for(size_t i = 0; i < points.size(); i++)
				{
					double x = points[i].first.x - centerPoint.x;
					double y = points[i].first.y - centerPoint.y;
					a += x * x;
					b += x * y;
					d += y * y;
				}
				double c = b;
				double T = a + d, D = a*d - b*c;
				double L1 = T/2 + sqrt(T*T/4-D);

				if(c != 0.0)
				{
					Point2D dir(L1-d, c);
					dir.Normalize();
					return dir;
				}
				else if(b != 0.0)
				{
					Point2D dir(b, L1-a);
					dir.Normalize();
					return dir;
				}
				else
				{
					return Point2D(1.0, 0.0);
				}
			}
		}

		std::string to_string() const {
			std::string poly_str = "{\n";
			poly_str += "\t#" + std::to_string(stroke_ind) + "\t" + std::to_string(group_ind) + "\n";
			poly_str += "\t@" + std::to_string(width) + "\n";

			for (auto const & p : points) {
				poly_str += p.first.to_string();

				std::ostringstream oss;
				oss << "\t" << p.second;

				poly_str += oss.str() + "\n";
			}

			poly_str += "}\n";

			return poly_str;
		}

		void from_string(std::string poly_str) {
			stroke_ind = group_ind = -1;

			poly_str = poly_str.substr(poly_str.find_first_not_of("{"));
			std::vector<std::string> strs = SketchUI::split_str(poly_str, "\n");

			width = epsilon_small;
			for (size_t i = 0; i < strs.size(); i++) {
				if (strs[i].empty())
					continue;
				if (strs[i].find("#") != std::string::npos) {
					strs[i] = strs[i].substr(strs[i].find_first_of("#"));
					assert(strs[i].length() > 1);
					strs[i] = strs[i].substr(1);

					std::vector<std::string> strs_digit = SketchUI::split_str(strs[i], "\t");
					stroke_ind = std::stoi(strs_digit[0]);
					group_ind = std::stoi(strs_digit[1]);

					continue;
				}

				if (strs[i].find("@") != std::string::npos) {
					strs[i] = strs[i].substr(strs[i].find_first_of("@"));
					assert(strs[i].length() > 1);
					strs[i] = strs[i].substr(1);

					width = std::stoi(strs[i]);

					continue;
				}

				Point2D point;
				point.from_string(strs[i]);

				uint64_t time;
				{
					strs[i] = strs[i].substr(strs[i].find_first_not_of("\t"));
					std::vector<std::string> strs_digit = SketchUI::split_str(strs[i], "\t");

					time = std::stoul(strs_digit[2]);
				}

				// Deduplicate
				if (points.size() == 0 ||
					(points.back().first - point).SqLength() > std::numeric_limits<float>::epsilon()) {
					Point2D prev;
					if (points.size() > 0)
						prev = points.back().first;

					points.emplace_back(std::make_pair(point, time));
					Point2D cur = points.back().first;
				}
			}

			sampling_interval = 0;

			for (size_t i = 0; i + 1 < points.size(); i++) {
				sampling_interval += (points[i + 1].first - points[i].first).Length();
			}

			sampling_interval = (points.size() >= 2) ? sampling_interval/points.size() : 0;
		}

		std::string to_svg_string() const {
			std::string svg_str = "";
			if (points.empty()) return svg_str;

			svg_str += "\t<path d=\"";

			for (size_t i = 0; i < points.size(); i++) {
				svg_str += (i == 0) ? "M " : "L ";
				svg_str += std::to_string(points[i].first.x) + " " + std::to_string(points[i].first.y);
				if (i + 1 < points.size())
					svg_str += " ";
			}

			svg_str += "\" fill=\"none\" stroke=\"black\" stroke-width=\"" + std::to_string(this->width) + "\" />\n";

			return svg_str;
		}
	};
}