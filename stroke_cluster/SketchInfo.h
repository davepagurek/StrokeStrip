#pragma once

#include "Polyline2D.h"

#include <memory>

#include <cassert>
#include <iostream>
#include <iostream>

namespace SketchUI
{
#define POLYLINESAMPLETHRES 1e-6 // the min squared length of polyline segments
#define PICKNODETHRES 1e-4

	enum SketchState
	{
		SS_UNCLASSIFIED, // the initial state, waiting for sketch input 
		SS_HEMLINE_WRINKLE, // waiting for picking up wrinkle edge end node and adjust
		SS_ADJUST_HEMLINE_WRINKLE_EDGES // adjusting wrinkle edge directions by draging the end node
		//SS_HEMLINE_WRINKLE_ING
	};

	struct Calibrator {
		std::vector<Point2D> standard_corners;
		std::vector<Point2D> raw_corners;
		Point2D translation;
		double scale;
	};

	class SketchInfo
	{
	protected: // data
		bool m_isSketching; // whether the user is drawing the sketching
		SketchState m_sketchState;
		std::vector<std::string> sketchStateString;
		Polyline2D curSketchPolyline; // the curretnly drawing sketch curve, in NDCS [-1,1]*[-1,1]

		int pickedPolylineI;

	public: // attribute
		std::vector<Polyline2D> sketchedPolylines; // previously drawn sketch curves

		bool isSketching(void) const { return m_isSketching; }
		SketchState sketchState(void) const { return m_sketchState; }

		size_t getCurSketchPolylineSize(void) const { return curSketchPolyline.points.size(); }
		const Point2D& getCurSketchPolyline(size_t i) const 
		{
			assert(i < curSketchPolyline.points.size());
			return curSketchPolyline.points[i].first; 
		}

		size_t getSketchedPolylineSize(void) const { return sketchedPolylines.size(); }
		size_t getSketchedPolylineSize(size_t i) const
		{
			assert(i < sketchedPolylines.size());
			return sketchedPolylines[i].points.size();
		}
		Polyline2D& getSketchedPolyline(size_t lineI)
		{
			assert(lineI < sketchedPolylines.size());
			return sketchedPolylines[lineI];
		}
		const Point2D& getSketchedPolyline(size_t lineI, size_t pointI) const
		{
			assert(lineI < sketchedPolylines.size());
			assert(pointI < sketchedPolylines[lineI].points.size());
			return sketchedPolylines[lineI].points[pointI].first;
		}

		double thickness;

	public: // constructor
		SketchInfo(void) :
			thickness(2)
		{
			sketchStateString.push_back("SS_UNCLASSIFIED");
			sketchStateString.push_back("SS_HEMLINE_WRINKLE");
			sketchStateString.push_back("SS_ADJUST_HEMLINE_WRINKLE_EDGES");
			//sketchStateString.push_back("SS_HEMLINE_WRINKLE_ING");
			init();
		}

		void init(void)
		{
			m_isSketching = false;
			changeStateTo(SS_UNCLASSIFIED);
			pickedPolylineI = -1;
			clearSketchPolylines();
		}

	public: // operations
		void changeStateTo(SketchState newState)
		{
			m_sketchState = newState;
			//printf("Sketch state changed to: %s\n", sketchStateString[m_sketchState].c_str());
		}

		// start sampling polyline points for the currently drawing sketch curve
		void startSampling(const Point2D& startPoint, std::int64_t timestamp)
		{
			assert((!m_isSketching) && (curSketchPolyline.points.size() == 0));

			m_isSketching = true;
			curSketchPolyline.points.push_back(std::make_pair(startPoint, timestamp));

			//puts("Start sketching...");
		}

		// sample polyline point p for the currently drawing sketch curve
		// if p is not too close (squared_dist > POLYLINESAMPLETHRES) to the previously sampled point
		void samplePolylinePoint(const Point2D& p, std::int64_t timestamp)
		{
			assert(m_isSketching && (curSketchPolyline.points.size() > 0));
			// Debug: fails when running remotely
			//assert(curSketchPolyline.points.back().second < timestamp);
			
			std::cout << "sampling: " << p.x << ", " << p.y << " - " << (p - curSketchPolyline.points.back().first).SqLength() << std::endl;
			if((p - curSketchPolyline.points.back().first).SqLength() >= POLYLINESAMPLETHRES)
			{
				curSketchPolyline.points.push_back(std::make_pair(p, timestamp));
			}
		}

		// end sampling polyline points for the currently drawing sketch curve
		// the polyline just drawn will be saved into sketchedPolylines,
		// curSketchPolyline will be cleared
		void endSampling(void)
		{
			m_isSketching = false;
			sketchedPolylines.push_back(curSketchPolyline);
			curSketchPolyline.points.resize(0);

			printf("Sampled %u points for polyline %u.\n", sketchedPolylines.back().points.size(), sketchedPolylines.size() - 1);
		}

		void deleteLastPolyline(void)
		{
			if(sketchedPolylines.size() > 0)
			{
				sketchedPolylines.pop_back();
				printf("Last polyline deleted, %u polylines remained.\n", sketchedPolylines.size());
			}
			else
			{
				puts("No polylines to delete.");
			}			
		}

		void clearSketchPolylines(void)
		{
			assert(!m_isSketching);
			curSketchPolyline.points.resize(0);
			sketchedPolylines.resize(0);

			//puts("Sketches cleared.");
		}

		bool pickEndNode(size_t polylineI, const Point2D& p)
		{
			assert(polylineI < sketchedPolylines.size());
			assert(sketchedPolylines[polylineI].points.size() > 0);

			if((sketchedPolylines[polylineI].points.back().first - p).SqLength() <= PICKNODETHRES)
			{
				pickedPolylineI = polylineI;
				return true;
			}
			else
			{
				return false;
			}
		}

		std::string to_string() const {
			std::string poly_str = "";

			// Thickness
			poly_str += "@" + std::to_string(thickness) + "\n";

			for (auto const &sketch : sketchedPolylines) {
				if (sketch.points.size() == 1) {
					std::cout << "Warning: Skipping sketch with only one point." << std::endl;
					continue;
				}
				poly_str += sketch.to_string();
			}

			return poly_str;
		}

		void from_string(std::string sketch_str) {
			std::string pattern("}");

			// Read thickness
			if (sketch_str.find_first_of("@") != std::string::npos) {
				std::string thickness_str = sketch_str.substr(sketch_str.find_first_of("@") + 1);
				std::vector<std::string> strs_digit = SketchUI::split_str(thickness_str, "\n");
				assert(strs_digit.size() >= 1);
				thickness = std::stod(strs_digit[0]);
			}

			std::vector<std::string> strs = split_str(sketch_str, pattern);
			for (size_t i = 0; i < strs.size(); i++) {
				if (strs[i].find_first_of("{") == std::string::npos)
					continue;

				Polyline2D poly;
				poly.from_string(strs[i].substr(strs[i].find_first_of("{")));
				sketchedPolylines.emplace_back(poly);
			}
		}

		std::string to_svg_string(size_t width, size_t height) const {
			std::string svg_str = "<?xml version=\"1.0\" encoding=\"utf-8\" ?>\n";
			svg_str += "<svg baseProfile=\"full\" height=\"" + std::to_string(height) + 
				"\" version=\"1.1\" viewBox=\"0,0," + std::to_string(width) + "," + std::to_string(height) + 
				"\" width=\"" + std::to_string(width) + 
				"\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:ev=\"http://www.w3.org/2001/xml-events\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

			for (auto const &sketch : sketchedPolylines) {
				if (sketch.points.size() == 1) {
					std::cout << "Warning: Skipping sketch with only one point." << std::endl;
					continue;
				}
				svg_str += sketch.to_svg_string();
			}
			svg_str += "</svg>";

			return svg_str;
		}
	};
}

typedef SketchUI::SketchInfo Capture;
typedef SketchUI::Polyline2D Sketch;
typedef SketchUI::Point2D SketchPoint;