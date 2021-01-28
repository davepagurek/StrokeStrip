#include "Renderer.h"

#include <iostream>
#include <random>
#include <algorithm>

#define NANOVG_GL3_IMPLEMENTATION
#include "nanovg/src/nanovg_gl.h"

#include "DrawLog.h"

Renderer renderer;
Capture capture;
extern DrawLog draw_log;
extern double epsilon_small;

extern std::vector<StrokeFamily> demo_families;

std::string hex_palette[] = {"#c583d9", "#9de993", "#ad5ba7", "#58f2c3", "#c94e76", 
"#57b870", "#d271bb", "#3d8e43", "#f375a0", "#00d4db", "#d55969", 
"#64e4ff", "#e68e54", "#54aeff", "#eed973", "#008bd4", "#fdb968", 
"#6f73bc", "#bee481", "#8f6aaf", "#8bb658", "#bbabff", "#658424", 
"#f1b9ff", "#989c37", "#bb568c", "#cfe395", "#c45372", "#008c6f", 
"#ff94af", "#488656", "#ffc0e1", "#6c8139", "#a9d0ff", "#a98c2c", 
"#0294bb", "#d7b555", "#4a7ea6", "#ffa372", "#96ddff", "#b06735", 
"#bae0f9", "#b76246", "#b0e4e9", "#ff9b9c", "#378673", "#ffa98e", 
"#6a79a2", "#ffcca9", "#53827e", "#ffcdcf", "#648168", "#e7d3fa", 
"#93744c", "#cde1bc", "#9b6a8a", "#efd8b1", "#928098", "#7c7b5c", "#a6686c"};

int Renderer::setup_gui(int _width, int _height) {
	GLenum err = 0;
    /*********************************************
     * GLFW SETUP
     *********************************************/
    err = glfwInit();
    if (!err)
    {
        fputs("Failed to load the GLFW library", stderr);
        exit(EXIT_FAILURE);
    }

#ifdef __APPLE__
     glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	 glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	 glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	 glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#endif
    
    /*********************************************
     * STATE SETUP (initialize gl context)
     *********************************************/
    // must be setup before glew so that a valid openGL
    // context exists (created with the window)

	 {
		width  = _width;
		height = _height;
		/* As of right now we only have one window */
		window = glfwCreateWindow(width, height, "StrokeAggregator", NULL, NULL);
		if (!window)
		{
			glfwTerminate();
			fputs("failed to initialize window", stderr);
			return 1; // error
		}
		glfwMakeContextCurrent(window);

		int fbWidth, fbHeight;

		glfwGetFramebufferSize(window, &fbWidth, &fbHeight);
		// Calculate pixel ration for hi-dpi devices.
		pxRatio = (float)fbWidth / (float)width;
	 }

    /*********************************************
     * GLEW SETUP
     *********************************************/
#ifdef __APPLE__
     glewExperimental = GL_TRUE;
#endif
    err = glewInit();
    if (err != GLEW_OK)
    {
        fputs("Failed to initialize the GLEW library", stderr);
        exit(EXIT_FAILURE);
    }

    /*********************************************
     * NANOVG SETUP
     *********************************************/
    vg = nvgCreateGL3(NVG_ANTIALIAS | NVG_STENCIL_STROKES);

	int font;
	font = nvgCreateFont(vg, "times", "C:/Windows/Fonts/times.ttf");
	if (font == -1) {
		printf("Could not add font regular.\n");
		exit(EXIT_FAILURE);
	}
	/*********************************************
     * SETUP IMGUI
     *********************************************/
    ImGui_ImplGlfwGL3_Init(window, false);

    /*********************************************
     * SET GL STATE
     *********************************************/
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	/*********************************************
     * SET CALLBACKS
     *********************************************/
	setup_label_callback();

	return 0;
}

void Renderer::cleanup_gui() {
	if (!vg)
		return;

	nvgDeleteGL3(vg);

	/*********************************************
     * CLEAN UP
     *********************************************/

	ImGui_ImplGlfwGL3_Shutdown();

    glfwTerminate();
}

// Camera functions
void Renderer::look_at(glm::dvec2 center) {
	glm::dvec2 in_center(width / 2, height / 2);
	translation = in_center - center;
}

void Renderer::zoom(double ratio) {
	scale *= ratio;
}
//

void Renderer::begin_draw() {
	if (!vg)
		return;

	glViewport(0, 0, width, height);
	glClearColor(1.f, 1.f, 1.f, 0.f);
	glClear(GL_COLOR_BUFFER_BIT);

	nvgBeginFrame(vg, width, height, pxRatio);
}

void Renderer::end_draw() {
	if (!vg)
		return;
	nvgEndFrame(vg);

	glfwSwapBuffers(window);
	glfwPollEvents();
}

glm::dvec2 Renderer::transform_point(glm::dvec2 p) const {
	glm::dvec2 in_center(width / 2, height / 2);
	p -= in_center;
	return glm::dvec2((p + translation).x * scale.x, (p + translation).y * scale.y) + in_center;
}

void Renderer::draw_capture(Capture capture, NVGcolor color) {
	if (!vg)
		return;

	// For debug visualization
	float alpha = 1.f;
	if (alpha_switch)
		alpha = 0.5;
	//

	std::vector<std::pair<glm::dvec2, glm::dvec2>> segments;
	for (size_t i = 0; i < capture.getSketchedPolylineSize(); i++) {
		for (size_t j = 0; j + 1 < capture.getSketchedPolylineSize(i); j++) {
			segments.emplace_back(glm::dvec2(capture.getSketchedPolyline(i, j).x, capture.getSketchedPolyline(i, j).y), 
				glm::dvec2(capture.getSketchedPolyline(i, j + 1).x, capture.getSketchedPolyline(i, j + 1).y));
		}
	}
	draw_segments(segments, nvgRGBAf(color.r, color.g, color.b, alpha));
}

void Renderer::draw_strokes(std::vector<Sketch> strokes, NVGcolor color) {
	if (!vg)
		return;

	// For debug visualization
	float alpha = 1.f;
	if (alpha_switch)
		alpha = 0;
	//

	std::vector<std::pair<glm::dvec2, glm::dvec2>> segments;
	for (size_t i = 0; i < strokes.size(); i++) {
		for (size_t j = 0; j + 1 < strokes[i].points.size(); j++) {
			segments.emplace_back(glm::dvec2(strokes[i].points[j].first.x, strokes[i].points[j].first.y), 
				glm::dvec2(strokes[i].points[j + 1].first.x, strokes[i].points[j + 1].first.y));
		}
	}
	draw_segments(segments, nvgRGBAf(color.r, color.g, color.b, alpha));
}

void Renderer::draw_stroke(Sketch stroke, NVGcolor color) {
	if (!vg)
		return;

	std::vector<std::pair<glm::dvec2, glm::dvec2>> segments;
	for (size_t j = 0; j + 1 < stroke.points.size(); j++) {
		segments.emplace_back(glm::dvec2(stroke.points[j].first.x, stroke.points[j].first.y), 
			glm::dvec2(stroke.points[j + 1].first.x, stroke.points[j + 1].first.y));
	}

	draw_segments(segments, color);
}

void Renderer::draw_dot(glm::dvec2 dot, NVGcolor color) {
	if (!vg)
		return;

	glm::dvec2 transform_dot = transform_point(dot);

	nvgBeginPath(vg);
	nvgCircle(vg, transform_dot.x, transform_dot.y, 3.f /** scale.x*/);
	nvgFillColor(vg, color);
	nvgFill(vg);
}

void Renderer::draw_dots(std::vector<glm::dvec2> dots, NVGcolor color) {
	for (glm::dvec2 const &d : dots) {
		draw_dot(d, color);
	}
}

void Renderer::draw_dots(std::vector<glm::dvec2> dots, std::vector<NVGcolor> colors) {
	assert(colors.size() >= dots.size());
	for (size_t i = 0; i < dots.size(); i++) {
		glm::dvec2 d = dots[i];
		NVGcolor color = colors[i];

		draw_dot(d, color);
	}
}

void move_seg(NVGcontext* vg, std::pair<glm::dvec2, glm::dvec2> seg) {
	nvgMoveTo(vg, seg.first.x, seg.first.y);
	nvgLineTo(vg, seg.second.x, seg.second.y);
}

void Renderer::draw_segment(std::pair<glm::dvec2, glm::dvec2> seg, NVGcolor color) {
	if (!vg)
		return;

	std::pair<glm::dvec2, glm::dvec2> transform_seg = std::make_pair(transform_point(seg.first),
		transform_point(seg.second));

	nvgBeginPath(vg);
	move_seg(vg, transform_seg);
	nvgStrokeColor(vg, color);
	nvgStrokeWidth(vg, epsilon_small);
	nvgStroke(vg);
}

void Renderer::draw_segments(std::vector<std::pair<glm::dvec2, glm::dvec2>> segments, NVGcolor color) {
	if (!vg)
		return;

	nvgBeginPath(vg);
	for (auto const &seg : segments) {
		std::pair<glm::dvec2, glm::dvec2> transform_seg = std::make_pair(transform_point(seg.first),
			transform_point(seg.second));
		move_seg(vg, transform_seg);
	}
	nvgStrokeColor(vg, color);
	nvgStrokeWidth(vg, epsilon_small);
	nvgStroke(vg);
}

void Renderer::draw_segments(std::vector<std::pair<glm::dvec2, glm::dvec2>> segments, std::vector<NVGcolor> colors) {
	assert(colors.size() >= segments.size());
	for (size_t i = 0; i < segments.size(); i++) {
		draw_segment(segments[i], colors[i]);
	}
}

void Renderer::draw_circle(glm::dvec3 circle, NVGcolor color) {
	if (!vg)
		return;

	glm::dvec2 transform_dot = transform_point(glm::dvec2(circle.x, circle.y));

	nvgBeginPath(vg);
	nvgCircle(vg, transform_dot.x, transform_dot.y, circle.z * scale.x);
	nvgStrokeColor(vg, color);
	nvgStroke(vg);
}

void Renderer::draw_circles(std::vector<glm::dvec3> circles, NVGcolor color) {
	for (auto const &c : circles) {
		draw_circle(c, color);
	}
}

void Renderer::setup_stroke_colors(size_t num_labels) {
	std::mt19937 generator(0x1234);

	size_t max_label = num_labels;
	for (size_t label = 0; label <= max_label; label++) {
		NVGcolor color;

		bool accept = false;
		while (!accept) {
			color = nvgRGB(generator(), generator(), generator());
			accept = true;
			
			for (auto const &c : stroke_colors) {
				if (((color.r - c.second.r) * (color.r - c.second.r) + 
					(color.g - c.second.g) * (color.g - c.second.g) + 
					(color.b - c.second.b) * (color.b - c.second.b)) < 0.1 * 0.1) {
					accept = false;
					break;
				}
			}
		}

		stroke_colors.emplace(label, color);
	}
}

NVGcolor getHexColor(std::string hex)
{
	NVGcolor color = nvgRGBf(0.0f, 0.0f, 0.0f);
    const int length = hex.length();
    if (!hex[0] == '#' || (length != 7 && length != 9)) {
        return color;
    }

    bool ok;
    unsigned char a,r,g,b;
    if (length == 7) {
        a = 255;
		r = std::stoi (hex.substr(1, 2),nullptr,16);
        g = std::stoi (hex.substr(3, 2),nullptr,16);
        b = std::stoi (hex.substr(5, 2),nullptr,16);
    } else {
        a = std::stoi (hex.substr(1, 2),nullptr,16);
        r = std::stoi (hex.substr(3, 2),nullptr,16);
        g = std::stoi (hex.substr(5, 2),nullptr,16);
        b = std::stoi (hex.substr(7, 2),nullptr,16);
    }
    color = nvgRGBA(r, g, b, a);

	return color;
}

void Renderer::setup_stroke_colors() {
	std::vector<std::string> palette_vec(hex_palette, hex_palette + sizeof(hex_palette) / sizeof(hex_palette[0]));
	for (size_t label = 0; label < 120; label++) {
		size_t ind = (label % palette_vec.size());
		stroke_colors.emplace(label, getHexColor(palette_vec[ind]));
	}
}

void Renderer::setup_stroke_colors(std::vector<uint16_t> const &vertex_labels) {
	std::mt19937 generator(0x1234);

	size_t max_label = *std::max_element(vertex_labels.cbegin(), vertex_labels.cend());
	for (size_t label = 0; label <= max_label; label++) {
		NVGcolor color;

		bool accept = false;
		while (!accept) {
			color = nvgRGB(generator(), generator(), generator());
			accept = true;
			
			for (auto const &c : stroke_colors) {
				if (((color.r - c.second.r) * (color.r - c.second.r) + 
					(color.g - c.second.g) * (color.g - c.second.g) + 
					(color.b - c.second.b) * (color.b - c.second.b)) < 0.1 * 0.1) {
					accept = false;
					break;
				}
			}
		}

		stroke_colors.emplace(label, color);
	}
}

/********************************************************/
/********************* Callbacks ************************/
/********************************************************/

void char_callback(GLFWwindow* window, unsigned int c) {
	if (ImGui::GetIO().WantTextInput) {
		ImGui_ImplGlfwGL3_CharCallback(window, c);
		return;
	}
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mode)
{
    if (ImGui::GetIO().WantTextInput) {
		ImGui_ImplGlfwGL3_KeyCallback(window, key, scancode, action, mode);
		return;
	}

	if (action != GLFW_RELEASE)
		return;

    switch(key)
    {
    case GLFW_KEY_LEFT:
        renderer.highlight_cluster = (renderer.highlight_cluster + 1) % renderer.num_clusters;
		renderer.to_delete = false;
        break;
    case GLFW_KEY_RIGHT:
        renderer.highlight_cluster = (renderer.highlight_cluster - 1 + renderer.num_clusters) % renderer.num_clusters;
		renderer.to_delete = false;
        break;
	case GLFW_KEY_D:
		renderer.to_delete = !renderer.to_delete;
        break;
	case GLFW_KEY_A:
		renderer.alpha_switch = !renderer.alpha_switch;
        break;
	case GLFW_KEY_R:
		renderer.translation = glm::dvec2(0, 0);
		renderer.scale = glm::dvec2(1, 1);
        break;
	case GLFW_KEY_EQUAL:
		renderer.scale += glm::dvec2(0.05, 0.05);
		break;
	case GLFW_KEY_MINUS:
		renderer.scale -= glm::dvec2(0.05, 0.05);
		break;
	case GLFW_KEY_C:
		renderer.to_color = !renderer.to_color;
        break;
	case GLFW_KEY_0:
	case GLFW_KEY_1:
	case GLFW_KEY_2:
	case GLFW_KEY_3:
	case GLFW_KEY_4:
		renderer.demo_selected = key - GLFW_KEY_0;
		break;
	default: break;
    }
    
}

glm::dvec2 screen_to_world(glm::dvec2 p, double width, double height, glm::dvec2 translation, glm::dvec2 scale) {
	glm::dvec2 in_center(width / 2, height / 2);
	glm::dvec2 world_pos = p - in_center;
	world_pos.x /= scale.x;
	world_pos.y /= scale.y;
	world_pos -= translation;
	world_pos += in_center;

	return world_pos;
}

static void mousePos_callback(GLFWwindow* win, double x, double y)
{
	if (ImGui::GetIO().WantCaptureMouse) {
		return;
	}

	glm::dvec2 screen_pos(x, y);
	glm::dvec2 world_pos = screen_to_world(screen_pos, renderer.width, renderer.height, renderer.translation, renderer.scale);

	double min_dist = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < capture.getSketchedPolylineSize(); i++) {
		double dist = std::numeric_limits<double>::infinity();
		for (size_t j = 0; j < capture.getSketchedPolyline(i).points.size(); j++) {
			double mouse_dist = glm::distance(world_pos, 
				glm::dvec2(capture.getSketchedPolyline(i).points[j].first.x, capture.getSketchedPolyline(i).points[j].first.y));
			if (mouse_dist < dist) {
				dist = mouse_dist;
			}
		}

		if (dist < min_dist) {
			min_dist = dist;
			renderer.highlight_stroke = i;
		}
	}

	// too far
	if (min_dist > 5) {
		renderer.highlight_stroke = -1;
	}
}

static void mouseButton_callback(GLFWwindow* win, int button, int action, int mods) {
	if (ImGui::GetIO().WantCaptureMouse) {
		return;
	}

	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
		double m_x, m_y;
		glfwGetCursorPos(win, &m_x, &m_y);

		glm::dvec2 screen_pos(m_x, m_y);
		glm::dvec2 world_pos = screen_to_world(screen_pos, renderer.width, renderer.height, renderer.translation, renderer.scale);

		renderer.look_at(world_pos);
	}

	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && renderer.to_demo) {
		if (ImGui::GetIO().WantCaptureMouse) {
			return;
		}

		double x, y;
		glfwGetCursorPos(win, &x, &y);
		glm::dvec2 screen_pos(x, y);
		glm::dvec2 world_pos = screen_to_world(screen_pos, renderer.width, renderer.height, renderer.translation, renderer.scale);

		double min_dist = std::numeric_limits<double>::infinity();
		for (auto g : demo_families[renderer.demo_selected].stroke_groups) {
			for (auto s : g.second) {
				double dist = std::numeric_limits<double>::infinity();
				for (size_t j = 0; j < s.points.size(); j++) {
					double mouse_dist = glm::distance(world_pos, 
						glm::dvec2(s.points[j].first.x, s.points[j].first.y));
					if (mouse_dist < dist) {
						dist = mouse_dist;
					}
				}

				if (dist < min_dist) {
					min_dist = dist;
					renderer.highlight_cluster = g.first;
				}
			}
		}

		// too far
		if (min_dist > 5) {
			renderer.highlight_cluster = -1;
		}
	}
}

static void scroll_callback(GLFWwindow *win, double x, double y) {
	if (ImGui::GetIO().WantCaptureMouse) {
		return;
	}

	renderer.scale += y * glm::dvec2(0.05, 0.05);
}

void Renderer::setup_label_callback() {
	glfwSetKeyCallback(window, key_callback);
	glfwSetCursorPosCallback(window, mousePos_callback);
	glfwSetCharCallback(window, char_callback);

	glfwSetMouseButtonCallback(window, mouseButton_callback);
	glfwSetScrollCallback(window, scroll_callback);
}

/********************************************************/
/****************** Visualizations **********************/
/********************************************************/

void Renderer::draw_capture_label(Capture capture, std::vector<uint16_t> const &vertex_labels) {
	if (!vg)
		return;

	unsigned char alpha = 255;
	if (alpha_switch)
		alpha = 50;

	if (stroke_colors.size() != vertex_labels.size()) {
		stroke_colors.clear();
		setup_stroke_colors(vertex_labels);
	}

	if (highlight_cluster < 0 && !vertex_labels.empty()) {
		num_clusters = *std::max_element(vertex_labels.cbegin(), vertex_labels.cend()) + 1;
		//setup_label_callback();
		highlight_cluster = 0;
	}

	int offset_total = width / 2.;
	int offset0 = -offset_total / 2.;
	int offset_step = (float)offset_total / stroke_colors.size();
	offset_step = offset0 = 0;

	for (size_t i = 0; i < capture.getSketchedPolylineSize(); i++) {
		NVGcolor color = nvgRGBA(20, 20, 20, alpha);//stroke_colors[vertex_labels[i]];
		//NVGcolor color = nvgRGBA(216, 76, 53, alpha);//stroke_colors[vertex_labels[i]];
		//color = stroke_colors[vertex_labels[i]];

		if (!vertex_labels.empty() &&
			vertex_labels[i] == highlight_cluster && to_delete) {
			continue;
		}

		std::vector<std::pair<glm::dvec2, glm::dvec2>> segments;
		for (size_t j = 0; j + 1 < capture.getSketchedPolylineSize(i); j++) {
			segments.emplace_back(glm::dvec2(capture.getSketchedPolyline(i, j).x, capture.getSketchedPolyline(i, j).y), 
				glm::dvec2(capture.getSketchedPolyline(i, j + 1).x, capture.getSketchedPolyline(i, j + 1).y));
		}

		draw_segments(segments, color);
	}

	if (!to_delete) {
		size_t count = 0;

		for (size_t i = 0; i < capture.getSketchedPolylineSize(); i++) {
			NVGcolor color = nvgRGB(252, 70, 0);//stroke_colors[vertex_labels[i]];
			NVGcolor select_color = nvgRGB(0, 70, 252);

			/*{
				draw_log.set_all_layers_filter(FilterAttribute(highlight_cluster, -1), false);
			}*/

			if (vertex_labels.empty() ||
				vertex_labels[i] != highlight_cluster) {
				continue;
			}

			{
				draw_log.set_all_layers_filter(FilterAttribute(highlight_cluster, -1), true);
			}

			if (renderer.highlight_stroke == i) {
				nvgFontSize(vg, 50.0f);
				nvgFontFace(vg, "times");
				nvgFillColor(vg, nvgRGBA(0, 70, 252, 255));
				nvgTextAlign(vg,NVG_ALIGN_LEFT|NVG_ALIGN_MIDDLE);
				nvgText(vg, 50, 50, (std::to_string(capture.getSketchedPolyline(i).stroke_ind) + " in " + std::to_string(vertex_labels[i])).c_str(), NULL);
			}

			std::vector<std::pair<glm::dvec2, glm::dvec2>> segments;
			for (size_t j = 0; j + 1 < capture.getSketchedPolylineSize(i); j++) {
				segments.emplace_back(glm::dvec2(capture.getSketchedPolyline(i, j).x, capture.getSketchedPolyline(i, j).y), 
					glm::dvec2(capture.getSketchedPolyline(i, j + 1).x, capture.getSketchedPolyline(i, j + 1).y));
			}
			draw_segments(segments, (renderer.highlight_stroke == i) ? select_color : color);

			nvgFontSize(vg, 20.0f);
			nvgFontFace(vg, "times");
			nvgFillColor(vg, nvgRGBA(252, 70, 0, 255));
			nvgTextAlign(vg,NVG_ALIGN_LEFT|NVG_ALIGN_MIDDLE);
			nvgText(vg, 10, 50 + 25 * count, (std::to_string(capture.getSketchedPolyline(i).stroke_ind) + " (" + std::to_string(i) + ")").c_str(), NULL);
			size_t j = capture.getSketchedPolylineSize(i) / 2.f;
			//nvgText(vg, capture.getSketchedPolyline(i, j).x + 10, capture.getSketchedPolyline(i, j).y + 10, std::to_string(i).c_str(), NULL);

			count++;
		}
	}
}

void Renderer::draw_capture_fit(Capture capture) {
	if (!vg)
		return;

	for (size_t i = 0; i < capture.getSketchedPolylineSize(); i++) {
		NVGcolor select_color = nvgRGB(0, 70, 252);

		if (renderer.highlight_stroke == i) {
			nvgFontSize(vg, 50.0f);
			nvgFontFace(vg, "times");
			nvgFillColor(vg, nvgRGBA(0, 70, 252, 255));
			nvgTextAlign(vg,NVG_ALIGN_LEFT|NVG_ALIGN_MIDDLE);
			nvgText(vg, 50, 50, 
				(std::to_string(capture.getSketchedPolyline(i).stroke_ind) + " in " + std::to_string(capture.getSketchedPolyline(i).group_ind)).c_str(), NULL);

			std::vector<std::pair<glm::dvec2, glm::dvec2>> segments;
			for (size_t j = 0; j + 1 < capture.getSketchedPolylineSize(i); j++) {
				segments.emplace_back(glm::dvec2(capture.getSketchedPolyline(i, j).x, capture.getSketchedPolyline(i, j).y), 
					glm::dvec2(capture.getSketchedPolyline(i, j + 1).x, capture.getSketchedPolyline(i, j + 1).y));
			}
			draw_segments(segments, select_color);
		}
	}
}

void Renderer::draw_capture_demo(std::vector<StrokeFamily> const &families) {
	if (!vg)
		return;

	unsigned char alpha = 255;
	if (alpha_switch)
		alpha = 50;

	if (stroke_colors.size() == 0) {
		setup_stroke_colors();
	}

	//if (highlight_cluster < 0) {
	//	num_clusters = families[demo_selected].stroke_groups.size();
	//	//setup_label_callback();
	//	highlight_cluster = 0;
	//}

	{
		nvgFontSize(vg, 40.f);
		nvgFontFace(vg, "times");
		nvgFillColor(vg, nvgRGBA(0, 70, 252, 255));
		nvgTextAlign(vg,NVG_ALIGN_LEFT|NVG_ALIGN_MIDDLE);
		std::string stage_name = "";

		switch (demo_selected)
		{
		case 0:
			stage_name = "Input";
			break;
		case 1:
			stage_name = "Angle";
			break;
		case 2:
			stage_name = "Sub-clusters";
			break;
		case 3:
			stage_name = "Final clusters";
			break;
		case 4:
			stage_name = "Final";
			break;
		default:
			break;
		}

		nvgText(vg, 20, 50, stage_name.c_str(), NULL);
	}

	for (auto g : families[demo_selected].stroke_groups) {
		NVGcolor color = nvgRGBA(20, 20, 20, alpha);//stroke_colors[vertex_labels[i]];
		//NVGcolor color = nvgRGBA(216, 76, 53, alpha);//stroke_colors[vertex_labels[i]];
		//color = stroke_colors[vertex_labels[i]];

		if (to_color) color = stroke_colors[g.first];

		for (auto s : g.second) {

			if (s.group_ind == highlight_cluster && to_delete) {
				continue;
			}

			std::vector<std::pair<glm::dvec2, glm::dvec2>> segments;
			for (size_t j = 0; j + 1 < s.points.size(); j++) {
				segments.emplace_back(glm::dvec2(s.points[j].first.x, s.points[j].first.y), 
					glm::dvec2(s.points[j + 1].first.x, s.points[j + 1].first.y));
			}

			draw_segments(segments, color);
		}
	}

	if (!to_delete) {
		for (auto g : families[demo_selected].stroke_groups) {
			NVGcolor color = nvgRGB(252, 70, 0);//stroke_colors[vertex_labels[i]];

			if (g.first != highlight_cluster) {
				continue;
			}

			/*{
				draw_log.set_all_layers_filter(FilterAttribute(highlight_cluster, -1), true);
			}*/

			std::vector<std::pair<glm::dvec2, glm::dvec2>> segments;
			for (auto s : g.second) {
				for (size_t j = 0; j + 1 < s.points.size(); j++) {
					segments.emplace_back(glm::dvec2(s.points[j].first.x, s.points[j].first.y), 
						glm::dvec2(s.points[j + 1].first.x, s.points[j + 1].first.y));
				}
			}
			draw_segments(segments, color);
		}
	}
}