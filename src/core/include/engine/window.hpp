#pragma once
#include <GL/glew.h>
#include "GLFW/glfw3.h"
#include "exceptions.hpp"


enum class Resolution {
	FHD,
	HD2K,
	UHD,
	UNKNOWN
};

ivec2 ires2(Resolution r);

struct WindowSettings {
	int width, height;
	string windowTitle;
	bool stickyKeys, stickyMouseButtons;
	bool open = true;

	WindowSettings(ivec2 resolution, const string& windowTitle, bool stickyKeys=true, bool stickyMouseButtons=true);
	WindowSettings(Resolution resolution, const string& windowTitle, bool stickyKeys=true, bool stickyMouseButtons=true);
};


class Window {
	GLFWwindow* window;
	WindowSettings settings;

public:
	explicit Window(WindowSettings settings);
	~Window();

	void destroy() const;
	void initViewport() const;

	bool isOpen() const;
	void renderFramebufferToScreen() const;
};

