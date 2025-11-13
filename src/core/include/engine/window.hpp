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

	WindowSettings(ivec2 resolution, const string& windowTitle);
	WindowSettings(Resolution resolution, const string& windowTitle);
};


class Window {
	GLFWwindow* window;
	WindowSettings settings;

public:
	explicit Window(WindowSettings settings);
	~Window();
	void destroy() const;

	void initViewport() const;

	void showCursor() const;
	void disableCursor() const;
	void hideCursorWithinWindow() const;
	void stickyKeys(bool sticky) const;
	void stickyMouseButtons(bool sticky) const;
	void setCallbacks(const GLFWkeyfun* keyCallback = nullptr, const GLFWcharfun* charCallback = nullptr, const GLFWmousebuttonfun* mouseButtonCallback = nullptr,
					  GLFWcursorposfun* cursorPosCallback = nullptr, GLFWcursorenterfun* cursorEnterCallback = nullptr, GLFWscrollfun* scrollCallback = nullptr,
					  GLFWdropfun* dropCallback = nullptr) const;
	bool isOpen() const;
	void renderFramebufferToScreen() const;
};

