#pragma once
#include "GLFW/glfw3.h"
#include "exceptions.hpp"


enum Resolution {
	FHD,
	HD2K,
	UHD,
	UNKNOWN
};

struct WindowSettings {
	ivec2 resolution;
	string windowTitle;

	WindowSettings(ivec2 resolution, const string& windowTitle);
	WindowSettings(Resolution resolution, const string& windowTitle);
};


class Window {
public:
	GLFWwindow* window;
	int width;
	int height;
	float aspectRatio;
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

