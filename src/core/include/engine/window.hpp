#pragma once
#include "exceptions.hpp"


enum class TypicalResolution : uint {
	p480 = 480,
	HD720 = 720,
	FHD = 1080,
	HD2K = 1440,
	UHD = 2160
};

struct Resolution {
	int width;
	int height;

	Resolution(int width, int height);
	explicit Resolution(TypicalResolution r);
};

struct WindowSettings {
	Resolution resolution;
	string windowTitle;
	bool stickyKeys;
	bool stickyMouseButtons;
	bool open = true;

	WindowSettings(Resolution resolution, string_cr windowTitle, bool stickyKeys=true, bool stickyMouseButtons=true);
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

