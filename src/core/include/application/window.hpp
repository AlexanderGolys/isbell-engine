#pragma once
#include "exceptions.hpp"
#include "GLFW/glfw3.h"


enum Resolution {
	FHD,
	HD2K,
	UHD,
	UNKNOWN
};


struct WindowSettings {
	string title;
	int width;
	int height;

	WindowSettings(const string& title, int width, int height);
	WindowSettings(const string& title, Resolution resolution);
};

class Window {
	WindowSettings settings;
	shared_ptr<GLFWwindow> window;

public:
	explicit Window(const WindowSettings& settings);
	const WindowSettings& getSettings() const;
	ivec2 getResolution() const;
	~Window();
	bool isOpen() const;
	void update() const;
};
