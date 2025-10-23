#pragma once
#include "exceptions.hpp"
#include "GLFW/glfw3.h"

struct WindowSettings {
	string title = "Application Window";
	int width = 1280;
	int height = 720;
};

class Window {
	WindowSettings settings;
	shared_ptr<GLFWwindow> window;
public:
	explicit Window(const WindowSettings &settings);
	const WindowSettings &getSettings() const;
	ivec2 getResolution() const;
	~Window();
	bool isOpen() const;
	void update() const;
};