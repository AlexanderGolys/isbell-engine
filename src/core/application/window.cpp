#include "window.hpp"

Window::Window(const WindowSettings& settings): settings(settings) {
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	window = shared_ptr<GLFWwindow>(
		glfwCreateWindow(settings.width, settings.height, settings.title.c_str(), nullptr, nullptr),
		[](GLFWwindow* ptr) { glfwDestroyWindow(ptr); }
	);
	if (!this->window)
		THROW(SystemError, "GLFW window creation failed");
	glfwMakeContextCurrent(this->window.get());
}

const WindowSettings& Window::getSettings() const {
	return settings;
}
ivec2 Window::getResolution() const {
	return ivec2(settings.width, settings.height);
}
Window::~Window() {
	glfwDestroyWindow(this->window.get());
}

bool Window::isOpen() const {
	return !glfwWindowShouldClose(window.get());
}

void Window::update() const {
	glfwSwapBuffers(this->window.get());
	glfwPollEvents();
}
