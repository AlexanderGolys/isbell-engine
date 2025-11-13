#include "window.hpp"
#include "exceptions.hpp"
#include "glCommand.hpp"


ivec2 ires2(Resolution r) {
	switch (r) {
		case Resolution::FHD: return ivec2(1920, 1080);
		case Resolution::HD2K: return ivec2(2560, 1440);
		case Resolution::UHD: return ivec2(3840, 2160);
	}
	THROW(SystemError, "Unknown resolution enum value");
}

WindowSettings::WindowSettings(ivec2 resolution, const string& windowTitle)
: width(resolution.x), height(resolution.y), windowTitle(windowTitle) {}

WindowSettings::WindowSettings(Resolution resolution, const string& windowTitle) : WindowSettings(ires2(resolution), windowTitle) {}

Window::Window(WindowSettings settings) : settings(settings) {
	window = GLFWCommand::createWindow(settings.width, settings.height, settings.windowTitle.c_str());
}

Window::~Window() {
	destroy();
}

void Window::destroy() const {
	GLFWCommand::destroyWindow(this->window);
	GLFWCommand::terminate();
}

void Window::renderFramebufferToScreen() const {
	glfwSwapBuffers(this->window);
}

void Window::showCursor() const {
	glfwSetInputMode(this->window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
}

void Window::disableCursor() const {
	glfwSetInputMode(this->window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
}

void Window::hideCursorWithinWindow() const {
	glfwSetInputMode(this->window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
}

void Window::stickyKeys(bool sticky) const {
	if (sticky)
		glfwSetInputMode(this->window, GLFW_STICKY_KEYS, GL_TRUE);
	else
		glfwSetInputMode(this->window, GLFW_STICKY_KEYS, GL_FALSE);
}

void Window::stickyMouseButtons(bool sticky) const {
	if (sticky)
		glfwSetInputMode(this->window, GLFW_STICKY_MOUSE_BUTTONS, GL_TRUE);
	else
		glfwSetInputMode(this->window, GLFW_STICKY_MOUSE_BUTTONS, GL_FALSE);
}

void Window::setCallbacks(const GLFWkeyfun* keyCallback, const GLFWcharfun* charCallback, const GLFWmousebuttonfun* mouseButtonCallback, GLFWcursorposfun* cursorPosCallback,
						  GLFWcursorenterfun* cursorEnterCallback, GLFWscrollfun* scrollCallback, GLFWdropfun* dropCallback) const {
	if (keyCallback != nullptr)
		glfwSetKeyCallback(this->window, *keyCallback);
	if (charCallback != nullptr)
		glfwSetCharCallback(this->window, *charCallback);
	if (mouseButtonCallback != nullptr)
		glfwSetMouseButtonCallback(this->window, *mouseButtonCallback);
	if (cursorPosCallback != nullptr)
		glfwSetCursorPosCallback(this->window, *cursorPosCallback);
	if (cursorEnterCallback != nullptr)
		glfwSetCursorEnterCallback(this->window, *cursorEnterCallback);
	if (scrollCallback != nullptr)
		glfwSetScrollCallback(this->window, *scrollCallback);
	if (dropCallback != nullptr)
		glfwSetDropCallback(this->window, *dropCallback);
}

bool Window::isOpen() const {
	return !glfwWindowShouldClose(this->window);
}

void Window::initViewport() const {
	glViewport(0, 0, settings.width, settings.height);
}
