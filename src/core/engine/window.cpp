#include "window.hpp"
#include "exceptions.hpp"


WindowSettings::WindowSettings(ivec2 resolution, const string& windowTitle): resolution(resolution), windowTitle(windowTitle) {}

WindowSettings::WindowSettings(Resolution resolution, const string& windowTitle): windowTitle(windowTitle) {
	switch (resolution) {
	case FHD:
		this->resolution = ivec2(1920, 1080);
		break;
	case HD2K:
		this->resolution = ivec2(2560, 1440);
		break;
	case UHD:
		this->resolution = ivec2(3840, 2160);
		break;
	default:
		THROW(SystemError, "Unknown resolution enum value");
	}
}

Window::Window(int width, int height, const char* title) {
	// glfwInit();
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	this->width = width;
	this->height = height;
	this->aspectRatio = (float)width / (float)height;
	this->window = glfwCreateWindow(width, height, title, nullptr, nullptr);
	if (!this->window) {
		// glfwTerminate(); // removed; handled by Renderer
		exit(2136);
	}
	glfwMakeContextCurrent(this->window);
	glfwGetFramebufferSize(window, &width, &height);
}

Window::Window(Resolution resolution, const char* title) {
	// glfwInit();
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	width = predefinedWidth(resolution);
	height = predefinedHeight(resolution);
	aspectRatio = (float)width / height;
	window = glfwCreateWindow(width, height, title, nullptr, nullptr);
	if (!window) {
		// glfwTerminate(); // removed; handled by Renderer
		throw SystemError("GLFW window creation failed", __FILE__, __LINE__);
	}
	glfwMakeContextCurrent(window);
}

Window::~Window() {
	destroy();
}

void Window::destroy() const {
	glfwDestroyWindow(this->window);
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
	glViewport(0, 0, this->width, this->height);
}
